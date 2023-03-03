
abstract type ITPFuncParams end

struct PolynomialEvalParams{T} <: ITPFuncParams
    c::Vector{T}
    t0::T # could get rid of this, since we're evaluating x:= t-t0, instead of t.
    # c is set up for variable x.
end

function evalITPobjective(x::T, p::PolynomialEvalParams{T})::T where T
    #return evaltaylor(p.c, x, p.t0)
    return evaltaylor(p.c, x, zero(T))
end


# https://doi-org.proxy.library.carleton.ca/10.1145/3423597
# assumes f(a) < 0 < f(b).
# f = xx->evalITPobjective(x, params)
function runITP(
    lb::T,
    ub::T,
    params::ITPFuncParams,
    config::ITPConfig{T},
    ) where T

    # check if we have a valid problem to solve. See equation 2.
    f_lb = evalITPobjective(lb, params)
    f_ub = evalITPobjective(ub, params)
    if !(f_lb*f_ub < 0) || lb > ub
        # cannot do root finding on this interval.
        return convert(T, NaN), false
    end

    # set up constants.
    f_tol = config.f_tol
    x_tol = config.x_tol
    n0 = config.n0

    k1 = config.k1
    @assert k1 > zero(T)
    
    k2 = config.k2
    @assert one(T) <= k2 < one(T) + MathConstants.golden
    #ϕ = 0.5*(1+sqrt(5)) # the golden ratio.
    #ϕ = MathConstants.golden
    #k2 = 1 + k2_percent*ϕ

    ϵ = x_tol/2
    n_bin = ceil(Int, log2((ub-lb)/x_tol)) # Theorem 1.1.
    n_max = n_bin + n0

    # # bracketing algorithm.
    k = 0
    a = lb
    b = ub
    f_a = f_lb
    f_b = f_ub
    C = 2^(n_max+1)

    while abs(b-a) > x_tol

        # interpolation.
        x_RF = (a*f_b - b*f_a)/(f_b - f_a) # equation 5.
        
        # truncation, x_t.
        x_bin = (a+b)/2 # equation 4

        σ = sign(x_bin - x_RF)
        δ = k1*abs(b-a)^k2
    
        x_t = x_bin
        if δ <= abs(x_bin - x_RF)
            x_t = x_RF + σ*δ
        end

        # minimax radius and interval.
        #r_k = ϵ*2^(n_bin-k) - (b-a)/2
        C = C/2
        r_k = ϵ*C - (b-a)/2
        
        #interval_k = [x_bin - r_k, x_bin + r_k]

        # projection onto interval_k.
        x_ITP = x_t
        if abs(x_t -x_bin) > r_k
            x_ITP = x_bin - σ*r_k
        end

        # choose a new position in the interval, and see if it is numerically close to being a root.
        w = x_ITP
        #w = (a + b)/2 # uncomment to overide ITP with bisection.
        
        f_w = evalITPobjective(w, params) #f(w)
        
        if abs(f_w) < f_tol
            return clamp(w, lb, ub), true
        end

        # Set up for the next iteration: bracket steps.
        if f_w*f_a > 0
            a = w
            f_a = f_w
            # b = b
        elseif f_w*f_b > 0
            #a = a
            b = w
            f_b = f_w
        else
            a = w
            b = w
        end

        k += 1 # first iteration is defined as k == 0. See the paragraph after equation 16.
    end

    return clamp((a + b)/2, lb, ub), false
end

# front end for intersection polynomials in `cs`
# cs is from BudanIntersectionBuffers' cs field.
function runITP(
    cs::Vector{Vector{T}}, # [constraints][order]
    t0::T,
    h::T,
    config::ITPConfig{T},
    ) where T <: AbstractFloat

    # set up.
    a = zero(T)
    b = h
    
    min_t = h
    found_root = false


    # solve roots.
    for m in eachindex(cs)
        t, status_flag = runITP(a, b, PolynomialEvalParams(cs[m], t0), config)
        #@show t, status_flag

        # the polynomial is usually very ill-conditioned, so f_tol is unlikely to be reached.
        # check explicitly for root in the calling routine.
        if t < min_t #&& status_flag
            min_t = t
            found_root = true
        end
    end

    return min_t, found_root
end


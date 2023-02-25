
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

struct ITPConfig{T}
    f_tol::T
    x_tol::T
    k1::T
    k2_percent::T
end

function ITPConfig(
    ::Type{T};
    f_tol::T = convert(T, 1e-8),
    x_tol::T = convert(T, 1e-15),
    k1::T = one(T),
    k2_percent::T = convert(T, 0.5),
    )::ITPConfig{T} where T

    return ITPConfig(f_tol, x_tol, k1, k2_percent)
end


# https://doi-org.proxy.library.carleton.ca/10.1145/3423597
function runITP(
    a::T,
    b::T,
    params::ITPFuncParams,
    config::ITPConfig{T},
    ) where T

    # parse.
    f_tol = config.f_tol
    x_tol = config.x_tol
    k1 = config.k1
    k2_percent = config.k2_percent

    @assert !(k1 < zero(T))
    @assert zero(T) <= k2_percent < one(T)
    ϵ = x_tol/2

    #ϕ = 0.5*(1+sqrt(5)) # the golden ratio.
    ϕ = MathConstants.golden
    k2 = 1 + k2_percent*ϕ

    x_bin = (a+b)/2 # equation 4
    n_bin = ceil(Int, log2((b-a)/x_tol)) # Theorem 1.1.

    f_a = evalITPobjective(a, params) # f(a)
    f_b = evalITPobjective(b, params) #f(b)
    x_RF = (b*f_a-a*f_b)/(f_a-f_b) # equation 5.
    
    σ = sign(x_bin - x_RF)
    δ = k1*abs(b-a)^k2

    x_t = x_bin
    if δ <= abs(x_bin - x_RF)
        x_t = x_RF + σ*δ
    end

    # # bracketing algorithm.
    f_lb = f_a
    f_ub = f_b

    a_k = a
    b_k = b
    k = 0
    C = 2^(n_bin+1)
    while abs(b_k - a_k) > x_tol

        # minimax radius and interval.
        #r_k = ϵ*2^(n_bin-k) - (b_k-a_k)/2
        r_k = ϵ*C/2 - (b_k-a_k)/2
        
        #interval_k = [x_bin - r_k, x_bin + r_k]

        # projection onto interval_k.
        x_ITP = x_t
        if abs(x_t -x_bin) > r_k
            x_ITP = x_bin - σ*r_k
        end

        # bracket steps:
        #w = x_ITP
        w = (a_k + b_k)/2 # overide with bisection for now. ITP isn't working.

        #@show w
        f_w = evalITPobjective(w, params) #f(w)

        if f_w*f_lb > 0
            a_k = w
            # b_k = b_k
        elseif f_w*f_ub > 0
            #a_k = a_k
            b_k = w
        end

        if abs(f_w) < f_tol
            return w, true
        end

        k += 1 # first iteration is defined as k == 0. See the paragraph after equation 16.
    end

    return (a_k + b_k)/2, false
end

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
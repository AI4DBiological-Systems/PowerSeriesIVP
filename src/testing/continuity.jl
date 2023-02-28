#################

function continuitycheck(
    sol::PiecewiseTaylorPolynomial{T};
    atol = eps(T)*10,
    ) where T
    
    N = getNpieces(sol)
    N_vars = getNvars(sol)
    for n in Iterators.take(eachindex(sol.coefficients), N-1)


        h = sol.steps[n]
        t0 = sol.expansion_points[n]
        t_next = t0 + h # the expansion point of the next solution piece, n+1.
        
        c_x = sol.coefficients[n].x
        c_u = sol.coefficients[n].u

        c_x_next = sol.coefficients[n+1].x
        c_u_next = sol.coefficients[n+1].u

        for d = 1:N_vars
            
            # # separate.
            pass_flag = initialconditioncheck(c_x[d], c_x_next[d], h, t0)
            pass_flag = pass_flag & initialconditioncheck(c_u[d], c_u_next[d], h, t0)
            
            # ignore parallel transport vectors for now.
            if !pass_flag
                return false
            end

            # # joint.
            if abs(c_x_next[d][begin+1] - c_u_next[d][begin]) > atol
                return false
            end
        end
    end

    return true
end

function initialconditioncheck(c_d::Vector{T}, c_d_next::Vector{T}, h, t0; atol = eps(T)*10 ) where T
    
    sol_t = evaltaylor(c_d, t0+h, t0)
    sol_next_t = c_d_next[begin] # 0-th order coefficient of x_next.

    if abs(sol_t - sol_next_t) > atol
        return false
    end

    return true
end

"""
differentiatepolynomial!(out::Vector{T}, c::Vector{T})::Nothing where T

Do not pass the same array as `out` and `c`.
Test code:
L = 9
c = randn(L+1)
dc = similar(c)
t0 = rand()

f = tt->PowerSeriesIVP.evaltaylor(c, tt, t0)
f_AD = tt->PowerSeriesIVP.evaltaylorAD(c, tt, t0)
df_AD = tt->ForwardDiff.derivative(f_AD, tt)
f_DR = tt->PowerSeriesIVP.evaltaylordirect(c, tt, t0)

PowerSeriesIVP.differentiatepolynomial!(dc, c)
df_AN = tt->PowerSeriesIVP.evaltaylor(dc, tt, t0)

t = t0 + 0.1
@show df_AN(t)
@show df_AD(t)
@show f(t)
@show f_DR(t)
@show f_AD(t)
"""
function differentiatepolynomial!(out::Vector{T}, c::Vector{T})::Nothing where T
    resize!(out, length(c)-1)

    for i in eachindex(out)
        out[i] = i*c[begin+i]
    end

    return nothing
end

function getderivativediscrepancies(
    sol::PiecewiseTaylorPolynomial{T},
    M::Integer, # we analyze up to this order, for x.
    ) where T
    
    N_vars = getNvars(sol)
    N = getNpieces(sol)
    N_boundaries = N-1 # excludeing the end points.

    # intermediate buffer.
    dx_coefficients = Vector{T}(undef, 0)
    du_coefficients = Vector{T}(undef, 0)
    
    # output.
    endpoint_err = Matrix{T}(undef, N_vars, N_boundaries)
    fill!(endpoint_err, convert(T, NaN))

    endpoint_err_x = collect( Matrix{T}(undef, N_vars, N_boundaries) for _ = 0:M )
    for i in eachindex(endpoint_err_x)
        fill!(endpoint_err_x[i], convert(T, NaN))
    end

    xu0_err = Matrix{T}(undef, N_vars, N_boundaries)
    fill!(xu0_err, convert(T, NaN))

    xuh_err = Matrix{T}(undef, N_vars, N_boundaries)
    fill!(xuh_err, convert(T, NaN))

    endpoint_err_u = collect( Matrix{T}(undef, N_vars, N_boundaries) for _ = 0:M )
    for i in eachindex(endpoint_err_u)
        fill!(endpoint_err_u[i], convert(T, NaN))
    end

    interval_err = zeros(T, N_vars, N_boundaries)

    orders = Vector{Int}(undef, N_boundaries)

    for n in Iterators.take(eachindex(sol.coefficients), N_boundaries)


        h = sol.steps[n]
        t0 = sol.expansion_points[n]
        t_next = t0 + h # the expansion point of the next solution piece, n+1.
        t_range = LinRange(t0, t_next, 100)

        c_x = sol.coefficients[n].x
        c_u = sol.coefficients[n].u

        c_x_next = sol.coefficients[n+1].x
        c_u_next = sol.coefficients[n+1].u

        for d = 1:N_vars
            
            # this should be essentially eps(), since the PSM was based on this relation.
            xu0_err[d,n] = abs(c_u[d][begin] - c_x[d][begin+1])

            # # taylor solution derivative vs. u0.
            differentiatepolynomial!(dx_coefficients, c_x[d])
            #@show c_u_next[d][begin]            
            endpoint_err[d,n] = abs(evaltaylor(dx_coefficients, t_next, t0) - c_u_next[d][begin])

            # shold be the same as endpoint_err[d,n]
            xuh_err[d,n] = abs(evaltaylor(dx_coefficients, t_next, t0) - evaltaylor(c_u[d], t_next, t0))

            # # interval error.
            for t in t_range
                # l-1 norm
                interval_err[d,n] += abs(evaltaylor(dx_coefficients, t, t0) - evaltaylor(c_u[d], t, t0))

                # l-2 norm.
                #interval_err[d,n] += (evaltaylor(dx_coefficients, t, t0) - evaltaylor(c_u[d], t, t0))^2
            end

            # # errors for the derivatives of x. Overwrites dx_coefficients
            # the 0-th order case.
            endpoint_err_x[begin][d,n] = abs(evaltaylor(c_x[d], t_next, t0) - c_x_next[d][begin])

            dx_coefficients = copy(c_x[d])
            for m = 1:M
                differentiatepolynomial!(dx_coefficients, copy(dx_coefficients))

                endpoint_err_x[begin+m][d,n] = abs(evaltaylor(dx_coefficients, t_next, t0) - c_x_next[d][begin+m])
            end

            # # similarly for u.
            # the 0-th order case.
            endpoint_err_u[begin][d,n] = abs(evaltaylor(c_u[d], t_next, t0) - c_u_next[d][begin])

            du_coefficients = copy(c_u[d])
            for m = 1:M
                differentiatepolynomial!(du_coefficients, copy(du_coefficients))

                endpoint_err_u[begin+m][d,n] = abs(evaltaylor(du_coefficients, t_next, t0) - c_u_next[d][begin+m])
            end

        end

        orders[n] = getorder(sol.coefficients[n])
    end

    return xu0_err, xuh_err, endpoint_err, endpoint_err_x, endpoint_err_u, interval_err, orders
end

function evaltaylorAD(c, x, a)
    τ = x-a
    return sum( c[i]*τ^(i-1) for i in eachindex(c) )
end

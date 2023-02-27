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


function differentiatepolynomial!(out::Vector{T}, c::Vector{T}) where T
    resize!(out, length(c)-1)

    for i in eachindex(out)
        out[i] = c[begin+i]
    end

    return nothing
end

function getderivativediscrepancies(
    sol::PiecewiseTaylorPolynomial{T},
    ) where T
    
    
    N_vars = getNvars(sol)
    N = getNpieces(sol)
    N_boundaries = N-1 # excludeing the end points.

    # intermediate buffer.
    dx_coefficients = Vector{T}(undef, 0)
    
    # output.
    discrepancies = Matrix{T}(undef, N_vars, N_boundaries)
    fill!(discrepancies, convert(T, NaN))

    orders = Vector{Int}(undef, N_boundaries)

    for n in Iterators.take(eachindex(sol.coefficients), N_boundaries)


        h = sol.steps[n]
        t0 = sol.expansion_points[n]
        t_next = t0 + h # the expansion point of the next solution piece, n+1.
        
        c_x = sol.coefficients[n].x
        #c_u = sol.coefficients[n].u

        #c_x_next = sol.coefficients[n+1].x
        c_u_next = sol.coefficients[n+1].u

        for d = 1:N_vars
 
            # # taylor solution derivative vs. u0.
            differentiatepolynomial!(dx_coefficients, c_x[d])
            #@show c_u_next[d][begin]            
            discrepancies[d,n] = evaltaylor(dx_coefficients, t_next, t0) - c_u_next[d][begin]            
        end

        orders[n] = getorder(sol.coefficients[n])
    end

    return discrepancies, orders
end
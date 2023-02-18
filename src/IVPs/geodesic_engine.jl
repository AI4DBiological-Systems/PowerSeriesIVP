# main routine for computing the power series method iterations for the generic geodesic IVP problem.


function getfirstorder!(p::GeodesicIVPBuffer, _::DisableParallelTransport)

    # variables
    initializeorder!(p.x, p.x0)
    initializeorder!(p.u, p.u0)

    # du/dt
    initializeorder!(p.θ, p.x.c, p.u.c)

    # variables
    increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
    increaseorder!(p.u, p.θ.c)
    
    return nothing
end


function increaseorder!(p::GeodesicIVPBuffer, _::DisableParallelTransport)

    # du/dt
    increaseorder!(p.θ, p.x.c, p.u.c)

    # variables
    increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
    increaseorder!(p.u, p.θ.c)

    return nothing
end


###### adaptive step

# eq:error_estimate
function computeerror(c::Vector{T}, h::T, N_analysis_terms::Integer)::T where T
    L_analysis = length(c) - 1
    @assert L_analysis > N_analysis_terms > 0

    L = L_analysis - N_analysis_terms

    return abs(sum( c[begin+n]*h^(n-1) for n = (L+1):L_analysis ))
end

# function computesumerror(cs::Vector{Vector{T}}, h::T, N_analysis_terms::Integer)::T where T
    
#     total_err = zero(T)
#     for d in eachindex(cs)
        
#         total_err += computeerror(cs[d], h, N_analysis_terms)
#     end

#     return total_err
# end

function computemaxerror!(
    error_buffer::Vector{T},
    cs::Vector{Vector{T}},
    h::T,
    N_analysis_terms::Integer,
    ) where T
    
    resize!(error_buffer, length(cs))

    for d in eachindex(cs)
        error_buffer[d] = computeerror(cs[d], h, N_analysis_terms)
    end

    max_error, ind = findmax(error_buffer)

    return max_error, ind
end

# eq:choose_ODE_step_size. An h_zero_error of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsize(ϵ::T, c::Vector{T}; h_zero_error::T = convert(T,Inf))::T where T
    L = length(c) - 1
    h = stepsizeformula(ϵ, c[end], L)
    if isfinite(h)
        return h
    end

    return h_zero_error
end

function stepsizeformula(ϵ::T, x::T, L::Integer)::T where T
    return (abs(ϵ/(2*x)))^(1/(L-1))
end

# take care of constraints in the piece-wise solution part, not this inner routine for an individual Taylor solution.
# final order is N_analysis_terms + L_test.
# Some functions at particular exapnsion points have Taylor series with vanishing coefficients. Example: https://www.wolframalpha.com/input?i=taylor+expansion+of+%281%2Bx%5E4%29%2F%282%2Bx%5E4%29&key=33rzp
# - N_analysis_terms > 1 attempts to mitigate this issue.
# - Since N_analysis_terms cannot be infinite, use h_zero_error to guard the case when zero estimated error is detected. h_zero_error = Inf means no guard, and that zero error means the solution polynomial is exactly the solution of the DE.
function computetaylorsolution!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT;
    ϵ::T = convert(T, 1e-6),
    h_initial::T = convert(T, NaN),
    L_test_max::Integer = 10,
    r_order::T = convert(T, 0.3),
    h_zero_error::T = convert(T, Inf),
    N_analysis_terms = 2, # make larger than 2 if the solution has many consecutive vanishing Taylor coefficients for certain orders.
    ) where {PT,T}

    ### set up
    
    p = prob
    error_threshold = ϵ/2

    err_record = ones(T, L_test_max) # first entry is for 1st-order, so on.
    fill!(err_record, Inf)

    # intermediate buffers.
    N_vars = length(prob.x0)
    error_across_variables = Vector{T}(undef, N_vars)

    ### get the highest compute order N_analysis_terms + 1. The +1 is since we want a minimum of a 1-st order (line) solution.
    getfirstorder!(prob, pt_trait)
    
    ### get to a high enough order so that we can start computing the error.
    for _ = 1:N_analysis_terms
        increaseorder!(prob, pt_trait)
    end

    # are we using the default initial h?
    h = h_initial
    if !isfinite(h)
        # use the following heurestic in this if-block for initial h:

        # use the first (arb. chosen) variable to analyze how far to step.
        expansion_factor = 100.0 # hard code for now.
        h = choosestepsize(ϵ, p.x.c[begin]; h_zero_error = h_zero_error)*expansion_factor
    end

    # error for the l-th order.
    #err_record[begin] = computeerror(p.x.c, h, N_analysis_terms)

    error_val, error_var_ind = computemaxerror!(error_across_variables, p.x.c, h, N_analysis_terms)
    err_record[begin] = error_val

    if error_val < error_threshold

        # error is already tolerable. no need for adaptation. find step and exit.
        return choosestepsize(ϵ, p.x.c[error_var_ind]; h_zero_error = h_zero_error)
    end
    

    ### start adaption strategy.
    for l = (1+N_analysis_terms):L_test_max # l is L_test
        # increase order.
        increaseorder!(prob, pt_trait)

        #err_record[l] = computeerror(p.x.c, h, N_analysis_terms) # error for 1st-order solution.
        error_val, error_var_ind = computemaxerror!(error_across_variables, p.x.c, h, N_analysis_terms)
        
        # eq:sufficient_error_reduction_order, eq:error_threshold_condition
        shrink_change = err_record[l]/err_record[l-1]

        if error_val < error_threshold || shrink_change < r_order
        
            # Stop increasing the order. Use highest computed order coefficient to get step size.
            return choosestepsize(ϵ, p.x.c[error_var_ind]; h_zero_error = h_zero_error)
        end

        err_record[l] = error_val # book keep.
    end

    # Reached max order. Use highest computed order coefficient to get step size.
    return choosestepsize(ϵ, p.x.c[error_var_ind]; h_zero_error = h_zero_error)
end

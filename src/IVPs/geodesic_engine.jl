# main routine for computing the power series method iterations for the generic geodesic IVP problem.

function getorder(p::GeodesicIVPBuffer)::Int
    return length(p.x.c[begin])-1
end

function getfirstorder!(p::GeodesicIVPBuffer, _::DisableParallelTransport)

    # variables
    initializeorder!(p.x, p.x0)
    initializeorder!(p.u, p.u0)

    # du/dt
    initializeorder!(p.θ, p.x.c, p.u.c)

    # variables: first update.
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

# only decrease x and u. Does not decrease intermediates variables and their buffers, such as  θ.
function decreaseorder!(p::GeodesicIVPBuffer, _::DisableParallelTransport, min_order::Integer)

    # if length(p.x.c) < 1
    #     println("Warning, decreaseorder!() received an polynomial less than order 1. Did not decrease order further.")
    #     return nothing
    # end

    # variables
    decreaseorder!(p.x, 1, min_order)
    decreaseorder!(p.u, 1, min_order)

    return nothing
end

function getfirstorder!(p::GeodesicIVPBuffer, _::EnableParallelTransport)

    pt = p.parallel_transport

    # set up the geodesic variables and related quantities.
    getfirstorder!(p, DisableParallelTransport())

    # set up the transport vector variables and quantities.
    for m in eachindex(pt)
        
        # variables
        initializeorder!(pt[m].v, pt[m].v0)

        # du/dt
        initializeorder!(pt[m].ζ, p.θ, pt[m].v.c)

        # variables: first update.
        increaseorder!(pt[m].v, pt[m].ζ.c)

        # increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
        # increaseorder!(p.u, p.θ.c)
    end

    return nothing
end


function increaseorder!(p::GeodesicIVPBuffer, _::EnableParallelTransport)

    pt = p.parallel_transport

    # set up the geodesic variables and related quantities.
    increaseorder!(p, DisableParallelTransport())

    # update the transport vector variables and quantities.
    for m in eachindex(pt)
    
        # dv/dt
        increaseorder!(pt[m].ζ, p.θ, pt[m].v.c)

        # variables: first update.
        increaseorder!(pt[m].v, pt[m].ζ.c)

        # increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
        # increaseorder!(p.u, p.θ.c)
    end
    
    return nothing
end

# only decrease x and u and {pt[m].v}_m. Does not decrease intermediates variables and their buffers, such as  θ and ζ.
function decreaseorder!(p::GeodesicIVPBuffer, _::EnableParallelTransport, min_order::Integer)

    # if length(p.x.c) < 1
    #     println("Warning, decreaseorder!() received an polynomial less than order 1. Did not decrease order further.")
    #     return nothing
    # end

    pt = p.parallel_transport

    # decrease the position p.x and geodesic velocity/vector field p.u.
    decreaseorder!(p, DisableParallelTransport(), min_order)

    # update the transport vector variables and quantities.
    for m in eachindex(pt)
        decreaseorder!(pt[m].v, 1, min_order)
    end

    return nothing
end

###### adaptive step

# eq:error_estimate backup.
# function computeerror(
#     c::Vector{T},
#     h::T,
#     N_analysis_terms::Integer,
#     )::T where T

#     L_analysis = length(c) - 1
#     #@show L_analysis, N_analysis_terms, 0
#     @assert L_analysis > N_analysis_terms > 0

#     L = L_analysis - N_analysis_terms

#     return abs(sum( c[begin+n]*h^(n-1) for n = (L+1):L_analysis ))
# end

# This is E(M,h) in my notes.
function computeerror(
    c::Vector{T},
    h::T,
    L::Integer,
    N_analysis_terms::Integer,
    )::T where T

    order = length(c) - 1
    @assert order > N_analysis_terms > 0

    L = order - N_analysis_terms

    return abs(sum( c[begin+n]*h^(n-1) for n = (L+1):order ))
end

# based on sum error formula in my notes.
function computeerrorratio(
    x::Vector{Vector{T}},
    ϵ::T,
    N_analysis_terms::Integer;
    h_zero_error = Inf,
    step_reduction_factor = 2,
    )::T where T
    
    order = length(x[begin])
    @assert order > N_analysis_terms > 0

    M = order - N_analysis_terms
    
    a = choosestepsize(
        ϵ,
        x;
        order = M,
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )

    b = choosestepsize(
        ϵ,
        x;
        order = M-1,
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )

    numerator = zero(T)
    denominator = zero(T)
    for d in eachindex(x)
        
        x_d = x[d]

        numerator += abs(sum( x_d[end-n+1]*a^(n-1) for n = 1:N_analysis_terms ))
        denominator += abs(sum( x_d[end-1-n]*b^(n-1) for n = 1:N_analysis_terms ))
    end
    error_ratio = (numerator*a^M)/(denominator*b^(M-1))

    return error_ratio
end

# function computemaxerror(
#     cs::Vector{Vector{T}},
#     h::T,
#     N_analysis_terms::Integer,
#     ) where T
    
#     max_error = convert(T, -Inf)
#     ind = 1 # at least a default valid index, in case all computed errors are -Inf.

#     for d in eachindex(cs)
#         current_error = computeerror(cs[d], h, N_analysis_terms)

#         if current_error > max_error
#             ind = d
#             max_error = current_error
#         end
#     end

#     return max_error, ind
# end

# function computemaxerror!(
#     error_buffer::Vector{T},
#     cs::Vector{Vector{T}},
#     h::T,
#     N_analysis_terms::Integer,
#     ) where T
    
#     resize!(error_buffer, length(cs))

#     for d in eachindex(cs)
#         error_buffer[d] = computeerror(cs[d], h, N_analysis_terms)
#     end

#     max_error, ind = findmax(error_buffer)

#     return max_error, ind
# end

# eq:choose_ODE_step_size. An h_zero_error of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsize(
    ϵ::T,
    c::Vector{T};
    h_zero_error = Inf,
    step_reduction_factor = 2,
    order = length(c) -1,
    )::T where T

    L = order
    total_order = length(c) - 1
    @assert 0 < L <= total_order
    
    h = stepsizeformula(ϵ, c[begin+L], L)
    if isfinite(h)
        # further shorten step size as a heuristic strategy to reduce error.
        return convert(T, h/step_reduction_factor)
    end

    return convert(T, h_zero_error/step_reduction_factor)
end

# see http://www.phys.uri.edu/nigh/NumRec/bookfpdf/f16-2.pdf for a constrast of the method from PSM 2019 with mainstream RK adaptive step size selection algorithms.
# we hardcoded such that only one analysis term for the step size determination.
function stepsizeformula(ϵ::T, x::T, L::Integer)::T where T
    return (abs(ϵ/(2*x)))^(1/(L-1))
end

function choosestepsize(
    ϵ::T,
    c_x::Vector{Vector{T}};
    order = length(c_x[begin])-1,
    h_zero_error = Inf,
    step_reduction_factor = 2,
    )::T where T

    min_h = convert(T, Inf)
    for d in eachindex(c_x)
        h = choosestepsize(
            ϵ,
            c_x[d];
            h_zero_error = h_zero_error,
            step_reduction_factor = step_reduction_factor,
            order = order,
        )
        min_h = min(min_h, h)
    end

    return min_h
end

# no constraints.
function computetaylorsolution!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT,
    #h_test::T,
    config::FixedOrderConfig{T},
    _::NoConstraints,
    ) where {T, PT}

    ### get to a high enough order so that we can start computing the error.
    #@show length(prob.x.c[1]), length(prob.u.c[1])
    # from the calling function, firstorder!() got prob.x and prob.u up to order 1 already. Start at order 2.
    for _ = 2:config.L
        increaseorder!(prob, pt_trait)
    end
    
    #error_val, error_var_ind = computemaxerror(prob.x.c, h_test, 1)
    # return choosestepsize(
    #     config.ϵ,
    #     prob.x.c[error_var_ind];
    #     h_zero_error = config.h_zero_error,
    #     step_reduction_factor = config.step_reduction_factor,
    # )
    h = choosestepsize(
        config.ϵ,
        prob.x.c;
        h_zero_error = config.h_zero_error,
        step_reduction_factor = config.step_reduction_factor,
    )

    return h, :continue_simulation
end

function computetaylorsolution!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT,
    config::AdaptOrderConfig{T},
    C::ConstraintType,
    ) where {T, PT}

    ### set up

    ϵ = config.ϵ
    L_min = config.L_min
    L_max = config.L_max
    r_order = config.r_order
    h_zero_error = config.h_zero_error
    N_analysis_terms = config.N_analysis_terms
    step_reduction_factor = config.step_reduction_factor
    p = prob

    #explicit_roots_buffer = constraints_container.explicit_roots_buffer
    #constraints = constraints_container.constraints

    ### get to a high enough order so that we can start computing the error.
    # from the calling function, firstorder!() got prob.x and prob.u up to order 1 already.
    # start from 2, but make sure we have at least an extra order number to do order-adaption's error analysis. Look into this later.
    L_start = getorder(prob) + 1
    for _ = L_start:L_min
        increaseorder!(prob, pt_trait)
    end

    # start to decide if we should exit.
    h = choosestepsize(
        ϵ,
        p.x.c;
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )
    
    h_new, constraint_ind = refinestep!(
        C,
        h,
        p.x.c,
    )

    if constraint_ind > 0
        # found intersection.
        return h_new, :stop_simulation
    end
    

    ### start adaption strategy.
    h_prev = h

    r = r_order + 1 # initialize to a value so we satisfy the loop condition the first time.
    while getorder(prob) < L_max && r > r_order

        # increase order.
        increaseorder!(prob, pt_trait)

        h = choosestepsize(
            ϵ,
            p.x.c;
            h_zero_error = h_zero_error,
            step_reduction_factor = step_reduction_factor,
        )

        h_new, constraint_ind = refinestepnumerical!(
            C, h, h_prev, p.x.c, 
        )

        #valid_h_new = (0 < h_new < h)
        valid_h_new = isfinite(h_new)

        # we avoid returning inside the if-else if we have no roots in [0, h_new] at the current order.
        if valid_h_new
            # if constraint_ind == 0
            #     # no roots on [0, h_new)

            #     return h_new, :continue_simulation
            # else
            #     # found intersection.

            #     return h_new, :stop_simulation
            # end

            if constraint_ind == 1
                # found intersection.
                return h_new, :stop_simulation
            end
        
        else 
            # the inconclusive case and unexpected error case, like root solve failed.
            # no intersection so far, up to and including the previous order.

            if getorder(prob) > L_min

                decreaseorder!(prob, pt_trait, L_min)
            end
            # don't need to decrease if we're already at order: L_min.

            return h_prev, :continue_simulation
        end

        # estimate the decrease in error by our increase in order.
        # see current_notes under `root find` folder.
        r = computeerrorratio(
            p.x.c,
            ϵ,
            N_analysis_terms;
            h_zero_error = h_zero_error,
            step_reduction_factor = step_reduction_factor,
        )

        # if r < r_order
        #     return h_new
        # end

        h_prev = h_new
    end

    return h_new, :continue_simulation
end

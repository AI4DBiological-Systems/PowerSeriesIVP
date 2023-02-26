# main routine for computing the power series method iterations for the generic geodesic IVP problem.


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
    cs::Vector{Vector{T}},
    N_analysis_terms::Integer;
    h_zero_error = Inf,
    step_reduction_factor = 2,
    )::T where T
    
    order = length(cs[begin])
    @assert order > N_analysis_terms > 0

    M = order - N_analysis_terms
    
    a = choosestepsize(
        ϵ,
        cs;
        order = M,
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )

    b = choosestepsize(
        ϵ,
        cs;
        order = M-1,
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )

    numerator = zero(T)
    denominator = zero(T)
    for d in eachindex(cs)
        
        c = cs[d]

        numerator += abs(sum( c[begin+M+n]*a^(n-1) for n = 1:N_analysis_terms ))
        denominator += abs(sum( c[begin+M-1+n]*b^(n-1) for n = 1:N_analysis_terms ))
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

function applyadaptionstrategy!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT,
    h_test::T,
    config::FixedOrderConfig{T},
    ) where {T, PT}

    ### get to a high enough order so that we can start computing the error.
    #@show length(prob.x.c[1]), length(prob.u.c[1])
    # from the calling function, firstorder!() got prob.x and prob.u up to order 1 already. Start at order 2.
    for _ = 2:config.L
        increaseorder!(prob, pt_trait)
    end
    #@show length(prob.x.c[1])
    #println()
    
    error_val, error_var_ind = computemaxerror(prob.x.c, h_test, 1)
    
    return choosestepsize(
        config.ϵ,
        prob.x.c[error_var_ind];
        h_zero_error = config.h_zero_error,
        step_reduction_factor = config.step_reduction_factor)
end


function applyadaptionstrategy!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT,
    h_initial::T,
    config::AdaptOrderConfig{T},
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

    error_threshold = ϵ/2

    err_record = ones(T, L_test_max) # first entry is for 1st-order, so on.
    fill!(err_record, Inf)

    ### get to a high enough order so that we can start computing the error.
    # from the calling function, firstorder!() got prob.x and prob.u up to order 1 already.
    # start from 2, but make sure we have at least an extra order number to do order-adaption's error analysis. Look into this later.
    for _ = 2:L_max
        increaseorder!(prob, pt_trait)
    end

    # start to decide if we should exit.
    h = choosestepsize(
        ϵ,
        p.x.c;
        h_zero_error = h_zero_error,
        step_reduction_factor = step_reduction_factor,
    )

    if checkroots()
        return h
    end
    

    ### start adaption strategy.
    h_prev = h

    for l = (1+N_analysis_terms):L_test_max # l is L_test
        # increase order.
        increaseorder!(prob, pt_trait)

        h = choosestepsize(
            ϵ,
            p.x.c;
            h_zero_error = h_zero_error,
            step_reduction_factor = step_reduction_factor,
        )

        x_right, status = checkroots()

        if status == :one_root
            
            h_new = runITP()
            return h_new
        
        elseif status == :no_root

            # do nothing.
        else 
            # the inconclusive case and unexpected error case, like root solve failed.

            decreaseorder!()

            return h_prev
        end

        # I am here. see current_notes under `root find` folder.
        r = computeerrorratio(
            cs,
            N_analysis_terms;
            h_zero_error = h_zero_error,
            step_reduction_factor = step_reduction_factor,
        )

        if r < r_order
            return h
        end

        h_prev = h
    end

    return h
end

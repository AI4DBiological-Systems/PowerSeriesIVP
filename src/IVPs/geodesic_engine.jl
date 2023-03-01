# main routine for computing the power series method iterations for the generic geodesic IVP problem.

function getorder(p::GeodesicIVPBuffer)::Int
    return length(p.x.c[begin])-1
end

function getfirstorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    )

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


function increaseorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    )

    # du/dt
    increaseorder!(p.θ, p.x.c, p.u.c)

    # variables
    increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
    increaseorder!(p.u, p.θ.c)

    return nothing
end

# only decrease x and u. Does not decrease intermediates variables and their buffers, such as  θ.
function decreaseorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    min_order::Integer,
    )

    # if length(p.x.c) < 1
    #     println("Warning, decreaseorder!() received an polynomial less than order 1. Did not decrease order further.")
    #     return nothing
    # end

    # variables
    decreaseorder!(p.x, 1, min_order)
    decreaseorder!(p.u, 1, min_order)

    return nothing
end

function getfirstorder!(
    p::GeodesicIVPBuffer,
    ::EnableParallelTransport,
    )

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


function increaseorder!(
    p::GeodesicIVPBuffer,
    ::EnableParallelTransport,
    )

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
function decreaseorder!(
    p::GeodesicIVPBuffer, 
    ::EnableParallelTransport,
    min_order::Integer,
    )

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


# no constraints.
function computetaylorsolution!(
    prob::GeodesicIVPBuffer,
    pt_trait::PT,
    #h_test::T,
    config::FixedOrderConfig{T},
    ::NoConstraints,
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
    #     h_max = config.h_max,
    #     step_reduction_factor = config.step_reduction_factor,
    # )
    h = choosestepsize(
        VelocityContinuityStep(),
        config.ϵ,
        prob;
        #GuentherWolfStep();
        h_max = config.h_max,
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
    h_max = config.h_max
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
        VelocityContinuityStep(),
        ϵ,
        p;
        #GuentherWolfStep();
        h_max = h_max,
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
            VelocityContinuityStep(),
            ϵ,
            p,
            #GuentherWolfStep();
            h_max = h_max,
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
            p,
            ϵ,
            N_analysis_terms;
            h_max = h_max,
            step_reduction_factor = step_reduction_factor,
        )

        # if r < r_order
        #     return h_new
        # end

        h_prev = h_new
    end

    return h_new, :continue_simulation
end

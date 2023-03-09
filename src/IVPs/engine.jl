
function solveIVP(
    ::Type{ST},
    prob_params::IVPStatement,
    pt_trait::IVPVariationTrait,
    t_start::T,
    t_fin::T,
    config::IVPConfig,
    constraints_info::ConstraintsTrait,
    )::Tuple{ST,Symbol} where {T, ST}

    sol = PiecewiseTaylorPolynomial(T, getsolpiecetype(prob_params))
    eval_buffer = allocatevariablecontainer(prob_params)

    exit_flag = solveIVP!(
        sol,
        eval_buffer,
        prob_params,
        pt_trait,
        t_start,
        t_fin,
        config,
        constraints_info,
    )

    return sol, exit_flag
end

function solveIVP!(
    sol::PiecewiseTaylorPolynomial,
    eval_buffer::VariableContainer,
    prob_params::IVPStatement,
    pt_trait::IVPVariationTrait,
    t_start::T,
    t_fin::T,
    config::IVPConfig,
    constraints_info::ConstraintsTrait,
    )::Symbol where T

    return solveIVP!(
        sol,
        eval_buffer,
        getIVPtrait(prob_params),
        prob_params,
        pt_trait,
        t_start,
        t_fin,
        config,
        constraints_info,
    )
end


# generate a piece-wise polynomials (seperately for each variable) that approximately solve an IVP of the form:
# - starts at t = 0, stop at t = t_fin,
# - h_initial used for adaption of the first polynomial.
# Subsequent h_initials are based on the solved step size for the previous polynomial segment.
function solveIVP!(
    sol::PiecewiseTaylorPolynomial,
    eval_buffer::VariableContainer,
    ::NumericalIVPTrait,
    prob_params::IVPStatement,
    pt_trait::IVPVariationTrait,
    t_start::T,
    t_fin::T,
    config::IVPConfig,
    constraints_info::ConstraintsTrait,
    )::Symbol where T

    # #set up.
    t_expansion = t_start

    # problem buffer.
    prob = getIVPbuffer(getbuffertype(prob_params), prob_params)

    # for adaptive step. Does not use parallel transport.
    test_IVP_buffer = copyIVPtestbuffer(getbuffertype(prob_params), prob_params)

    #resetsolution!(sol, prob_params.x0, prob_params.u0, prob_params.v0_set)
    resetsolution!(sol)

    # # solve for the first solution piece.
    t_next, instruction = solvesegmentIVP!(
        sol,
        test_IVP_buffer,
        eval_buffer,
        prob,
        pt_trait,
        t_expansion,
        config,
        constraints_info,
    )
    #h_initial = sol.steps[end] * 10 # heurestic, so that we have similar step sizes?

    # check stopping condition
    if t_next > t_fin || instruction != :continue_simulation
        return mapexitflag(t_next, t_fin, instruction)
    end

    # each iteration is a piece.
    for _ = 2:config.max_pieces 

        # set up new IVP problem for the current expansion time.
        prob_current = getIVPbuffer(
            getbuffertype(prob_params),
            prob_params,
            eval_buffer,
        )
        t_expansion = t_next

        # solve for the current solution piece.
        t_next, instruction = solvesegmentIVP!(
            sol,
            test_IVP_buffer,
            eval_buffer,
            prob_current,
            pt_trait,
            t_expansion,
            config,
            constraints_info,
        )
    
        # check stopping condition.
        if t_next > t_fin || instruction != :continue_simulation
            
            return mapexitflag(t_next, t_fin, instruction)
        end
    
    end

    return :max_piece_reached
end

function mapexitflag(t::Real, t_fin::Real, instruction::Symbol)::Symbol
    
    if t < t_fin
        return instruction
    end

    return :success
end

function wasvalidsession(exit_flag::Symbol)::Bool
    if exit_flag == :success || exit_flag == :intersection_found
        return true
    end

    return false
end

# exits with eval_buffer holding the solution evaluated at t_next = t_expansion + h.
function solvesegmentIVP!(
    sol::PiecewiseTaylorPolynomial,
    test_IVP_buffer::IVPBuffer,
    eval_buffer::VariableContainer,
    prob::IVPBuffer,
    pt_trait::IVPVariationTrait,
    t_expansion::T,
    config::IVPConfig,
    C::ConstraintsTrait,
    )::Tuple{T,Symbol} where T

    getfirstorder!(prob, pt_trait) # this brings the solution to order 1.

    # # solve for an appropriate step size, increasing the order according to the adaptation strategy in strategy_config.
    h, instruction = computetaylorsolution!(
        prob, 
        test_IVP_buffer,
        eval_buffer,
        pt_trait,
        t_expansion,
        config,
        C,
    )

    # # update solution, and prepare for the next IVP.
    t_next = t_expansion + h
    storesolutionpiece!(sol, eval_buffer, pt_trait, prob, t_expansion, h)
    
    return t_next, instruction
end

# see post_methods.jl for each family of IVPs for storesolutionpiece!() 

###### adaptive step


# no constraints, fixed order.
function computetaylorsolution!(
    prob::IVPBuffer,
    test_IVP_buffer::IVPBuffer,
    eval_buffer::VariableContainer,
    pt_trait::IVPVariationTrait,
    t0::T,
    config::FixedOrderConfig{T},
    ::NoConstraints,
    )::Tuple{T,Symbol} where T

    ### get to a high enough order so that we can start computing the error.
    #@show length(prob.x.c[1]), length(prob.u.c[1])
    # from the calling function, firstorder!() got prob.x and prob.u up to order 1 already. Start at order 2.
    for _ = 2:config.L
        increaseorder!(prob, pt_trait)
    end
    
    # h = choosestepsize(
    #     #VelocityContinuityStep(),
    #     GuentherWolfStep(),
    #     config.step_config.ϵ,
    #     prob;
    #     #GuentherWolfStep();
    #     h_max = config.step_config.h_max,
    #     step_reduction_factor = config.step_config.reduction_factor,
    # )
    h = choosestepsize!(
        test_IVP_buffer,
        eval_buffer,
        config.step_config.strategy,
        prob,
        t0,
        config.step_config,
    )

    return packagereturnstep(h, config.min_step)
end

function computetaylorsolution!(
    prob::IVPBuffer,
    test_IVP_buffer::IVPBuffer,
    eval_buffer::VariableContainer,
    pt_trait::IVPVariationTrait,
    t0::T,
    config::AdaptOrderConfig{T},
    C::ConstraintsTrait,
    )::Tuple{T,Symbol} where T

    ### set up

    #ϵ = config.ϵ
    L_min = config.L_min
    L_max = config.L_max
    order_increase_factor = config.order_increase_factor
    min_step = config.min_step
    #r_order = config.r_order
    #h_max = config.h_max
    #N_analysis_terms = config.N_analysis_terms
    #step_reduction_factor = config.reduction_factor
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
    # h = choosestepsize(
    #     VelocityContinuityStep(),
    #     ϵ,
    #     p;
    #     #GuentherWolfStep();
    #     h_max = h_max,
    #     step_reduction_factor = step_reduction_factor,
    # )
    h = choosestepsize!(
        test_IVP_buffer,
        eval_buffer,
        config.step_config.strategy,
        prob,
        t0,
        config.step_config,
    )

    h, constraint_ind = refinestep!(
        C,
        h,
        prob,
    )
    #@show constraint_ind

    if constraint_ind > 0
        # found intersection.
        return h, :intersection_found
    end
    
    #return h, :continue_simulation # debug.

    ### start adaption strategy.
    h_prev = zero(T)

    #r = r_order + 1 # initialize to a value so we satisfy the loop condition the first time.
    while getorder(prob) < L_max && h_prev*order_increase_factor < h #&& r > r_order

        h_prev = h

        # increase order.
        increaseorder!(prob, pt_trait)

        h = choosestepsize!(
            test_IVP_buffer,
            eval_buffer,
            config.step_config.strategy,
            prob,
            t0,
            config.step_config,
        )

        h, constraint_ind = refinestepnumerical!(
            C, h, h_prev, prob, 
        )

        #valid_h = (0 < h < h)
        valid_h = isfinite(h)

        # we avoid returning inside the if-else if we have no roots in [0, h] at the current order.
        if valid_h
            # if constraint_ind == 0
            #     # no roots on [0, h)

            #     return h, :continue_simulation
            # else
            #     # found intersection.

            #     return h, :stop_simulation
            # end

            if constraint_ind > 0
                # found intersection.
                return h, :intersection_found
            end
        
        else 
            # the inconclusive case and unexpected error case, like root solve failed.
            # no intersection so far, up to and including the previous order.

            if getorder(prob) > L_min

                decreaseorder!(prob, pt_trait, L_min)
            end
            # don't need to decrease if we're already at order: L_min.

            return packagereturnstep(h_prev, min_step)
        end

        # estimate the decrease in error by our increase in order.
        # see current_notes under `root find` folder.
        # r = computeerrorratio(
        #     p,
        #     ϵ,
        #     N_analysis_terms;
        #     h_max = h_max,
        #     step_reduction_factor = step_reduction_factor,
        # )

        # if r < r_order
        #     return h
        # end
    end

    return packagereturnstep(h, min_step)
end

function packagereturnstep(h::T, min_step::T)::Tuple{T,Symbol} where T
    if h > min_step
        return h, :continue_simulation
    end

    return h, :step_too_small
end


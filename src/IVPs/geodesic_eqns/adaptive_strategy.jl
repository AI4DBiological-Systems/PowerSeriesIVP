
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


############ step size


# eq:choose_ODE_step_size. An h_default of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsizeunivariate(
    ::VelocityContinuityStep,
    ϵ::T,
    p::GeodesicIVPBuffer,
    d::Integer;
    h_default = one(T),
    step_reduction_factor = 2,
    )::T where T

    c_u = p.u.c
    h = stepsizeformulaxu(ϵ, c_u[d][end], length(c_u) - 1)

    if isfinite(h)
        # further shorten step size as a heuristic strategy to reduce error.
        return convert(T, h/step_reduction_factor)
    end

    return convert(T, h_default/step_reduction_factor)
end

# also use the conservative factor of 1/2.
function stepsizeformulaxu(ϵ, x, order)
    return (ϵ/abs(x))^(1/order)
end

# eq:choose_ODE_step_size. An h_default of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsizeunivariate(
    ::GuentherWolfStep,
    ϵ::T,
    p::GeodesicIVPBuffer,
    d::Integer;
    h_default = one(T),
    step_reduction_factor = 2,
    )::T where T

    c = p.x.c[d]

    h = stepsizeformulaGuentherWolf(ϵ, c[end], length(c) -1)
    if isfinite(h)
        # further shorten step size as a heuristic strategy to reduce error.
        return convert(T, h/step_reduction_factor)
    end

    return convert(T, h_default/step_reduction_factor)
end

# see http://www.phys.uri.edu/nigh/NumRec/bookfpdf/f16-2.pdf for a constrast of the method from PSM 2019 with mainstream RK adaptive step size selection algorithms.
# we hardcoded such that only one analysis term for the step size determination.
# conservative factor of 1/2.
function stepsizeformulaGuentherWolf(ϵ::T, x::T, L::Integer)::T where T
    return (abs(ϵ/(2*x)))^(1/(L-1))
end

function choosestepsize(
    step_trait::NonProbingStep,
    ϵ::T,
    p::GeodesicIVPBuffer;
    h_default = one(T),
    step_reduction_factor = 2,
    )::T where T

    c_x = p.x.c

    min_h = convert(T, Inf)
    for d in eachindex(c_x)
        h = choosestepsizeunivariate(
            step_trait,
            ϵ,
            p,
            d;
            h_default = h_default,
            step_reduction_factor = step_reduction_factor,
        )
        min_h = min(min_h, h)
    end

    return min_h
end


function choosestepsize(
    ::AllNonProbingStep,
    ϵ::T,
    p::GeodesicIVPBuffer;
    h_default = one(T),
    step_reduction_factor = 2,
    )::T where T

    h1 = choosestepsize(
        GuentherWolfStep(),
        ϵ,
        p;
        h_default = h_default,
        step_reduction_factor = step_reduction_factor,
    )

    h2 = choosestepsize(
        GuentherWolfStep(),
        ϵ,
        p;
        h_default = h_default,
        step_reduction_factor = step_reduction_factor,
    )

    if !isfinite(h1)
        h1 = convert(T, Inf)
    end

    if !isfinite(h2)
        h2 = convert(T, Inf)
    end

    min_h = min(h1,h2)

    if isfinite(min_h)
        return min_h
    end

    return zero(T)
end

########## general.

# does not actually mutate any of the inputs.
function choosestepsize!(
    ::GeodesicIVPBuffer, # unused.
    ::GeodesicEvaluation{T}, # unused.
    ::NonProbingStep,
    prob::GeodesicIVPBuffer,
    t0,
    config::StepConfig,
    ) where T

    h = choosestepsize(
        config.strategy,
        config.ϵ,
        prob,
        h_default = config.h_default,
        step_reduction_factor = config.reduction_factor,
    )

    return h
end

########## continuity error

# mutates test_IVp and eval_buffer.
# if continuity conditions fail, mutates sol and eval_buffer.
function choosestepsize!(
    test_IVP::GeodesicIVPBuffer,
    eval_buffer::GeodesicEvaluation{T},
    ::DerivativeContinuityStep,
    prob::GeodesicIVPBuffer,
    t0::T,
    config::StepConfig,
    ) where T

    # set up.
    ϵ = config.ϵ
    discount_factor = config.discount_factor
    c_x = prob.x.c
    c_u = prob.u.c
    step_strategy = config.strategy

    # first try: use continuity of dx/dt at t=h vs. u(h).
    h = choosestepsize(
        step_strategy.first_step_strategy,
        ϵ,
        prob,
        h_default = config.h_default,
        step_reduction_factor = config.reduction_factor,
    )

    #return h

    # get test piece.
    t0_next = h + t0
    #evalsolution!(eval_buffer, DisableParallelTransport(), sol.coefficients[end], t0_next, t0)
    evalsolution!(eval_buffer, c_x, c_u, t0_next, t0)

    resetbuffer!(
        test_IVP,
        DisableParallelTransport(),
        eval_buffer.position,
        eval_buffer.velocity,
    )
    getfirstorder!(test_IVP, DisableParallelTransport()) # this brings the solution of test_IVP to order 1.
    increaseorder!(test_IVP, DisableParallelTransport()) # bring to order 2.

    # get error.
    err = getcontinuityerror(step_strategy, h, c_x, test_IVP.x.c)
    #@assert 3==4

    while err > ϵ
        # redo the test solution with a smaller current step.
        h = h * discount_factor
        
        # get test piece.
        t0_next = h + t0
        evalsolution!(eval_buffer, c_x, c_u, t0_next, t0)

        resetbuffer!(
            test_IVP,
            DisableParallelTransport(),
            eval_buffer.position,
            eval_buffer.velocity,
        )
        getfirstorder!(test_IVP, DisableParallelTransport()) # this brings the solution of test_IVP to order 1.
        increaseorder!(test_IVP, DisableParallelTransport()) # bring to order 2.

        # get error.
        err = getcontinuityerror(step_strategy, h, c_x, test_IVP.x.c)
    end

    return h
end


function getcontinuityerror(
    ::ContinuityFirstDerivative,
    h::T,
    c::Vector{Vector{T}},
    c_next::Vector{Vector{T}},
    )::T where T
    
    return getscalederror(h, c, c_next, 1)
end

function getcontinuityerror(
    ::ContinuitySecondDerivative,
    h::T,
    c::Vector{Vector{T}},
    c_next::Vector{Vector{T}},
    )::T where T
    
    err1 = getscalederror(h, c, c_next, 1)
    err2 = getscalederror(h, c, c_next, 1)
    return err1 + err2*2
end

function getcontinuityerror(
    step_strategy::ContinuityHigherDerivative,
    h::T,
    c::Vector{Vector{T}},
    c_next::Vector{Vector{T}},
    )::T where T
    
    return sum( getscalederror(h, c, c_next, m)*factorial(m) for m = 1:step_strategy.max_order )
end

# continuity for higher-order, m, derivatives.
function getscalederror(
    h::T,
    c::Vector{Vector{T}},
    c_next::Vector{Vector{T}},
    m::Integer,
    ) where T
    
    @assert length(c) == length(c_next)
    @assert !isempty(c)

    total_order = length(c[begin])-1
    @assert m < total_order

    max_scaled_err = zero(T)

    for d in eachindex(c)

        scaled_difference = c[d][begin+m]-c_next[d][begin+m]
        for n = m+1:total_order
            scaled_difference += binomial(n,m)*c[d][begin+n]*h^(n-m) # TODO: use bino_mat look up table here
        end

        #@show abs(scaled_difference)
        max_scaled_err = max(abs(scaled_difference), max_scaled_err)
    end

    return max_scaled_err
end

# # continuity for derivative.
# function continuitycheck1st(
#     next_u0::Vector{T},
#     next_x::Vector{Vector{T}},
#     ) where T

#     # joint checks: 1st order of x should be 0th order of u. If the first-order IVP we solve was derived from a higher-order IVP, then we need to do this type of check.
#     #err = maximum( abs(sol_eval.velocity[d] - next_piece.x[d][begin+1]) for d in eachindex(sol_eval.velocity) )
#     err = maximum( abs(next_u0[d] - next_x[d][begin+1]) for d in eachindex(next_u0) )
    
#     return err
# end
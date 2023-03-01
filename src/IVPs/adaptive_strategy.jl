
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
    step_trait::StepStrategyTrait,
    p::GeodesicIVPBuffer,
    ϵ::T,
    N_analysis_terms::Integer;
    h_max = one(T),
    step_reduction_factor = 2,
    )::T where T
    
    order = length(p.x.c[begin])
    @assert order > N_analysis_terms > 0

    M = order - N_analysis_terms
    
    a = choosestepsize(
        step_trait,
        #GuentherWolfStep(),
        ϵ,
        p;
        order = M, # I am here.
        h_max = h_max,
        step_reduction_factor = step_reduction_factor,
    )

    b = choosestepsize(
        step_trait,
        #GuentherWolfStep(),
        ϵ,
        x;
        order = M-1, # I am here
        h_max = h_max,
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



############ step size


# eq:choose_ODE_step_size. An h_max of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsizeunivariate(
    ::VelocityContinuityStep,
    ϵ::T,
    p::GeodesicIVPBuffer,
    d::Integer;
    h_max = one(T),
    step_reduction_factor = 2,
    )::T where T

    c_u = p.u.c
    h = stepsizeformulaxu(ϵ, c_u[d][end], length(c_u) - 1)

    if isfinite(h)
        # further shorten step size as a heuristic strategy to reduce error.
        return convert(T, h/step_reduction_factor)
    end

    return convert(T, h_max/step_reduction_factor)
end

# also use the conservative factor of 1/2.
function stepsizeformulaxu(ϵ, x, order)
    return (ϵ/abs(x))^(1/order)
end

# eq:choose_ODE_step_size. An h_max of Inf means if no higher-order errors are computed, then we assume the model is exact, thus any step is valid. This means the maximum step for which the Taylor polynomial is valid is Inf.
function choosestepsizeunivariate(
    ::GuentherWolfStep,
    ϵ::T,
    p::GeodesicIVPBuffer,
    d::Integer;
    h_max = one(T),
    step_reduction_factor = 2,
    )::T where T

    c = p.x.c[d]

    h = stepsizeformulaGuentherWolf(ϵ, c[end], length(c) -1)
    if isfinite(h)
        # further shorten step size as a heuristic strategy to reduce error.
        return convert(T, h/step_reduction_factor)
    end

    return convert(T, h_max/step_reduction_factor)
end

# see http://www.phys.uri.edu/nigh/NumRec/bookfpdf/f16-2.pdf for a constrast of the method from PSM 2019 with mainstream RK adaptive step size selection algorithms.
# we hardcoded such that only one analysis term for the step size determination.
# conservative factor of 1/2.
function stepsizeformulaGuentherWolf(ϵ::T, x::T, L::Integer)::T where T
    return (abs(ϵ/(2*x)))^(1/(L-1))
end

function choosestepsize(
    step_trait::StepStrategyTrait,
    ϵ::T,
    p::GeodesicIVPBuffer;
    h_max = one(T),
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
            h_max = h_max,
            step_reduction_factor = step_reduction_factor,
        )
        min_h = min(min_h, h)
    end

    return min_h
end


########### continuity error

# if continuity conditions fail, mutates sol and eval_buffer.
function continuitycheck!(
    sol::PiecewiseTaylorPolynomial{T},
    eval_buffer::GeodesicEvaluation{T},
    #prob::GeodesicIVPBuffer,
    pt_trait::PT,
    metric_params::MT,
    config::ContinuityConfig{T},
    ) where {T,PT,MT}

    min_h = config.min_h
    discount_factor = config.discount_factor

    # get initial conditions for the next piece.
    t0_current = sol.t_expansion[end]
    h_current = sol.steps[end]
    t0_next = t0_current + h_current
    evalsolution!(eval_buffer, pt_trait, sol.coefficients[end], t0_next, t0_current)

    # solve the test solution for the next piece.
    prob = getivpbuffer(
        metric_params,
        eval_buffer.position,
        eval_buffer.velocity,
        eval_buffer.vector_fields,
    )
    getfirstorder!(prob, pt_trait) # this brings the solution to order 1.
    
    # check if x_dot_current(t0_next) and u0 agrees.
    pass_flag = continuitycheck(eval_buffer.velocity, prob.x.c, config)
    
    ϵ = config.ϵ
    err = convert(T, Inf)
    while err > ϵ && h_current > min_h
        # redo the test solution with a smaller current step.

        h_current = h_current * discount_factor
        
        prob = getivpbuffer(
            metric_params,
            eval_buffer.position,
            eval_buffer.velocity,
            eval_buffer.vector_fields,
        )
        getfirstorder!(prob, pt_trait)
        
        #pass_flag = continuitycheckxuunsimplified(eval_buffer.velocity, prob.x.c, config)
        err = continuitymaxerrorhigher(h, c, c_next, m)
        # I am here.
        
    end

    if h_current < min_h
        return h_urrent, prob, :try_again_with_decreased_order_need_to_get_starting_h_again
    end

    return h_current, prob, :clear_to_proceed
end



# function continuitycheckxuunsimplified(
#     next_u0::Vector{T},
#     next_x::Vector{Vector{T}},
#     config::ContinuityConfig{T},
#     ) where T

#     # joint checks: 1st order of x should be 0th order of u. If the first-order IVP we solve was derived from a higher-order IVP, then we need to do this type of check.
#     #err = maximum( abs(sol_eval.velocity[d] - next_piece.x[d][begin+1]) for d in eachindex(sol_eval.velocity) )
#     err = maximum( abs(next_u0[d] - next_x[d][begin+1]) for d in eachindex(next_u0) )
#     if err > config.ϵ
#         return false
#     end

#     return true
# end

# for higher-order, m, derivatives.
function continuitymaxerrorhigher(
    h::T,
    c::Vector{Vector{T}}
    c_next::Vector{Vector{T}},
    m::Integer,
    ) where T
    
    @assert length(x) == length(x_next)

    total_order = length(c[d])-1
    @assert m < total_order

    max_err = zero(T)

    for d in eachindex(x)

        err = c[d][begin+m]-c_next[d][begin+m]
        for n = m+1:total_order
            err += binomial(n,m)*c[d][begin+n]*h^(n-m) # TODO: use bino_mat look up table here
        end

        err = abs(err)
        max_err = max(err, min_err)
    end

    return max_err
end

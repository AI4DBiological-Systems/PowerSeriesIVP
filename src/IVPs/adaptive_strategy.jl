
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



############ step size

# # try all methods, and choose smallest.
# function choosestepsize(
#     ϵ::T,
#     p::GeodesicIVPBuffer,
#     ::MinAllStep;
#     h_max = one(T),
#     step_reduction_factor = 2,
#     )::T where T

#     h1 = choosestepsize(
#         ϵ,
#         p,
#         XUStepStragey();
#         h_max = h_max,
#         step_reduction_factor = step_reduction_factor,
#     )

#     h2 = choosestepsize(
#         ϵ,
#         p,
#         GuentherWolfStep();
#         h_max = h_max,
#         step_reduction_factor = step_reduction_factor,
#     )

#     h = min(h2, h1)

#     return h
# end

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
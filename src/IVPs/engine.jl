
# struct IVPConfig{T}
#     ϵ::T
#     h_initial::T
#     L_test_max::Int
#     r_order::T
#     h_zero_error::T
#     step_reduction_factor::T
#     max_pieces::Int
#     N_analysis_terms::Int
# end

# function IVPConfig(
#     ::Type{T};
#     ϵ::T = convert(T, 1e-6),
#     h_initial = one(T),
#     L_test_max::Int = convert(Int, 10),
#     r_order = convert(T, 0.3),
#     h_zero_error = convert(T, Inf),
#     step_reduction_factor = convert(T, 2),
#     max_pieces::Int = typemax(Int),
#     N_analysis_terms::Int = convert(Int, 2),
#     ) where T

#     return IVPConfig(
#         convert(T, ϵ),
#         convert(T, h_initial),
#         convert(Int, L_test_max),
#         convert(T, r_order),
#         convert(T, h_zero_error),
#         convert(T, step_reduction_factor),
#         convert(Int, max_pieces),
#         convert(Int, N_analysis_terms),
#     )
# end

struct GeodesicPiece{T}
    
    # [variable index][order index]
    x::Vector{Vector{T}} # coefficients for the solution state.
    u::Vector{Vector{T}} # coefficients for first-derivative of state.
    vs::Vector{Vector{Vector{T}}} # i-th entry contain the coefficients for the i-th parallel vector field.
end

function GeodesicPiece(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    parallel_transport::Vector,
    )::GeodesicPiece{T} where T
    
    vs::Vector{Vector{Vector{T}}} = collect( parallel_transport[m].v.c for m in eachindex(parallel_transport) )

    return GeodesicPiece(x, u, vs)
end

function GeodesicPiece(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    )::GeodesicPiece{T} where T
    
    return GeodesicPiece(x, u, Vector{Vector{Vector{T}}}(undef, 0))
end

function getNvars(c::GeodesicPiece)::Int
    return length(c.x)
end

struct PiecewiseTaylorPolynomial{T}

    # [piece index]
    coefficients::Vector{GeodesicPiece{T}}

    # [piece index]
    expansion_points::Vector{T}
    steps::Vector{T} # diagnostic information. Use this to determine if t_fin was actually reached by our solver algorithm for a given ODE solution.

    # the interval for the simulation. Do expansion_points + steps for the actual simulated time.
    #t_start::T
    #t_fin::T
    # 3-element 1-D array:
    #time_endpoints::Vector{T} # first entry is start time, last is finish time.
    # the second interval is the simulated last time.

    # the initial conditions.
    starting_position::Vector{T}
    starting_velocity::Vector{T}
    starting_vectors::Vector{Vector{T}} # transport vectors.
end

function getNtransports(A::PiecewiseTaylorPolynomial)
    return length(A.coefficients[begin].vs)
end

function getendtime(sol::PiecewiseTaylorPolynomial{T})::T where T
    return sol.expansion_points[end] + sol.steps[end]
end

function getstarttime(sol::PiecewiseTaylorPolynomial{T})::T where T
    return sol.expansion_points[begin]
end

function getNvars(A::PiecewiseTaylorPolynomial)::Int
    return length(A.coefficients[begin].x)
end

function resetsolution!(
    A::PiecewiseTaylorPolynomial{T},
    x0::Vector{T},
    u0::Vector{T},
    v0_set::Vector{Vector{T}},
    ) where T

    resize!(A.coefficients, 0)
    resize!(A.expansion_points, 0)
    resize!(A.steps, 0)

    A.starting_position[:] = x0
    A.starting_velocity[:] = u0

    resize!(A.starting_vectors, length(v0_set))
    for m in eachindex(v0_set)
        A.starting_vectors[m] = v0_set[m]
    end

    return nothing
end

# straight line as a PiecewiseTaylorPolynomial{T,GeodesicPiece{T}}
function createline(position::Vector{T}, velocity::Vector{T}, t_start::T, t_fin::T) where T <: AbstractFloat
    
    D = length(position)
    @assert length(velocity) == D

    # single piece.

    # assemble coefficients
    x = collect( [position[d]; velocity[d]] for d = 1:D ) # starting position, velocity.
    u = collect( [velocity[d]; zero(T)] for d = 1:D ) # velocity, acceleration.

    coefficients = Vector{GeodesicPiece{T}}(undef, 0)
    push!(coefficients, GeodesicPiece(x,u))

    # time.
    expansion_points = [t_start;]
    steps::Vector{T} = [ (t_fin-t_start)*1.1 ;] # 1.1 instead of 1.0 so that the end point t_fin is included in the interval of validity for the solution piece.

    return PiecewiseTaylorPolynomial(
        coefficients,
        expansion_points,
        steps,
        position,
        velocity,
        Vector{Vector{T}}(undef, 0),
    )
end

struct GeodesicEvaluation{T}
    position::Vector{T}
    velocity::Vector{T}
    vector_field::Vector{Vector{T}}
end

function GeodesicEvaluation(::Type{T}, N_vars::Integer, N_transports::Integer)::GeodesicEvaluation{T} where T
    return GeodesicEvaluation(
        ones(T, N_vars) .* NaN,
        ones(T, N_vars) .* NaN,
        collect( ones(T, N_vars) .* NaN for _ = 1:N_transports ),
    )
end

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::GeodesicEvaluation{T},
    _::DisableParallelTransport,
    c::GeodesicPiece{T},
    t,
    a,
    ) where T

    @assert length(c.x) == length(c.u) == length(out.position) == length(out.velocity)

    for d in eachindex(c.x)
        out.position[d] = evaltaylor(c.x[d], t, a)
        out.velocity[d] = evaltaylor(c.u[d], t, a)
    end

    return nothing
end

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::GeodesicEvaluation{T},
    _::EnableParallelTransport,
    c::GeodesicPiece{T},
    t,
    a,
    ) where T

    evalsolution!(out, DisableParallelTransport(), c, t, a)

    for m in eachindex(c.vs)
        for d in eachindex(c.vs[m])
            out.vector_field[m][d] = evaltaylor(c.vs[m][d], t, a)
        end
    end

    return nothing
end

# handles the selection of a solution piece from the piece-wise solution.
function evalsolution!(
    out::GeodesicEvaluation{T},
    pt_trait::PT,
    A::PiecewiseTaylorPolynomial,
    t::T,
    ) where {PT<:ParallelTransportTrait, T}

    expansion_points = A.expansion_points
    t_start = getstarttime(A)
    t_fin = getendtime(A)

    if !(t_start <= t <= t_fin)

        return false
    end

    for k in Iterators.drop(eachindex(expansion_points), 1)
    
        if t < expansion_points[k]
            
            evalsolution!(out, pt_trait, A.coefficients[k-1], t, expansion_points[k-1])
            return true
        end
    end

    if t < expansion_points[end] + A.steps[end]

        evalsolution!(out, pt_trait, A.coefficients[end], t, expansion_points[end])
        return true
    end

    # case: our solver algorithm did not reach t_fin, and t is beyond the last solution piece's estimated interval of validity.
    return false
end


function evalsolution(
    pt_trait::PT,
    A::PiecewiseTaylorPolynomial{T},
    t,
    )::Tuple{GeodesicEvaluation{T}, Bool} where {PT,T}

    out = GeodesicEvaluation(T, getNvars(A), getNtransports(A))
    status_flag = evalsolution!(out, pt_trait, A, t)

    return out, status_flag
end

function batchevalsolution!(
    status_flags::BitVector,
    positions_buffer::Vector{Vector{T}},
    velocities_buffer::Vector{Vector{T}},
    A::PiecewiseTaylorPolynomial{T},
    ts,
    ) where T

    @assert length(positions_buffer) == length(ts) == length(velocities_buffer) == length(status_flags)
    out = GeodesicEvaluation(T, getNvars(A), getNtransports(A))

    for n in eachindex(ts)
        status_flags[n] = evalsolution!(out, DisableParallelTransport(), A, ts[n]) # TODO, return error flags or which evals were valid.
        
        positions_buffer[n][:] = out.position
        velocities_buffer[n][:] = out.velocity
    end

    return nothing
end

function batchevalsolution!(
    status_flags::BitVector,
    positions_buffer::Vector{Vector{T}}, # [eval_index][dimension].
    velocities_buffer::Vector{Vector{T}}, # [eval_index][dimension].
    vector_fields_buffer::Vector{Vector{Vector{T}}}, # [eval index][vector field index][dimension]
    A::PiecewiseTaylorPolynomial{T},
    ts,
    ) where T

    @assert length(positions_buffer) == length(ts) == length(velocities_buffer) == length(status_flags)
    out = GeodesicEvaluation(T, getNvars(A), getNtransports(A))

    for n in eachindex(ts)
        status_flags[n] = evalsolution!(out, EnableParallelTransport(), A, ts[n]) # TODO, return error flags or which evals were valid.
        
        positions_buffer[n][:] = out.position
        velocities_buffer[n][:] = out.velocity
        
        for m in eachindex(out.vector_field)
            vector_fields_buffer[n][m][:] = out.vector_field[m]
        end
    end

    return nothing
end

###################

function tautology(args...)::Bool
    return true
end

function contradiction(args...)::Bool
    return false
end

function solveIVP(
    prob_params::GeodesicIVPProblem{MT,T},
    pt_trait::PT,
    t_start::T,
    t_fin::T,
    config::IVPConfig;
    shouldstop = contradiction, # this function takes at least an input that is of type GeodesicEvaluation{T}, and returns true if we should stop simulating the forward trajectory of the IVP.
    h_initial::T = convert(T, (t_fin-t_start)/10),
    )::PiecewiseTaylorPolynomial{T} where {T, MT<:MetricParams, PT<:ParallelTransportTrait}

    sol = PiecewiseTaylorPolynomial(
        Vector{GeodesicPiece{T}}(undef,0),
        Vector{T}(undef,0),
        Vector{T}(undef,0),
        prob_params.x0,
        prob_params.u0,
        prob_params.v0_set,
    )

    next_conditions = GeodesicEvaluation(
        T,
        getNvars(prob_params),
        getNtransports(prob_params),
    )

    solveIVP!(
        sol,
        next_conditions,
        prob_params,
        pt_trait,
        t_start,
        t_fin,
        config;
        shouldstop = shouldstop,
        h_initial = h_initial,
    )

    return sol
end

# generate a piece-wise polynomials (seperately for each variable) that approximately solve an IVP of the form:
# - starts at t = 0, stop at t = t_fin,
# - h_initial used for adaption of the first polynomial.
# Subsequent h_initials are based on the solved step size for the previous polynomial segment.
function solveIVP!(
    sol::PiecewiseTaylorPolynomial{T},
    next_conditions::GeodesicEvaluation{T},
    prob_params::GeodesicIVPProblem{MT,T},
    pt_trait::PT,
    t_start::T,
    t_fin::T,
    config::IVPConfig;
    shouldstop = contradiction, # this function takes at least an input that is of type GeodesicEvaluation{T}, and returns true if we should stop simulating the forward trajectory of the IVP.
    h_initial::T = convert(T, (t_fin-t_start)/10),
    ) where {T, MT<:MetricParams, PT<:ParallelTransportTrait}

    # #set up.
    t_expansion = t_start

    # a = prob_params.a
    # b = prob_params.b
    metric_params = prob_params.metric_params
    prob = getivpbuffer(
        metric_params,
        prob_params.x0,
        prob_params.u0,
        prob_params.v0_set,
    )

    resetsolution!(sol, prob_params.x0, prob_params.u0, prob_params.v0_set)

    # # solve for the first solution piece.
    t_next = solvesegmentIVP!(
        sol,
        next_conditions,
        prob,
        pt_trait,
        h_initial,
        t_expansion,
        config,
    )
    #h_initial = sol.steps[end] * 10 # heurestic, so that we have similar step sizes?

    # check stopping condition
    if t_next > t_fin || shouldstop(sol, next_conditions)
        return sol
    end

    # each iteration is a piece.
    for _ = 2:config.max_pieces 

        # set up new IVP problem for the current expansion time.
        prob_current = getivpbuffer(
            metric_params,
            next_conditions.position,
            next_conditions.velocity,
            next_conditions.vector_field,
        )
        t_expansion = t_next

        # solve for the current solution piece.
        t_next = solvesegmentIVP!(
            sol,
            next_conditions,
            prob_current,
            pt_trait,
            h_initial, # more direct control over how large the step sizes could be. set this to be large to encourage taking a big step and high-order.
            t_expansion,
            #sol.steps[end], # use the step from the last solution piece as the initial segment. actually, might not make sense for the geodesic setting, where extreme accuracy might not be needed.
            config,
        )
    
        # check stopping condition.
        if t_next > t_fin || shouldstop(sol, next_conditions)
            return sol
        end
    
    end

    return sol
end


# mutates sol and next_conditions.
function storesolutionpiece!(
    sol::PiecewiseTaylorPolynomial,
    next_conditions::GeodesicEvaluation{T},
    pt_trait::PT,
    prob::GeodesicIVPBuffer,
    t_expansion::T,
    h::T,
    ) where {PT,T}

    # add the coefficients for the solution piece.
    new_coefficients = GeodesicPiece(prob.x.c, prob.u.c, prob.parallel_transport)
    push!(sol.coefficients, new_coefficients)
    push!(sol.expansion_points, t_expansion)
    push!(sol.steps, h)

    # get the initial conditions for the next IVP that the next solution piece solves.
    t_next = t_expansion + h
    evalsolution!(next_conditions, pt_trait, new_coefficients, t_next, t_expansion)
    
    return nothing
end

# exits with next_conditions holding the solution evaluated at t_next = t_expansion + h.
function solvesegmentIVP!(
    sol::PiecewiseTaylorPolynomial{T},
    next_conditions::GeodesicEvaluation{T},
    prob::GeodesicIVPBuffer,
    pt_trait::ParallelTransportTrait,
    h_initial::T,
    t_expansion::T,
    config,
    ) where T

    # # solve for an appropriate step size, increasing the order according to the adaptation strategy in strategy_config.
    # h = computetaylorsolution!(
    #     prob,
    #     pt_trait,
    #     config.h_initial;
    #     ϵ = config.ϵ,
    #     L_test_max = config.L_test_max,
    #     r_order = config.r_order,
    #     h_zero_error = config.h_zero_error,
    #     N_analysis_terms = config.N_analysis_terms,
    # )
    
    getfirstorder!(prob, pt_trait) # this brings the solution to order 1.
    h = applyadaptionstrategy!(prob, pt_trait, h_initial, config)

    # # shorten step size as a heuristic strategy to reduce error.
    h = h/config.step_reduction_factor

    # # update solution, and prepare for the next IVP.
    t_next = t_expansion + h
    storesolutionpiece!(sol, next_conditions, pt_trait, prob, t_expansion, h)
    
    return t_next
end


####  front end methods for use with other parts of this library.

# ::Type{BT} is for use with trait function getbuffertype(), which is defined in Levi-Civita metric scripts like RQ22.jl
function getIVPbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    )::BT where BT <: GeodesicIVPBuffer
    
    return getIVPbuffer(
        prob_params.metric_params,
        prob_params.x0,
        prob_params.u0,
        prob_params.v0_set,
    )
end

function getIVPbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    initial_conditions::GeodesicEvaluation,
    )::BT where BT <: GeodesicIVPBuffer

    return getIVPbuffer(
        prob_params.metric_params,
        # copy(initial_conditions.position),
        # copy(initial_conditions.velocity),
        # collect( copy(initial_conditions.vector_fields[i]) for i in eachindex(initial_conditions.vector_fields) )
        initial_conditions.position,
        initial_conditions.velocity,
        initial_conditions.vector_fields,
    )
end

function copyIVPtestbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    )::BT where BT <: GeodesicIVPBuffer
    
    return getIVPbuffer(
        prob_params.metric_params,
        copy(prob_params.x0),
        copy(prob_params.u0),
        empty(prob_params.v0_set), # no parallel transport. test is on x, which only needs x and u to be updated.
    )
end

# used only for continuity test on x, which doesn't need the computation involving the parallel transported vector field solutions.
# therefore, we ignore resetting the parallel transport field to reduce computation.
function resetbuffer!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    x0::Vector{T},
    u0::Vector{T},
    ) where T <: AbstractFloat

    resetbuffer!(p.θ)
    resetbuffer!(p.u)
    resetbuffer!(p.x)

    p.x0[:] = x0
    p.u0[:] = u0

    return nothing
end

############### solution (PiecewiseTaylorPolynomials) methods

# mutates sol and eval_buffer.
function storesolutionpiece!(
    sol::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    eval_buffer::GeodesicEvaluation{T},
    pt_trait::ParallelTransportTrait,
    prob::GeodesicIVPBuffer,
    t_expansion::T,
    h::T,
    ) where T

    # add the coefficients for the solution piece.
    new_coefficients = GeodesicPiece(prob.x.c, prob.u.c, prob.parallel_transport)
    push!(sol.coefficients, new_coefficients)
    push!(sol.expansion_points, t_expansion)
    push!(sol.steps, h)

    # get the initial conditions for the next IVP that the next solution piece solves.
    t_next = t_expansion + h
    evalsolution!(eval_buffer, pt_trait, new_coefficients, t_next, t_expansion)
    
    return nothing
end

############ constraint intersection-related.

function refinestep!(
    C::ConstraintsContainer{T},
    h::T,
    p::GeodesicIVPBuffer,
    )::Tuple{T,Int} where T <: AbstractFloat

    return refinestep!(C, h, p.x.c)
end

function refinestepnumerical!(
    C::ConstraintsContainer{T},
    h::T,
    h_prev::T, # previously known good step size with no intersections.    
    p::GeodesicIVPBuffer,
    )::Tuple{T,Int} where T <: AbstractFloat

    return refinestepnumerical!(C, h, h_prev, p.x.c)
end



######### evaluate solution.

# same as evalsolution!(), DisableParallelTransport() version.
function evalsolution!(
    out::GeodesicEvaluation{T},
    c_x::Vector{Vector{T}},
    c_u::Vector{Vector{T}},
    t,
    a,
    ) where T

    @assert length(c_x) == length(c_u) == length(out.position) == length(out.velocity)

    for d in eachindex(c_x)
        out.position[d] = evaltaylor(c_x[d], t, a)
        out.velocity[d] = evaltaylor(c_u[d], t, a)
    end

    return nothing
end

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::GeodesicEvaluation{T},
    ::DisableParallelTransport,
    c::GeodesicPiece{T},
    t,
    a,
    ) where T

    # @show length(c.x)
    # @show length(out.position)
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
    ::EnableParallelTransport,
    c::GeodesicPiece{T},
    t,
    a,
    ) where T

    evalsolution!(out, DisableParallelTransport(), c, t, a)

    for m in eachindex(c.vs)
        for d in eachindex(c.vs[m])
            out.vector_fields[m][d] = evaltaylor(c.vs[m][d], t, a)
        end
    end

    return nothing
end


######### solution analysis.

# used in solution analysis.
function extractcoefficients(
    A::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    )::Tuple{Vector{Vector{T}},Vector{Vector{T}}} where T

    c_x = collect( A.coefficients[k].x for k in eachindex(A.coefficients) )
    c_u = collect( A.coefficients[k].u for k in eachindex(A.coefficients) )

    return c_x, c_u
end

function getNtransports(A::PiecewiseTaylorPolynomial{T, GeodesicPiece{T}}) where T
    return getNtransports(A.coefficients[begin])
end
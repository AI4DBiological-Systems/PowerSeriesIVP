
############## methods for the entire solution container.

function resetsolution!(
    A::PiecewiseTaylorPolynomial,
    # x0::Vector{T},
    # u0::Vector{T},
    # v0_set::Vector{Vector{T}},
    )

    resize!(A.coefficients, 0)
    resize!(A.expansion_points, 0)
    resize!(A.steps, 0)

    # A.starting_position[:] = x0
    # A.starting_velocity[:] = u0

    # resize!(A.starting_vectors, length(v0_set))
    # for m in eachindex(v0_set)
    #     A.starting_vectors[m] = v0_set[m]
    # end

    return nothing
end

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

################# methods for evalating solution trajectory.

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
    pt_trait::ParallelTransportTrait,
    A::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    t,
    )::Tuple{GeodesicEvaluation{T}, Bool} where T

    out = GeodesicEvaluation(T, getdim(A), getNtransports(A))
    status_flag = evalsolution!(out, pt_trait, A, t)

    return out, status_flag
end

function batchevalsolution!(
    status_flags::BitVector,
    positions_buffer::Vector{Vector{T}},
    velocities_buffer::Vector{Vector{T}},
    A::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    ts,
    ) where T

    @assert length(positions_buffer) == length(ts) == length(velocities_buffer) == length(status_flags)
    out = GeodesicEvaluation(T, getdim(A), getNtransports(A))

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
    A::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    ts,
    ) where T

    @assert length(positions_buffer) == length(ts) == length(velocities_buffer) == length(status_flags)
    out = GeodesicEvaluation(T, getdim(A), getNtransports(A))

    for n in eachindex(ts)
        status_flags[n] = evalsolution!(out, EnableParallelTransport(), A, ts[n]) # TODO, return error flags or which evals were valid.
        
        positions_buffer[n][:] = out.position
        velocities_buffer[n][:] = out.velocity
        
        for m in eachindex(out.vector_fields)
            vector_fields_buffer[n][m][:] = out.vector_fields[m]
        end
    end

    return nothing
end
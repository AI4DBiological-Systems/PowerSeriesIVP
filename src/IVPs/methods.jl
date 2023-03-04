
############## methods for the entire solution container.

function resetsolution!(A::PiecewiseTaylorPolynomial)

    resize!(A.coefficients, 0)
    resize!(A.expansion_points, 0)
    resize!(A.steps, 0)

    return nothing
end


################# methods for evalating solution trajectory.

# no generic version for evalsolution!(), but make a setupevalbuffer() generic function.

# I am here. fix this up to use traits, and make it generic.
function evalsolution(
    pt_trait::ParallelTransportTrait,
    A::PiecewiseTaylorPolynomial{T,GeodesicPiece{T}},
    t,
    )::Tuple{GeodesicEvaluation{T}, Bool} where T

    out = GeodesicEvaluation(T, getdim(A), getNtransports(A))
    status_flag = evalsolution!(out, pt_trait, A, t)

    return out, status_flag
end

# I am here. make a generic version of this.
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

# handles the selection of a solution piece from the piece-wise solution.
function evalsolution!(
    out::VariableContainer,
    pt_trait::IVPVariationTrait,
    A::PiecewiseTaylorPolynomial,
    t,
    )::Bool

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

############## methods for the entire solution container.

function resetsolution!(A::PiecewiseTaylorPolynomial)

    resize!(A.coefficients, 0)
    resize!(A.expansion_points, 0)
    resize!(A.steps, 0)

    return nothing
end


################# methods for evalating solution trajectory.

# no generic version for evalsolution!(), but make a setupevalbuffer() generic function.

function allocatevariablecontainer(::Type{VC}, X)::VC where VC <: VariableContainer
    return allocatevariablecontainer(X)
end

function evalsolution(
    ::Type{VT},
    pt_trait::IVPVariationTrait,
    A::PiecewiseTaylorPolynomial,
    t,
    )::Tuple{VT, Bool} where VT <: VariableContainer

    out = allocatevariablecontainer(A)
    status_flag = evalsolution!(out, pt_trait, A, t)

    return out, status_flag
end

function batchevalsolution!(
    status_flags::BitVector,
    out::Vector{VT}, # [eval_index]
    variation_trait::IVPVariationTrait,
    A::PiecewiseTaylorPolynomial,
    ts::AbstractArray,
    ) where VT <: VariableContainer

    N_evals = length(ts)
    resize!(out, N_evals)
    resize!(status_flags, N_evals)

    for n in eachindex(ts)
        status_flags[n] = evalsolution!(out[n], variation_trait, A, ts[n])
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
    t_fin = t_start + getsimulationinterval(A)

    if !(t_start <= t <= t_fin)
        
        return false
    end

    for k in Iterators.drop(eachindex(expansion_points), 1)
    
        if t < expansion_points[k]
            
            evalsolution!(out, pt_trait, A.coefficients[k-1], t, expansion_points[k-1])
            return true
        end
    end

    if t > expansion_points[end] + A.steps[end]
        # case: our solver algorithm did not reach t_fin, and t is beyond the last solution piece's estimated interval of validity.
        return false
    end

    evalsolution!(out, pt_trait, A.coefficients[end], t, expansion_points[end])
    return true
end
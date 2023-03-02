
# outputs allocated multivariate coefficients for a curve, [dimension][order-1 index].
function allocatecoefficients(::Type{T}, N_vars::Integer)::Vector{Vector{T}} where T
    
    return collect( Vector{T}(undef, 0) for _ = 1:N_vars )
end

function resetcoefficients!(x::Vector{Vector{T}})::Nothing where T

    for d in eachindex(x)
        resize!(x[d], 0)
    end

    return nothing
end
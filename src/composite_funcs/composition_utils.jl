function initializecoefficients(::Type{T}, N_vars::Integer)::Vector{Vector{T}} where T
    return collect( Vector{T}(undef, 0) for _ = 1:N_vars )
end
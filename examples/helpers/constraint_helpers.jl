
function generateHyperplaneConstraints(
    N_constraints::Integer,
    )
    as = collect( randn(N_vars) for _ = 1:N_constraints )
    bs = randn(N_constraints) .- 25
    
    for n in eachindex(as)
        if as[n][end] > 0
            as[n] = -as[n]
            bs[n] = -bs[n]
        end
    end

    return as, bs
end

function generateHyperplaneConstraintscase1(::Type{T}) where T

    as = Vector{Vector{T}}(undef, 2)
    as[1] = [-0.7262044632757468, -0.8138894370909177, -0.6104189261116074]
    as[2] = [-2.0501294946950264, -0.18095967976913974, -0.9406747232875855]

    bs = Vector{T}(undef, 2)
    bs[1] = 23.95929560119815
    bs[2] = 25.147764931652375

    return as, bs
end

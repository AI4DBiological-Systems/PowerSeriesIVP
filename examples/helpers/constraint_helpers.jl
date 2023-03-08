
function generateBoundConstraints(::Type{T}, D::Integer) where T

    lbs = ones(T, D)
    fill!(lbs, convert(T,randn()))

    ubs = ones(T, D)
    fill!(ubs, convert(T,randn()))

    if lbs[begin] > ubs[begin]
        lbs, ubs = ubs, lbs
    end

    return lbs, ubs
end

# make sure the generated hyperplanes form a convex polyhedron.
# must contain the origin of the coordinate system.
function generatecvxpolyhedron(::Type{T}, N::Integer, D::Integer; interior_pt = zeros(T,D)) where T

    as = Vector{Vector{T}}(undef, N)
    bs = Vector{T}(undef, N)

    for m in eachindex(as)
        
        as[m] = rand(D)
        bs[m] = randn()

        while dot(as[m], interior_pt) > bs[m]
            as[m] = rand(D)
            bs[m] = randn()
        end
    end

    return as, bs
end

#### specific constraints for testing.

function generateHyperplaneConstraintscase1(::Type{T}) where T

    as = Vector{Vector{T}}(undef, 2)
    as[1] = [-0.7262044632757468, -0.8138894370909177, -0.6104189261116074]
    as[2] = [-2.0501294946950264, -0.18095967976913974, -0.9406747232875855]

    bs = Vector{T}(undef, 2)
    bs[1] = 23.95929560119815
    bs[2] = 25.147764931652375

    return as, bs
end

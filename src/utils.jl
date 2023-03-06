
# function tautology(args...)::Bool
#     return true
# end

# function contradiction(args...)::Bool
#     return false
# end

#= # straight line.
function createline(position::Vector{T}, velocity::Vector{T}, t_start::T, t_fin::T) where T <: AbstractFloat
    
    D = length(position)
    @assert length(velocity) == D

    # single piece.

    # assemble coefficients
    x = collect( [position[d]; velocity[d]] for d = 1:D ) # starting position, velocity.
    u = collect( [velocity[d]; zero(T)] for d = 1:D ) # velocity, acceleration.

    coefficients = Vector{GeodesicPowerSeries{T}}(undef, 0)
    push!(coefficients, GeodesicPowerSeries(x,u))

    # time.
    expansion_points = [t_start;]
    steps::Vector{T} = [ (t_fin-t_start)*1.1 ;] # 1.1 instead of 1.0 so that the end point t_fin is included in the interval of validity for the solution piece.

    return PiecewiseTaylorPolynomial(
        coefficients,
        expansion_points,
        steps,
        # position,
        # velocity,
        # Vector{Vector{T}}(undef, 0),
    )
end =#

# function findbindingconstraint(
#     sol::PiecewiseTaylorPolynomial,
#     constraints::BoundConstraints
#     )
#     # TODO: do this for basic error checking on the intersection with many constraints.
#     # since the 3D plots are getting unwiedly for many constraints and 3D.
# end

# returns the first sampled time that violates the constraint. Returns `t_fin` if no constraint violations.
function forwardseekintersection(
    sol::PiecewiseTaylorPolynomial,
    constraints::ConstraintType,
    t_start,
    t_fin,
    pt_trait::IVPVariationTrait;
    N_samples = 10000,
    )

    x_t = allocatevariablecontainer(sol)
    
    t_range = LinRange(t_start, t_fin, N_samples)
    feasible_flags = trues(getNconstraints(constraints))

    # find interval for first constraint violation.
    for t in t_range
        
        evalsolution!(
            x_t,
            pt_trait,
            sol,
            t,
        )
        #@show x_t.position

        checkconstraints!(feasible_flags, x_t.position, constraints)
        if !(all(feasible_flags))
            return t, step(t_range)
        end
    end

    return t_fin, step(t_range)
end

function checkconstraints!(
    feasible_flags::BitVector, # m-th entry is true if x0 satisfies the m-th constraint in `constraints`.
    x0::Vector{T},
    constraints::AffineConstraints,
    ) where T

    N_affine = getNconstraints(constraints.hyperplane)
    
    checkconstraints!(feasible_flags, x0, constraints.hyperplane; ind_offset = 0)
    checkconstraints!(feasible_flags, x0, constraints.bound; ind_offset = N_affine)

    return nothing
end

function checkconstraints!(
    feasible_flags::BitVector,
    x0::Vector{T},
    constraints::HyperplaneConstraints;
    ind_offset::Integer = 0,
    ) where T

    as = constraints.normals
    bs = constraints.offsets

    for m in eachindex(as)
        feasible_flags[ind_offset+m] = ( dot(as[m], x0) <= bs[m] )
    end

    return nothing
end

function checkconstraints!(
    feasible_flags::BitVector,
    x0::Vector{T},
    constraints::BoundConstraints;
    ind_offset::Integer = 0,
    ) where T

    lb_dims = constraints.lb_dims
    lbs = constraints.lbs

    ub_dims = constraints.ub_dims
    ubs = constraints.ubs

    k = length(ub_dims)
    #resize!(feasible_flags, k*2)

    for m in eachindex(lb_dims)
        d = lb_dims[m]
        feasible_flags[ind_offset+m] = ( lbs[m] <= x0[d] )
    end

    for m in eachindex(ub_dims)
        d = ub_dims[m]
        feasible_flags[ind_offset+m+k] = ( x0[d] <= ubs[m] )
    end

    return nothing
end
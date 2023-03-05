
# function tautology(args...)::Bool
#     return true
# end

# function contradiction(args...)::Bool
#     return false
# end

#= # straight line. # I am here.
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

function findbindingconstraint(
    sol::PiecewiseTaylorPolynomial,
    constraints::BoundConstraints
    )
    # I am here. do this for basic error checking on the intersection with many constraints.
    # since the 3D plots are getting unwiedly for many constraints and 3D.
end
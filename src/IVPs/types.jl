# this script contain types and basic setter/getter methods.
# more complex methods for each type are in methods.jl.

###### step size strategies specific for geodesic IVPs. Encode as traits.

struct VelocityContinuityStep <: NonProbingStep end # dx/dt at t = h vs. u_next(0). h is step size. This is equivalent to the abs() of the last Taylor coefficient of u(t), since u(h) = u_next(0) by initial conditions of the next IVP piece. See notes for details.


############ solution piece of an interval from the geodesic IVP.

struct GeodesicPiece{T}
    
    # [variable index][order index]
    x::Vector{Vector{T}} # coefficients for the solution state.
    u::Vector{Vector{T}} # coefficients for first-derivative of state.
    vs::Vector{Vector{Vector{T}}} # i-th entry contain the coefficients for the i-th parallel vector field.
end

function GeodesicPiece(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    parallel_transport::Vector,
    )::GeodesicPiece{T} where T
    
    vs::Vector{Vector{Vector{T}}} = collect( parallel_transport[m].v.c for m in eachindex(parallel_transport) )

    return GeodesicPiece(x, u, vs)
end

function GeodesicPiece(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    )::GeodesicPiece{T} where T
    
    return GeodesicPiece(x, u, Vector{Vector{Vector{T}}}(undef, 0))
end

function getNvars(c::GeodesicPiece)::Int
    return length(c.x)
end

function getorder(C::GeodesicPiece)::Int
    return length(C.x[begin])-1
end

########### evaluation container.

struct GeodesicEvaluation{T}
    position::Vector{T}
    velocity::Vector{T}
    vector_fields::Vector{Vector{T}}
end

function GeodesicEvaluation(::Type{T}, N_vars::Integer, N_transports::Integer)::GeodesicEvaluation{T} where T
    return GeodesicEvaluation(
        ones(T, N_vars) .* NaN,
        ones(T, N_vars) .* NaN,
        collect( ones(T, N_vars) .* NaN for _ = 1:N_transports ),
    )
end


########### the entire solution container of a geodesic IVP.

struct PiecewiseTaylorPolynomial{T}

    # [piece index]
    coefficients::Vector{GeodesicPiece{T}}

    # [piece index]
    expansion_points::Vector{T}
    steps::Vector{T} # diagnostic information. Use this to determine if t_fin was actually reached by our solver algorithm for a given ODE solution.

    # the initial conditions.
    starting_position::Vector{T}
    starting_velocity::Vector{T}
    starting_vectors::Vector{Vector{T}} # transport vectors.
end

function getendtime(sol::PiecewiseTaylorPolynomial{T})::T where T
    return sol.expansion_points[end] + sol.steps[end]
end

function getstarttime(sol::PiecewiseTaylorPolynomial{T})::T where T
    return sol.expansion_points[begin]
end

function getNpieces(A::PiecewiseTaylorPolynomial)::Int
    return length(A.coefficients)
end

function getNvars(A::PiecewiseTaylorPolynomial{T})::Int where T
    return length(A.coefficients[begin].x)
end

function getNtransports(A::PiecewiseTaylorPolynomial{T}) where T
    return length(A.coefficients[begin].vs)
end

function getpieceorders(A::PiecewiseTaylorPolynomial)::Vector{Int}
    return collect( getorder(A.coefficients[i]) for i in eachindex(A.coefficients) )
end
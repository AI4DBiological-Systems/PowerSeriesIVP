

##### Geodesic equations could have parallel transport, at a higher computational cost. Use a trait to signify on or off.

abstract type ParallelTransportTrait <: IVPVariationTrait end

struct DisableParallelTransport <: ParallelTransportTrait end
struct EnableParallelTransport <: ParallelTransportTrait end

####  the super type that is used in specific geodesic equations' composition functions.
# e.g., used in RQ22.jl.
abstract type GeodesicIVPBuffer <: IVPBuffer end


#### parameters to the geodesic equation IVPs.

abstract type MetricParams <: IVPParameters end

struct RQ22Metric{T} <: MetricParams
    a::T
    b::T
end

# use getbuffertype(::MetricParams) = GeodesicIVPBuffer.
# i.e. getbuffertype(::RQ22Metric) = RQ22IVPBuffer is defined in RQ22.jl.


#### problem statement: parameters, and initial conditions.

struct GeodesicIVPStatement{MT,T} <: IVPStatement
    metric_params::MT
    x0::Vector{T}
    u0::Vector{T}
    v0_set::Vector{Vector{T}}
end

function GeodesicIVPStatement(
    metric_params::MT,
    x0::Vector{T},
    u0::Vector{T},
    )::GeodesicIVPStatement{MT,T} where {MT,T}

    return GeodesicIVPStatement(metric_params, x0, u0, Vector{Vector{T}}(undef, 0))
end

function getdim(A::GeodesicIVPStatement)
    return length(A.x0)
end

function getNtransports(A::GeodesicIVPStatement)
    return length(A.v0_set)
end

getsoltype(::GeodesicIVPStatement{MT,T}) where {MT,T} = PiecewiseTaylorPolynomial{T,GeodesicPiece{T}}



# this script contain types and basic setter/getter methods.
# more complex methods for each type are in methods.jl.

###### step size strategies specific for geodesic IVPs. Encode as traits.

struct VelocityContinuityStep <: NonProbingStep end # dx/dt at t = h vs. u_next(0). h is step size. This is equivalent to the abs() of the last Taylor coefficient of u(t), since u(h) = u_next(0) by initial conditions of the next IVP piece. See notes for details.


############ solution piece of an interval from the geodesic IVP.

struct GeodesicPiece{T} <: SolutionPiece
    
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

getsolpiecetype(::GeodesicIVPStatement{MT,T}) where {MT,T} = GeodesicPiece{T}

function getdim(C::GeodesicPiece)::Int
    return length(C.x)
end

function getorder(C::GeodesicPiece)::Int
    return length(C.x[begin])-1
end

function getNtransports(C::GeodesicPiece)::Int
    return length(C.vs)
end

########### evaluation container.

struct GeodesicEvaluation{T} <: VariableContainer
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

function getvariablecontainer(
    prob_params::GeodesicIVPStatement{MT,T},
    )::GeodesicEvaluation{T} where {T, MT}
    
    return GeodesicEvaluation(
        T,
        getdim(prob_params),
        getNtransports(prob_params),
    )
end
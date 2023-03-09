

##### Geodesic equations could have parallel transport, at a higher computational cost. Use a trait to signify on or off.

abstract type ParallelTransportTrait <: IVPVariationTrait end

struct DisableParallelTransport <: ParallelTransportTrait end
struct EnableParallelTransport <: ParallelTransportTrait end

struct EnableGeodesicLine{PT<: ParallelTransportTrait} <: IVPVariationTrait
    pt_trait::PT
end

####  the super type that is used in specific geodesic equations' composition functions.
# e.g., used in RQ22.jl.
abstract type GeodesicIVPBuffer <: IVPBuffer end

abstract type GeodesicSolutionPiece <: SolutionPiece end

#### parameters to the geodesic equation IVPs.

abstract type MetricParams <: IVPParameters end
abstract type NonEuclideanMetric <: MetricParams end

struct RQ22Metric{T} <: NonEuclideanMetric
    a::T
    b::T
end

struct EuclideanMetric <: MetricParams end


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




# this script contain types and basic setter/getter methods.
# more complex methods for each type are in methods.jl.

###### step size strategies specific for geodesic IVPs. Encode as traits.

struct VelocityContinuityStep <: NonProbingStep end # dx/dt at t = h vs. u_next(0). h is step size. This is equivalent to the abs() of the last Taylor coefficient of u(t), since u(h) = u_next(0) by initial conditions of the next IVP piece. See notes for details.


############ solution piece of an interval from the geodesic IVP.

struct GeodesicPowerSeries{T} <: GeodesicSolutionPiece
    
    # [dimension index][order index]
    x::Vector{Vector{T}} # coefficients for the position.
    u::Vector{Vector{T}} # coefficients for first-derivative of position.
    vs::Vector{Vector{Vector{T}}} # [vec field][dim][order] i-th entry contain the coefficients for the i-th parallel vector field.
end

function GeodesicPowerSeries(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    parallel_transport::Vector,
    )::GeodesicPowerSeries{T} where T
    
    vs::Vector{Vector{Vector{T}}} = collect( parallel_transport[m].v.c for m in eachindex(parallel_transport) )

    return GeodesicPowerSeries(x, u, vs)
end

function GeodesicPowerSeries(
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    )::GeodesicPowerSeries{T} where T
    
    return GeodesicPowerSeries(x, u, Vector{Vector{Vector{T}}}(undef, 0))
end

getsoltype(::GeodesicIVPStatement{MT,T}) where {MT,T} = PiecewiseTaylorPolynomial{T,GeodesicPowerSeries{T}}
getsolpiecetype(::GeodesicIVPStatement{MT,T}) where {MT,T} = GeodesicPowerSeries{T}

function getdim(C::GeodesicSolutionPiece)::Int
    return length(C.x)
end

function getNtransports(C::GeodesicSolutionPiece)::Int
    return length(C.vs)
end

function getorder(C::GeodesicPowerSeries)::Int
    return length(C.x[begin])-1
end

########## the line case.
# see line_methods.jl in the IVP/geodesic_eqns folder for methods.

# GeodesicLine replaces GeodesicPowerSeries when we solve the geodesic IVP with the Euclidean Riemannian metric.
struct GeodesicLine{T} <: GeodesicSolutionPiece
    
    # [dimension index]
    x::Vector{T} # position.
    u::Vector{T} # velocity.
    vs::Vector{Vector{T}} # [m-th vector field][d-th dimension] set of vector fields.
end

getsoltype(::GeodesicIVPStatement{EuclideanMetric,T}) where T = PiecewiseTaylorPolynomial{T,GeodesicLine{T}}
getsolpiecetype(::GeodesicIVPStatement{EuclideanMetric,T}) where T = GeodesicLine{T}
getvariabletype(::PiecewiseTaylorPolynomial{T,GeodesicLine{T}}) where T = GeodesicEvaluation{T}


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

getvariabletype(::GeodesicIVPStatement{MT,T}) where {T,MT} = GeodesicEvaluation{T}
getvariabletype(::PiecewiseTaylorPolynomial{T,GeodesicPowerSeries{T}}) where T = GeodesicEvaluation{T}

function allocatevariablecontainer(
    p::GeodesicIVPStatement{MT,T},
    )::GeodesicEvaluation{T} where {T, MT}
    
    return GeodesicEvaluation(T, getdim(p), getNtransports(p))
end

function allocatevariablecontainer(
    A::PiecewiseTaylorPolynomial{T,SP},
    )::GeodesicEvaluation{T} where {T,SP<:GeodesicSolutionPiece}
    
    return GeodesicEvaluation(T, getdim(A), getNtransports(A))
end


###### other traits.

getIVPtrait(::GeodesicIVPStatement{MT,T}) where {MT,T} = NumericalIVPTrait() # default, catch-all.
getIVPtrait(::GeodesicIVPStatement{EuclideanMetric,T}) where T = LineIVPTrait() # an analytical case: solve without invoking numerical methods for the ODE part.


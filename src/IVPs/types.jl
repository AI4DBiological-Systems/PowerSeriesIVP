
##### for solveIVP!()

abstract type MetricParams end

struct RQ22Metric{T} <: MetricParams
    a::T
    b::T
end

struct GeodesicIVPProblem{MT,T}
    metric_params::MT
    x0::Vector{T}
    u0::Vector{T}
    v0_set::Vector{Vector{T}}
end

function GeodesicIVPProblem(
    metric_params::MT,
    x0::Vector{T},
    u0::Vector{T},
    )::GeodesicIVPProblem{MT,T} where {MT,T}

    return GeodesicIVPProblem(metric_params, x0, u0, Vector{Vector{T}}(undef, 0))
end

function getNvars(A::GeodesicIVPProblem)
    return length(A.x0)
end

function getNtransports(A::GeodesicIVPProblem)
    return length(A.v0_set)
end

abstract type ParallelTransportTrait end

# struct ParallelTransport{T} <: ParallelTransportTrait
#     v0_set::Vector{Vector{T}} # v at t=0.
#     v_set::Vector{IntegralSequence{T}}
# end

struct DisableParallelTransport <: ParallelTransportTrait end
struct EnableParallelTransport <: ParallelTransportTrait end

#####

abstract type GeodesicIVPBuffer end
# sutypes of GeodesicIVPBuffer are defined in RQ22.jl

##### configs.

abstract type IVPConfig end

# still adapts step.
struct FixedOrderConfig{T} <: IVPConfig
    ϵ::T
    h_zero_error::T
    L::Int

    # IVP-related
    step_reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    max_pieces::Int
end

function FixedOrderConfig(
    ::Type{T};
    ϵ::Real = 1e-6,
    h_zero_error::Real = Inf,
    L::Integer = 10,
    max_pieces::Integer = typemax(Int),
    step_reduction_factor = 2,
    ) where T

    return FixedOrderConfig(
        convert(T, ϵ),
        convert(T, h_zero_error),
        convert(Int, L),
        convert(T, step_reduction_factor),
        convert(Int, max_pieces),
    )
end

# adapt order and step size.
struct AdaptOrderConfig{T} <: IVPConfig
    ϵ::T
    L_min::Int
    L_max::Int # L_test_max
    r_order::T
    h_zero_error::T
    N_analysis_terms::Int

    # IVP-related
    step_reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    max_pieces::Int

    # buffers
    Es::Vector{T} # errors, length N_analysis_terms.
end

function AdaptOrderConfig(
    ::Type{T};
    ϵ::T = convert(T, 1e-6),
    L_min::Int = 4,
    L_max::Int = convert(Int, 10),
    r_order = convert(T, 0.3),
    h_zero_error = convert(T, Inf),
    N_analysis_terms::Int = convert(Int, 2),
    max_pieces::Integer = typemax(Int),
    step_reduction_factor = 2,
    ) where T

    @assert L_min > N_analysis_terms
    @assert L_min < L_max
    @assert N_analysis_terms > 0
    @assert max_pieces > 0
    @assert step_reduction_factor > 0

    return AdaptOrderConfig(
        convert(T, ϵ),
        convert(Int, L_min),
        convert(Int, L_max),
        convert(T, r_order),
        convert(T, h_zero_error),
        convert(Int, N_analysis_terms),
        convert(T, step_reduction_factor),
        convert(Int, max_pieces),
    )
end



############# constraints

struct PositiveRealTrait end

struct AffineConstraints{T}
    normals::Vector{Vector{T}}
    offsets::Vector{T}
end

struct IntersectionBuffer{T<:AbstractFloat}
    intersection_coefficients::Vector{Vector{T}} # [constraints][order]
    all_roots::Vector{Complex{T}} # [order]
    smallest_positive_roots::Vector{T} # [constraints], real roots.
    zero_tol::T # for deciding wheather a complex number variable is a real number.
end

function IntersectionBuffer(zero_tol::T, order::Integer, N_constraints::Integer)::IntersectionBuffer{T} where T
    return IntersectionBuffer(
        collect( zeros(T, order) for _ = 1:N_constraints ),
        Vector{Complex{T}}(undef, order),
        Vector{T}(undef, N_constraints),
        zero_tol,
    )
end

struct BudanIntersectionBuffers{T}
    # [constraints][order]
    cs::Vector{Vector{T}}
    #cs_left::Vector{Vector{T}}
    cs_right::Vector{Vector{T}}
    #cs_center::Vector{Vector{T}}
end

function BudanIntersectionBuffers(::Type{T}, N_constraints::Integer, order::Integer)::BudanIntersectionBuffers{T} where T
    
    return BudanIntersectionBuffers(
        collect( zeros(T, order) for _ = 1:N_constraints ),
        #collect( zeros(T, order) for _ = 1:N_constraints ),
        collect( zeros(T, order) for _ = 1:N_constraints ),
        #collect( zeros(T, order) for _ = 1:N_constraints ),
    )
end

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

##### adaptive step
abstract type StepStrategyTrait end

struct VelocityContinuityStep <: StepStrategyTrait end # for geodesic problems.
struct GuentherWolfStep <: StepStrategyTrait end # for any ODE IVP.
#struct MinAllStep <: StepStrategyTrait end # for geodesic problems.

##### configs.

abstract type IVPConfig end

# still adapts step.
struct FixedOrderConfig{T} <: IVPConfig
    ϵ::T
    h_max::T
    L::Int

    # IVP-related
    step_reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    max_pieces::Int
end

function FixedOrderConfig(
    ::Type{T};
    ϵ::Real = 1e-6,
    h_max::Real = Inf,
    L::Integer = 10,
    max_pieces::Integer = typemax(Int),
    step_reduction_factor = 2,
    ) where T

    return FixedOrderConfig(
        convert(T, ϵ),
        convert(T, h_max),
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
    h_max::T
    N_analysis_terms::Int

    # IVP-related
    step_reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    max_pieces::Int

    # buffers
    #Es::Vector{T} # errors, length N_analysis_terms.
    #explicit_roots_buffer::IntersectionBuffer{T}
end

function AdaptOrderConfig(
    ::Type{T};
    ϵ::T = convert(T, 1e-6),
    L_min::Int = 4,
    L_max::Int = convert(Int, 10),
    r_order = convert(T, 0.3),
    h_max = convert(T, 1),
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
        convert(T, h_max),
        convert(Int, N_analysis_terms),
        convert(T, step_reduction_factor),
        convert(Int, max_pieces),
    )
end

############ continuity step refinement.

# abstract type ContinuityConfig end
# # 1-st order continuity conditions for x and u.
# struct ContinuityConfigXU{T} <: ContinuityConfig
#     stop_order::Int # >= 0.
#     zero_tol::T

#     # buffer.
#     discrepancies::Matrix{T} # [order coefficient index][variable dimension]
# end
struct ContinuityConfig{T}
    zero_tol::T # above 0.
    min_h::T # above 0.
    discount_factor::T # in the interval (0,1).
end

############# constraints step refinement.

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
    cs_right::Vector{Vector{T}}
    
    ubs::Vector{T} # buffer, [constraints].

    #cs_left::Vector{Vector{T}}
    #cs_center::Vector{Vector{T}}
end

function BudanIntersectionBuffers(::Type{T}, N_constraints::Integer, order::Integer)::BudanIntersectionBuffers{T} where T
    
    return BudanIntersectionBuffers(
        collect( zeros(T, order) for _ = 1:N_constraints ),
        collect( zeros(T, order) for _ = 1:N_constraints ),
        ones(T, N_constraints) .* NaN,
    )
end


#### numerical solver types (put into separate library later.)
struct ITPConfig{T}
    f_tol::T
    x_tol::T
    k1::T
    k2::T
    n0::Int
end

function ITPConfig(
    ::Type{T};
    f_tol::T = convert(T, 1e-8),
    x_tol::T = convert(T, 1e-15),
    k1::T = convert(T, 0.1),
    k2::T = convert(T, 0.98*(1+MathConstants.golden)), # see equation 24.
    n0::Int = convert(Int, 0),
    )::ITPConfig{T} where T

    @assert k1 > zero(T)
    @assert one(T) <= k2 < one(T) + MathConstants.golden
    @assert n0 >= 0
    @assert x_tol > 0
    @assert f_tol > 0

    return ITPConfig(f_tol, x_tol, k1, k2, n0)
end

#### constraints front end.

abstract type ConstraintType end

struct NoConstraints <: ConstraintType end


struct AffineConstraintsContainer{T} <: ConstraintType
    #
    constraints::AffineConstraints{T}

    # buffers
    explicit_roots_buffer::IntersectionBuffer{T}
    upperbound_buffer::BudanIntersectionBuffers{T}
    bino_mat::Matrix{Int}
    max_divisions::Int

    # configs
    solver_config::ITPConfig{T}
end
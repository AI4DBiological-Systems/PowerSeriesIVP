
##### PSM elementary operation super types and traits.
abstract type SingleVariableTrait end # for PSM elementary operation containers that have a single field of type Vector{Vector{T<:Real}}.

function resetbuffer!(A::SingleVariableTrait)    
    resetcoefficients!(A.c)

    return nothing
end

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



#####


# sutypes of GeodesicIVPBuffer are defined in RQ22.jl

##### adaptive step
abstract type StepStrategyTrait end

### subtypes of NonProbingStep are step strategies that only use the current solution piece. It does not probe for a next piece. They are faster to compute, but less error guarantees, as they are based on Taylor coefficients. The coefficients might be zero, which poses a problem.

abstract type NonProbingStep <: StepStrategyTrait end

struct GuentherWolfStep <: NonProbingStep end # for any ODE IVP, using the last Taylor coefficient of the primary variable.

abstract type AllNonProbingStep <: StepStrategyTrait end # try all the strategies under NonProbingStep, and pick the smallest step.

### these are more stringent than GuentherWolfStep for the same level of tolerance.
# the higher the derivative order, the more stringent (i.e. yields smaller step) the condition imposes.

abstract type DerivativeContinuityStep <: StepStrategyTrait end

struct ContinuityFirstDerivative{ST<:NonProbingStep} <: DerivativeContinuityStep
    first_step_strategy::ST
end

# condition based on the first two derivatives.
struct ContinuitySecondDerivative{ST<:NonProbingStep} <: DerivativeContinuityStep
    first_step_strategy::ST
end

# conditions based on the first max_order derivatives.
struct ContinuityHigherDerivative{ST<:NonProbingStep} <: DerivativeContinuityStep 
    first_step_strategy::ST
    max_order::Int
end


##### configs.

struct StepConfig{T,ST<:StepStrategyTrait}
    ϵ::T
    h_max::T
    reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    discount_factor::T
    strategy::ST
end

function StepConfig(
    ::Type{T},
    strategy::ST;
    ϵ::Real = 1e-6,
    h_max::Real = Inf,
    reduction_factor = 2,
    discount_factor = 0.9,
    )::StepConfig{T,ST} where {T,ST<:StepStrategyTrait}

    return StepConfig(
        convert(T, ϵ),
        convert(T, h_max),
        convert(T, reduction_factor),
        convert(T, discount_factor),
        strategy,
    )
end

abstract type IVPConfig end

# still adapts step.
struct FixedOrderConfig{T,ST<:StepStrategyTrait} <: IVPConfig

    L::Int # the fixed order.
    max_pieces::Int
    
    #stopping_conditions::StopConditions{T}
    min_step::T
    step_config::StepConfig{T,ST}
end

function FixedOrderConfig(
    step_config::StepConfig{T,ST};
    L = 4,
    max_pieces = 100000,
    min_step = eps(T)*10,
    )::FixedOrderConfig{T,ST} where {T,ST<:StepStrategyTrait}

    return FixedOrderConfig(
        L,
        max_pieces,
        convert(T, min_step),
        step_config,
    )
end

# adapt order and step size.
struct AdaptOrderConfig{T,ST<:StepStrategyTrait} <: IVPConfig
    
    L_min::Int
    L_max::Int # L_test_max
    order_increase_factor::T
    # r_order::T
    # N_analysis_terms::Int

    # IVP-related
    max_pieces::Int

    #stopping_conditions::StopConditions{T}
    min_step::T
    step_config::StepConfig{T,ST}
end

function AdaptOrderConfig(
    step_config::StepConfig{T,ST};
    L_min = 4,
    L_max = 13,
    order_increase_factor = 1.35,
    max_pieces = 100000,
    min_step = eps(T)*10,
    )::AdaptOrderConfig{T,ST} where {T,ST<:StepStrategyTrait}

    return AdaptOrderConfig(
        L_min,
        L_max, 
        convert(T, order_increase_factor),
        max_pieces,
        convert(T, min_step),
        step_config,
    )
end


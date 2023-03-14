
##### PSM elementary operation super types and traits.
abstract type SingleVariableTrait end # for PSM elementary operation containers that have a single field of type Vector{Vector{T<:Real}}.

function resetbuffer!(A::SingleVariableTrait)    
    resetcoefficients!(A.c)

    return nothing
end

#### trait functions.
# scattered in various source files.
# getsoltype()
# getbuffertype()

##### IVP

# ## data types for the IVP definition.

abstract type IVPParameters end

abstract type IVPBuffer end

abstract type IVPStatement end



# to toggle between similar IVPs that share common computation features. Will be used in the future if we handle non-geodesic equations IVPs.
abstract type IVPVariationTrait end

# to decide whether to use the power series method to numerically solve the IVP, or use an analytic solution.
# search getIVPtrait() for the trait function.
abstract type IVPTrait end
struct NumericalIVPTrait <: IVPTrait end # the IVP requires the power series method, a numerical solution.
struct LineIVPTrait <: IVPTrait end # the IVP has an analytic solution: a line.

##### variable types, which defines the IVP.

abstract type VariableContainer end

########### the entire solution container of a geodesic IVP.

abstract type SolutionPiece end

struct PiecewiseTaylorPolynomial{T, VT <: SolutionPiece}

    # [piece index]
    coefficients::Vector{VT}

    # [piece index]
    expansion_points::Vector{T}
    steps::Vector{T} # diagnostic information. Use this to determine if t_fin was actually reached by our solver algorithm for a given ODE solution.

    # # the initial conditions.
    # starting_position::Vector{T}
    # starting_velocity::Vector{T}
    # starting_vectors::Vector{Vector{T}} # transport vectors.
    #initial_conditions::Vector{VT}
end

# use the trait function getsoltype(::IVPStatement) to get the sxact PiecewiseTaylorPolynomial{T,VT} data type. See geodesic_types.jl.

function PiecewiseTaylorPolynomial(::Type{T}, ::Type{VT}) where {T, VT<: SolutionPiece}
    return PiecewiseTaylorPolynomial(
        Vector{VT}(undef, 0),
        Vector{T}(undef, 0),
        Vector{T}(undef, 0),
    )
end

function getendtime(sol::PiecewiseTaylorPolynomial{T,VT})::T where {T,VT}
    return sol.expansion_points[end] + sol.steps[end]
end

function getstarttime(sol::PiecewiseTaylorPolynomial{T,VT})::T where {T,VT}
    return sol.expansion_points[begin]
end

function getNpieces(A::PiecewiseTaylorPolynomial)::Int
    return length(A.coefficients)
end

function getpieceorders(A::PiecewiseTaylorPolynomial)::Vector{Int}
    return collect( getorder(A.coefficients[i]) for i in eachindex(A.coefficients) )
end

function getdim(A::PiecewiseTaylorPolynomial)::Int
    return getdim(A.coefficients[begin])
end

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
    h_default::T
    reduction_factor::T # a step reduction factor to be applied in addition to the step size formula. larger or equal to 1.
    discount_factor::T
    strategy::ST
end

function StepConfig(
    ::Type{T},
    strategy::ST;
    ϵ::Real = 1e-6,
    h_default::Real = Inf,
    reduction_factor = 2,
    discount_factor = 0.9,
    )::StepConfig{T,ST} where {T,ST<:StepStrategyTrait}

    return StepConfig(
        convert(T, ϵ),
        convert(T, h_default),
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



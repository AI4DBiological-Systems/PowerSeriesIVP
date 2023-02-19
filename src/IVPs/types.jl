
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

#####

struct RQ22ζ{T}

    # stage 1
    # delta::InterVariableDifference{T}
    ## W4::ΔSumCol{T}
    #Z8::ΔSumCol{T}
    ## A::AddConstant{T}
    # B::AddConstant{T}
    # R::ScaledReciprocal{T}

    # stage 2
    σ::Quotient{T} # v/B
    Z0::Product{T} # η*σ
    Z9::ΔSumColProduct{T}

    # stage 3
    Z3::Product{T}
    Z5::ΔSumColProduct{T}
    Z1::SubtractFromScaled{T}

    # stage 4
    Z6::Product{T}
    Z7::Product{T}
    Z2::Addition{T}

    # du/dt.
    ζ::Product{T}
    c::Vector{Vector{T}} # this is ζ.c, the coefficients for ζ.

    # constants
    a_sq::T
    b_sq::T
    a_sq_m_b_sq::T
end

function RQ22ζ(a::T, b::T, N::Integer)::RQ22ζ{T} where T
    @assert a > zero(T)
    @assert b > zero(T)

    ζ = Product(T,N)

    return RQ22ζ(
        #ΔSumCol(T,N),
        
        # stage 2
        Quotient(T,N),
        Product(T,N),
        ΔSumColProduct(T,N),

        # stage 3
        Product(T,N),
        ΔSumColProduct(T,N),
        SubtractFromScaled(T,N),

        # stage 4
        Product(T,N),
        Product(T,N),
        Addition(T,N),

        # du/dt.
        ζ, ζ.c,

        # constants.
        a^2, b^2, a^2 - b^2,
    )
end

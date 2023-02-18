# θ := du/dt.
struct RQ22θ{T}

    # stage 1
    delta::InterVariableDifference{T}
    W4::ΔSumCol{T}
    W8::ΔSumCol{T}
    A::AddConstant{T}
    B::AddConstant{T}
    R::ScaledReciprocal{T}

    # stage 2
    η::Quotient{T}
    C::Squared{T}
    W9::ΔSumColProduct{T}

    # stage 3
    W3::Product{T}
    W5::ΔSumColProduct{T}
    W1::SubtractFromScaled{T}

    # stage 4
    W6::Product{T}
    W7::Product{T}
    W2::Addition{T}

    # du/dt.
    θ::Product{T}
    c::Vector{Vector{T}} # this is θ.c, the coefficients for ζ.

    # constants
    a_sq::T
    b_sq::T
    a_sq_m_b_sq::T
end

function RQ22θ(a::T, b::T, N::Integer)::RQ22θ{T} where T
    @assert a > zero(T)
    @assert b > zero(T)

    θ = Product(T,N)

    return RQ22θ(
        InterVariableDifference(T,N),
        ΔSumCol(T,N),
        ΔSumCol(T,N),
        AddConstant(T,N),
        AddConstant(T,N),
        ScaledReciprocal(T,N),

        # stage 2
        Quotient(T,N),
        Squared(T,N),
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
        θ, θ.c,

        # constants.
        a^2, b^2, a^2 - b^2,
    )
end

function initializeorder!(
    p::RQ22θ{T},
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    ) where T
    
    @assert length(x) == length(u)

    # stage 1
    initializeorder!(p.delta, x)
    initializeorder!(p.W4, p.delta.Δ_sq)
    initializeorder!(p.W8, p.delta.Δ)
    initializeorder!(p.A, p.W4.c, p.a_sq)
    initializeorder!(p.B, p.W4.c, p.b_sq)
    initializeorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    initializeorder!(p.η, u, p.B.c)
    initializeorder!(p.C, p.η.c)
    initializeorder!(p.W9, p.delta.Δ, p.C.c) # ΔSumColProduct

    # stage 3
    initializeorder!(p.W3, p.W8.c, u)
    initializeorder!(p.W5, p.delta.Δ, u)
    initializeorder!(p.W1, p.W5.c, p.W3.c, convert(T, 2))

    # stage 4
    initializeorder!(p.W6, p.η.c, p.W1.c)
    initializeorder!(p.W7, p.B.c, p.W9.c)
    initializeorder!(p.W2, p.W6.c, p.W7.c)

    # du/dt.
    initializeorder!(p.θ, p.R.c, p.W2.c) # Product

    return nothing
end

function increaseorder!(
    p::RQ22θ{T},
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    ) where T

    # stage 1
    increaseorder!(p.delta, x)
    increaseorder!(p.W4, p.delta.Δ_sq)
    increaseorder!(p.W8, p.delta.Δ)
    increaseorder!(p.A, p.W4.c, p.a_sq)
    increaseorder!(p.B, p.W4.c, p.b_sq)
    increaseorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    increaseorder!(p.η, u, p.B.c)
    increaseorder!(p.C, p.η.c)
    increaseorder!(p.W9, p.delta.Δ, p.C.c) # ΔSumColProduct

    # stage 3
    increaseorder!(p.W3, p.W8.c, u)
    increaseorder!(p.W5, p.delta.Δ, u)
    increaseorder!(p.W1, p.W5.c, p.W3.c, convert(T, 2))

    # stage 4
    increaseorder!(p.W6, p.η.c, p.W1.c)
    increaseorder!(p.W7, p.B.c, p.W9.c)
    increaseorder!(p.W2, p.W6.c, p.W7.c)

    # RHS of du/dt.
    increaseorder!(p.θ, p.R.c, p.W2.c)

    return nothing
end



#################### geodesic IVP problem.

struct RQ22ParallelTransport{T}
    ζ::Vector{RQ22ζ{T}}
    v_set::Vector{IntegralSequence{T}}
    v0_set::Vector{Vector{T}}
end

struct RQ22IVPBuffer{T} <: GeodesicIVPBuffer # formally RQGeodesicBuffer{T}
    
    # Its contents used only in parallel transport.
    #parallel_transport::PT
    parallel_transport::Vector{RQ22ParallelTransport{T}} # used only in parallel transport.

    # RHS of du/dt.
    θ::RQ22θ{T}

    # variables
    u::IntegralSequence{T} # dx/dt.
    x::IntegralSequence{T}

    # initial conditions
    x0::Vector{T} # x at t = 0.
    u0::Vector{T} # dx/dt at t = 0.

end

function getivpbuffer(
    #pt::PT,
    metric_params::RQ22Metric{T},
    x0::Vector{T},
    u0::Vector{T},
    )::RQ22IVPBuffer{T} where T

    a = metric_params.a
    b = metric_params.b

    @assert length(x0) == length(u0)
    @assert a > zero(T)
    @assert b > zero(T)
    N = length(x0)

    return RQ22IVPBuffer(
        #pt,
        Vector{RQ22ParallelTransport{T}}(undef,0),

        # RHS of du/dt.
        RQ22θ(a,b,N),

        # variables
        IntegralSequence(T,N),
        IntegralSequence(T,N),
        
        # initial conditions
        x0, u0,
    )
end


################### parallel transport

# I am here.
function initializeorder!(
    p::RQ22ζ{T},
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    ) where T
    
    @assert length(x) == length(u)

    # stage 1
    initializeorder!(p.delta, x)
    initializeorder!(p.W4, p.delta.Δ_sq)
    initializeorder!(p.W8, p.delta.Δ)
    initializeorder!(p.A, p.W4.c, p.a_sq)
    initializeorder!(p.B, p.W4.c, p.b_sq)
    initializeorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    initializeorder!(p.η, u, p.B.c)
    initializeorder!(p.C, p.η.c)
    initializeorder!(p.W9, p.delta.Δ, p.C.c) # ΔSumColProduct

    # stage 3
    initializeorder!(p.W3, p.W8.c, u)
    initializeorder!(p.W5, p.delta.Δ, u)
    initializeorder!(p.W1, p.W5.c, p.W3.c, convert(T, 2))

    # stage 4
    initializeorder!(p.W6, p.η.c, p.W1.c)
    initializeorder!(p.W7, p.B.c, p.W9.c)
    initializeorder!(p.W2, p.W6.c, p.W7.c)

    # du/dt.
    initializeorder!(p.θ, p.R.c, p.W2.c) # Product

    return nothing
end

function increaseorder!(
    p::RQ22ζ{T},
    x::Vector{Vector{T}},
    u::Vector{Vector{T}},
    ) where T

    # stage 1
    increaseorder!(p.delta, x)
    increaseorder!(p.W4, p.delta.Δ_sq)
    increaseorder!(p.W8, p.delta.Δ)
    increaseorder!(p.A, p.W4.c, p.a_sq)
    increaseorder!(p.B, p.W4.c, p.b_sq)
    increaseorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    increaseorder!(p.η, u, p.B.c)
    increaseorder!(p.C, p.η.c)
    increaseorder!(p.W9, p.delta.Δ, p.C.c) # ΔSumColProduct

    # stage 3
    increaseorder!(p.W3, p.W8.c, u)
    increaseorder!(p.W5, p.delta.Δ, u)
    increaseorder!(p.W1, p.W5.c, p.W3.c, convert(T, 2))

    # stage 4
    increaseorder!(p.W6, p.η.c, p.W1.c)
    increaseorder!(p.W7, p.B.c, p.W9.c)
    increaseorder!(p.W2, p.W6.c, p.W7.c)

    # RHS of du/dt.
    increaseorder!(p.θ, p.R.c, p.W2.c)

    return nothing
end
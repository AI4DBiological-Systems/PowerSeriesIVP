



### stage 1
struct RQGeodesicBuffer{T}

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

    # RHS of du/dt.
    θ::Product{T}

    # variables
    u::IntegralSequence{T} # dx/dt.
    x::IntegralSequence{T}

    # initial conditions
    x0::Vector{T} # x at t = 0.
    u0::Vector{T} # dx/dt at t = 0.

    # constants
    a_sq::T
    b_sq::T
    a_sq_m_b_sq::T
end

function RQGeodesicBuffer(a::T, b::T, x0::Vector{T}, u0::Vector{T})::RQGeodesicBuffer{T} where T
    @assert length(x0) == length(u0)
    @assert a > zero(T)
    @assert b > zero(T)
    N = length(x0)

    return RQGeodesicBuffer(
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

        # RHS of du/dt.
        Product(T,N),

        # variables
        IntegralSequence(T,N),
        IntegralSequence(T,N),

        # initial conditions and constants.
        x0, u0, a^2, b^2, a^-b^2,
    )
end

function initializeorder!(p::RQGeodesicBuffer{T}) where T

    # stage 1
    initializeorder!(p.delta, p.x.c)
    initializeorder!(p.W4, p.delta.Δ_sq)
    initializeorder!(p.W8, p.delta.Δ)
    initializeorder!(p.A, p.a_sq)
    initializeorder!(p.B, p.b_sq)
    initializeorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    initializeorder!(p.η, p.u.c, p.B.c)
    initializeorder!(p.C, p.η.c)
    initializeorder!(p.W9, p.delta.Δ, p.C.c)

    # stage 3
    initializeorder!(p.W3, p.W8.c, p.u.c)
    initializeorder!(p.W5, p.delta.Δ, p.u_c)
    initializeorder!(p.W1, p.W3.c, p.W5.c, convert(T, 2))

    # stage 4
    initializeorder!(p.W6, p.η.c, p.W1.c)
    initializeorder!(p.W7, p.B.c, p.W9.c)
    initializeorder!(p.W2, p.W6.c, p.W7.c)

    # RHS of du/dt.
    initializeorder!(p.θ, p.R.c, p.W2.c)

    # variables
    initializeorder!(p.u, p.u0)
    initializeorder!(p.x, p.x0)

    return nothing
end

function increaseorder!(p::RQGeodesicBuffer{T}) where T

    # stage 1
    increaseorder!(p.delta, p.x.c)
    increaseorder!(p.W4, p.delta.Δ_sq)
    increaseorder!(p.W8, p.delta.Δ)
    increaseorder!(p.A, p.a_sq)
    increaseorder!(p.B, p.b_sq)
    increaseorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    increaseorder!(p.η, p.u.c, p.B.c)
    increaseorder!(p.C, p.η.c)
    increaseorder!(p.W9, p.delta.Δ, p.C.c)

    # stage 3
    increaseorder!(p.W3, p.W8.c, p.u.c)
    increaseorder!(p.W5, p.delta.Δ, p.u_c)
    increaseorder!(p.W1, p.W3.c, p.W5.c, convert(T, 2))

    # stage 4
    increaseorder!(p.W6, p.η.c, p.W1.c)
    increaseorder!(p.W7, p.B.c, p.W9.c)
    increaseorder!(p.W2, p.W6.c, p.W7.c)

    # RHS of du/dt.
    increaseorder!(p.θ, p.R.c, p.W2.c)

    # variables
    increaseorder!(p.u, p.θ.c)
    increaseorder!(p.x, p.u.c)

    return nothing
end

# # parse
# delta = ivp.delta
# W4 = ivp.W4
# W8 = ivp.W8
# B = ivp.B
# A = ivp.A
# R = ivp.R

# # stage 2
# η = ivp.η
# C = ivp.C
# W9 = ivp.W9

# # stage 3
# W3 = ivp.W3
# W5 = ivp.W5
# W1 = ivp.W1

# # stage 4
# W6 = ivp.W6
# W7 = ivp.W7
# W2 = ivp.W2

# # variables
# u.c = ivp.u.c
# x.c = ivp.x.c

# # initial conditions
# x0 = ivp.x0
# u0 = ivp.u0

# # constants
# a = ivp.a
# b = ivp.b
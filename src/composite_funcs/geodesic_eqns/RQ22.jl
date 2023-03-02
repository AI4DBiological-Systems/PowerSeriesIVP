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

# does not reset A.a_sq, A.b_sq, A.a_sq_m_b_sq, which are constants to the IVP.
function resetbuffer!(p::RQ22θ)
    # stage 1.
    resetbuffer!(p.delta)
    resetbuffer!(p.W4)
    resetbuffer!(p.W8)
    resetbuffer!(p.A)
    resetbuffer!(p.B)
    resetbuffer!(p.R)

    # stage 2
    resetbuffer!(p.η)
    resetbuffer!(p.C)
    resetbuffer!(p.W9)

    # stage 3
    resetbuffer!(p.W3)
    resetbuffer!(p.W5)
    resetbuffer!(p.W1)

    # stage 4
    resetbuffer!(p.W6)
    resetbuffer!(p.W7)
    resetbuffer!(p.W2)

    # du/dt
    resetbuffer!(p.θ)
    # A.c is A.θ.c, and should be reset by the previous line.

    return nothing
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


################### parallel transport

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

function initializeorder!(
    ζ::RQ22ζ{T},
    θ::RQ22θ{T},
    v::Vector{Vector{T}},
    ) where T
    
    #@assert length(v) == length(u)

    # # stage 1
    # initializeorder!(p.delta, x)
    # initializeorder!(p.W4, p.delta.Δ_sq)
    # initializeorder!(p.W8, p.delta.Δ)
    # initializeorder!(p.A, p.W4.c, p.a_sq)
    # initializeorder!(p.B, p.W4.c, p.b_sq)
    # initializeorder!(p.R, p.A.c, p.a_sq_m_b_sq)

    # stage 2
    initializeorder!(ζ.σ, v, θ.B.c)
    initializeorder!(ζ.Z0, ζ.σ.c, θ.η.c)
    initializeorder!(ζ.Z9, θ.delta.Δ, ζ.Z0.c) # ΔSumColProduct

    # stage 3
    initializeorder!(ζ.Z3, θ.W8.c, v)
    initializeorder!(ζ.Z5, θ.delta.Δ, v)
    initializeorder!(ζ.Z1, ζ.Z5.c, ζ.Z3.c, convert(T, 2))

    # stage 4
    initializeorder!(ζ.Z6, θ.η.c, ζ.Z1.c)
    initializeorder!(ζ.Z7, θ.B.c, ζ.Z9.c)
    initializeorder!(ζ.Z2, ζ.Z6.c, ζ.Z7.c)

    # du/dt.
    initializeorder!(ζ.ζ, θ.R.c, ζ.Z2.c) # Product

    return nothing
end


function increaseorder!(
    ζ::RQ22ζ{T},
    θ::RQ22θ{T},
    v::Vector{Vector{T}},
    ) where T

    # stage 2
    increaseorder!(ζ.σ, v, θ.B.c)
    increaseorder!(ζ.Z0, ζ.σ.c, θ.η.c)
    increaseorder!(ζ.Z9, θ.delta.Δ, ζ.Z0.c) # ΔSumColProduct

    # stage 3
    increaseorder!(ζ.Z3, θ.W8.c, v)
    increaseorder!(ζ.Z5, θ.delta.Δ, v)
    increaseorder!(ζ.Z1, ζ.Z5.c, ζ.Z3.c, convert(T, 2))

    # stage 4
    increaseorder!(ζ.Z6, θ.η.c, ζ.Z1.c)
    increaseorder!(ζ.Z7, θ.B.c, ζ.Z9.c)
    increaseorder!(ζ.Z2, ζ.Z6.c, ζ.Z7.c)

    # RHS of du/dt.
    increaseorder!(ζ.ζ, θ.R.c, ζ.Z2.c)

    return nothing
end


#################### geodesic IVP problem.

struct RQ22ParallelTransport{T}

    # RHS of dv/dt.
    ζ::RQ22ζ{T}

    # variable.
    v::IntegralSequence{T} 
    
    # initial condition.
    v0::Vector{T}
end

# does not create a copy of v0.
function RQ22ParallelTransport(
    a::T,
    b::T,
    v0::Vector{T},
    )::RQ22ParallelTransport{T} where T

    N = length(v0)

    return RQ22ParallelTransport(

        # RHS of du/dt.
        RQ22ζ(a, b, N),

        # variables
        IntegralSequence(T, N),
        
        # initial conditions
        v0,
    )
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

# does not create a copy of x0, u0, v0_set.
function getivpbuffer(
    metric_params::RQ22Metric{T},
    x0::Vector{T},
    u0::Vector{T},
    v0_set::Vector{Vector{T}},
    )::RQ22IVPBuffer{T} where T

    a = metric_params.a
    b = metric_params.b

    @assert length(x0) == length(u0)
    @assert a > zero(T)
    @assert b > zero(T)
    N = length(x0)

    return RQ22IVPBuffer(
        #pt,
        collect( RQ22ParallelTransport(a, b, v0_set[m]) for m = 1:length(v0_set) ),

        # RHS of du/dt.
        RQ22θ(a, b, N),

        # variables
        IntegralSequence(T,N),
        IntegralSequence(T,N),
        
        # initial conditions
        x0, u0,
        #copy(x0), copy(u0), # the supplied initial conditions might be overwritten later. make a copy is safer.
    )
end


# used only for continuity test on x, which doesn't need the computation involving the parallel transported vector field solutions.
# therefore, we ignore resetting the parallel transport field to reduce computation.
function resetbuffer!(
    p::RQ22IVPBuffer{T},
    ::DisableParallelTransport,
    x0::Vector{T},
    u0::Vector{T},
    ) where T

    resetbuffer!(p.θ)
    resetbuffer!(p.u)
    resetbuffer!(p.x)

    p.x0[:] = x0
    p.u0[:] = u0

    return nothing
end
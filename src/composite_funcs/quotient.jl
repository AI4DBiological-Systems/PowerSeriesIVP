##
# use the transforms from elementary.jl to make more complex transforms.

struct Reciprocal{T<:AbstractFloat}
    # [variable index][order index].
    c::Vector{Vector{T}} # power series coefficients.
    b_tilde::Vector{Vector{T}} # buffer, b_tilde in my notes.
    α::Vector{Vector{T}} # buffer, index d is the α(k,0) from (Rodriguez-Bermúdez, 2022) for variable d.
end

function Reciprocal(::Type{T}, N::Integer)::Reciprocal{T} where T

    return Reciprocal(initializecoefficients(T, N), initializecoefficients(T,N), initializecoefficients(T,N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::Reciprocal{T},
    x::Vector{Vector{T}},
    ) where T

    c = A.c
    b_tilde = A.b_tilde
    α = A.α
    
    @assert !isempty(x)
    
    for i in eachindex(c)

        resize!(c[i], 1)
        c[i][begin] = one(T)/x[i][begin]

        resize!(b_tilde[i], 0)
        resize!(α[i], 1)
        α[i][begin] = one(T)
    end

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::Reciprocal{T},
    x::Vector{Vector{T}},
    ) where T

    c = A.c
    b_tilde = A.b_tilde
    α = A.α
    
    @assert !isempty(x)
    #L = length(x[begin])

    for i in eachindex(c)
        @assert (length(c[i]) + 1 == length(x[i]))

        appendreciprocalcoeff!(c[i], x[i], b_tilde[i], α[i])
    end

    return nothing
end

struct ScaledReciprocal{T}
    buf_reciprocal::Reciprocal{T}
    buf_multiply_constant::MultiplyConstant{T}
    c::Vector{Vector{T}}
end

function ScaledReciprocal(::Type{T}, N::Integer)::ScaledReciprocal{T} where T
    r = Reciprocal(T,N)
    m = MultiplyConstant(T,N)
    return ScaledReciprocal(r, m, m.c)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::ScaledReciprocal{T},
    x::Vector{Vector{T}},
    s::T,
    ) where T

    initializeorder!(A.buf_reciprocal, x)
    initializeorder!(A.buf_multiply_constant, A.buf_reciprocal.c, s)

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::ScaledReciprocal{T},
    x::Vector{Vector{T}},
    s::T,
    ) where T

    increaseorder!(A.buf_reciprocal, x)
    increaseorder!(A.buf_multiply_constant, A.buf_reciprocal.c, s)

    return nothing
end

########

struct Quotient{T}
    buf_reciprocal::Reciprocal{T}
    buf_product::Product{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function Quotient(::Type{T}, N::Integer) where T
    
    buf_product = Product(T, N)
    return Quotient(
        Reciprocal(T,N),
        buf_product,
        buf_product.c,
    )
end

# get the sequence for x/y, given the sequence for x and the sequence for y.
function initializeorder!(
    A::Quotient{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    ) where T

    initializeorder!(A.buf_reciprocal, y)
    initializeorder!(A.buf_product, x, A.buf_reciprocal.c)

    return nothing
end

function increaseorder!(
    A::Quotient{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    ) where T

    increaseorder!(A.buf_reciprocal, y)
    increaseorder!(A.buf_product, x, A.buf_reciprocal.c)

    return nothing
end

############

struct ScaledQuotient{T}
    buf_quotient::Quotient{T}
    buf_multiply_constant::MultiplyConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function ScaledQuotient(::Type{T}, N::Integer) where T
    m = MultiplyConstant(T,N)
    return ScaledQuotient(Quotient(T,N), m, m.c)
end

# get the sequence for x/y, given the sequence for x and the sequence for y.
function initializeorder!(
    A::ScaledQuotient{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s::T,
    ) where T

    initializeorder!(A.buf_quotient, x, y)
    initializeorder!(A.buf_multiply_constant, A.buf_quotient.c, s)

    return nothing
end

function increaseorder!(
    A::ScaledQuotient{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s::T,
    ) where T

    increaseorder!(A.buf_quotient, x, y)
    increaseorder!(A.buf_multiply_constant, A.buf_quotient.c, s)

    return nothing
end


############
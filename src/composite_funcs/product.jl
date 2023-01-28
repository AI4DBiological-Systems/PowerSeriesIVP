
struct MultiplyConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function MultiplyConstant(::Type{T}, N::Integer) where T
    return MultiplyConstant(initializecoefficients(T,N))
end

function initializeorder!(
    A::MultiplyConstant{T},
    x::Vector{Vector{T}},
    a,
    ) where T

    return increaseorder!(A, x, a)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::MultiplyConstant{T},
    x::Vector{Vector{T}},
    a::T,
    ) where T

    c = A.c

    for i in eachindex(c)
        #@assert length(x[i]) == length(c[i]) + 1

        push!(c[i], x[i][end]*a)
    end

    return nothing
end

function increaseorder!(
    A::MultiplyConstant{T},
    x::Vector{Vector{T}},
    a::Vector{T},
    ) where T

    c = A.c

    for i in eachindex(c)
        #@assert length(x[i]) == length(c[i]) + 1

        push!(c[i], x[i][end]*a[i])
    end

    return nothing
end


#################

struct Product{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function Product(::Type{T}, N::Integer)::Product{T} where T
    return Product(initializecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::Product{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    ) where T

    return increaseorder!(A, x, y)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::Product{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    ) where T

    #@assert !isempty(x)
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], conv(x[i], y[i]))
    end

    return nothing
end



################

struct Squared{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function Squared(::Type{T}, N::Integer)::Squared{T} where T
    return Squared(initializecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::Squared{T},
    x::Vector{Vector{T}},
    ) where T

    return increaseorder!(A, x)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::Squared{T},
    x::Vector{Vector{T}},
    ) where T

    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], conv(x[i]))
    end

    return nothing
end

#################

struct ScaledProduct{T}
    buf_product::Product{T}
    buf_multiply_constant::MultiplyConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function ScaledProduct(::Type{T}, N::Integer)::ScaledProduct{T} where T
    m = MultiplyConstant(T,N)
    return ScaledProduct(Product(T,N), m, m.c)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::ScaledProduct{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s,
    ) where T

    initializeorder!(A.buf_product, x, y)
    initializeorder!(A.buf_multiply_constant, A.buf_product.c, s)

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::ScaledProduct{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s,
    ) where T

    increaseorder!(A.buf_product, x, y)
    increaseorder!(A.buf_multiply_constant, A.buf_product.c, s)

    return nothing
end
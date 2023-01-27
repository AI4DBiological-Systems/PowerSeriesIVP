struct AdditionConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function AdditionConstant(::Type{T}, N::Integer) where T
    return AdditionConstant(initializecoefficients(T,N))
end

function initializeorder!(
    A::AdditionConstant{T},
    x::Vector{Vector{T}},
    a::Vector{T},
    ) where T

    c = A.c

    for i in eachindex(x)
        resize!(c, 1)
        c[i][begin] = a[i] + x[i][begin]
    end

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::AdditionConstant{T},
    args...
    ) where T

    return nothing
end

struct MultiplyConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function MultiplyConstant(::Type{T}, N::Integer) where T
    return MultiplyConstant(initializecoefficients(T,N))
end

function initializeorder!(
    A::MultiplyConstant{T},
    x::Vector{Vector{T}},
    a::T,
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
        @assert length(x[i]) == length(c[i]) + 1

        push!(c[i], x[i][end]*a)
    end

    return nothing
end

#############
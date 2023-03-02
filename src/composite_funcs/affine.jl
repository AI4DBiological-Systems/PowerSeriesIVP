struct AddConstant{T} <: SingleVariableTrait
    c::Vector{Vector{T}} # [variable index][order index].
end

function AddConstant(::Type{T}, N::Integer) where T
    return AddConstant(allocatecoefficients(T,N))
end

function initializeorder!(
    A::AddConstant{T},
    x::Vector{Vector{T}},
    a::Vector{T},
    ) where T

    c = A.c

    for i in eachindex(x)
        push!(c[i], a[i] + x[i][begin])
    end

    return nothing
end

function initializeorder!(
    A::AddConstant{T},
    x::Vector{Vector{T}},
    a::T,
    ) where T

    c = A.c

    for i in eachindex(x)
        push!(c[i], a + x[i][begin])
    end

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::AddConstant{T},
    x::Vector{Vector{T}},
    args...
    ) where T

    c = A.c

    for i in eachindex(x)
        push!(c[i], x[i][end])
    end

    return nothing
end

################## 

struct Subtraction{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function Subtraction(::Type{T}, N::Integer)::Subtraction{T} where T
    return Subtraction(allocatecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(A::Subtraction{T}, x::Vector{Vector{T}}, y::Vector{Vector{T}}) where T

    return increaseorder!(A,x,y)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(A::Subtraction{T}, x::Vector{Vector{T}}, y::Vector{Vector{T}}) where T
    #@assert !isempty(x)
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], x[i][end] - y[i][end])
    end

    return nothing
end

#############

struct Addition{T} <: SingleVariableTrait
    c::Vector{Vector{T}} # [variable index][order index].
end

function Addition(::Type{T}, N::Integer)::Addition{T} where T
    return Addition(allocatecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(A::Addition{T}, x::Vector{Vector{T}}, y::Vector{Vector{T}}) where T

    return increaseorder!(A,x,y)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(A::Addition{T}, x::Vector{Vector{T}}, y::Vector{Vector{T}}) where T
    #@assert !isempty(x)
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], x[i][end] + y[i][end])
    end

    return nothing
end

#############

struct LinearOperation{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function LinearOperation(::Type{T}, N::Integer)::LinearOperation{T} where T
    return LinearOperation(allocatecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::LinearOperation{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x,
    s_y,
    ) where T

    return increaseorder!(A, x, y, s_x, s_y)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::LinearOperation{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x::T,
    s_y::T,
    ) where T

    @assert !isempty(x)
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], s_x*x[i][end] + s_y*y[i][end])
    end

    return nothing
end

function increaseorder!(
    A::LinearOperation{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x::Vector{T},
    s_y::Vector{T},
    ) where T

    @assert !isempty(x)
    @assert length(x) == length(y) == length(s_x) == length(x_y)

    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], s_x[i]*x[i][end] + s_y[i]*y[i][end])
    end

    return nothing
end

###

struct SubtractFromScaled{T} <: SingleVariableTrait
    c::Vector{Vector{T}} # [variable index][order index].
end

function SubtractFromScaled(::Type{T}, N::Integer)::SubtractFromScaled{T} where T
    return SubtractFromScaled(allocatecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::SubtractFromScaled{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x,
    ) where T

    return increaseorder!(A, x, y, s_x)
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::SubtractFromScaled{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x::T,
    ) where T

    @assert !isempty(x)
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], s_x*x[i][end] - y[i][end])
    end

    return nothing
end

function increaseorder!(
    A::SubtractFromScaled{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    s_x::Vector{T},
    ) where T

    @assert !isempty(x)
    @assert length(x) == length(y) == length(s_x) == length(x_y)
    
    c = A.c

    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], s_x[i]*x[i][end] - y[i][end])
    end

    return nothing
end
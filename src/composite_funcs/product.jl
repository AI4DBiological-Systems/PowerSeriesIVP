

##############

struct Product{T}
    # [variable index][order index].
    c::Vector{Vector{T}} # power series coefficients.
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

    c = A.c
    
    @assert !isempty(x)
    
    for i in eachindex(c)

        resize!(c[i], 1)
        c[i][begin] = x[i][begin]*y[i][begin]
    end

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::Product{T},
    x::Vector{Vector{T}},
    y::Vector{Vector{T}},
    ) where T

    c = A.c
    
    @assert !isempty(x)
    #L = length(x[begin])

    for i in eachindex(c)
        @assert (length(c[i]) + 1 == length(x[i]))

        push!(c[i], conv(x[i], y[i]))
    end

    return nothing
end


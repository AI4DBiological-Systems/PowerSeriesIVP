## 
# given at least the 0-th order term, compute and append the L-th coefficient of each composition transformation.


# this is an internal data structure, not a forward-facing Taylor polynomial type.
# That is why it doesn't have the field c.
# The fields Δ and Δ_sq are two different coefficient fields.
# Internal data structures for used for defining forward-facing Taylor polynomial types.
struct InterVariableDifference{T}
    # not using Array{T,3} since we might resize the order on-the-fly.
    Δ::Matrix{Vector{T}} # i x j, order.
    Δ_sq::Matrix{Vector{T}} # i x j, order.
end

function InterVariableDifference(::Type{T}, N::Integer) where T
    return InterVariableDifference(
        collect( zeros(T, 0) for _ = 1:N, _ = 1:N),
        collect( zeros(T, 0) for _ = 1:N, _ = 1:N),
    )
end

function initializeorder!(
    A::InterVariableDifference{T},
    c::Vector{Vector{T}},
    ) where T

    return increaseorder!(A, c)
end

function increaseorder!(
    A::InterVariableDifference{T},
    c::Vector{Vector{T}},
    ) where T

    Δ = A.Δ
    Δ_sq = A.Δ_sq

    # traverse lower triangle.
    for j in eachindex(c)
        for i = j+1:length(c)
            #@show (i,j)

            @assert length(Δ[i,j]) == length(Δ[j,i]) == length(c[i]) - 1 == length(c[j]) - 1

            tmp = c[i][end] - c[j][end]
            
            push!(Δ[i,j], tmp)
            push!(Δ[j,i], -tmp)

            # product: conv.
            @assert length(Δ_sq[i,j]) == length(Δ[i,j]) - 1

            tmp_sq = conv(Δ[i,j])
            push!(Δ_sq[i,j], tmp_sq)
            push!(Δ_sq[j,i], tmp_sq)
        end
    end

    return nothing
end


###################

struct ΔSumCol{T} # W4. Skips Δ[i,i] in the sum.
    c::Vector{Vector{T}} # [variable index][order index].
end

function ΔSumCol(::Type{T}, N::Integer)::ΔSumCol{T} where T
    return ΔSumCol(initializecoefficients(T,N))
end

function initializeorder!(A::ΔSumCol{T}, S::Matrix{Vector{T}}) where T
    return increaseorder!(A, S)
end

function increaseorder!(A::ΔSumCol{T}, S::Matrix{Vector{T}}) where T

    c = A.c
    for i in axes(S,2)
        #@assert length(c[i]) == length(S[k,i]) - 1

        # sum from beginning up to index i.
        # tmp = sum( S[k,i][end] for k in Iterators.take(axes(S,1), i-1) )
        # tmp += sum( S[k,i][end] for k in Iterators.drop(axes(S,1), i) )

        # explicit loop since the iterator might be empty, which sum() cannot handle.
        tmp = zero(T)
        for k in Iterators.take(axes(S,1), i-1)
            tmp += S[k,i][end]
        end
        
        # sum from index i+1 to end
        for k in Iterators.drop(axes(S,1), i)
            tmp += S[k,i][end]
        end

        # slow version. Slow due to the if statement that is checked in every iteration.
        #tmp = sum( S[k,i][end] for k in axes(S,1) if i != k)

        push!(c[i], tmp)
    end

    return nothing
end


###### metric problem specific.

struct ΔSumColProduct{T}
    buf_product_mat::Matrix{Vector{T}}
    c::Vector{Vector{T}} # [variable index][order index].
end

function ΔSumColProduct(::Type{T}, N::Integer)::ΔSumColProduct{T} where T
    return ΔSumColProduct(
    collect( zeros(T, 0) for _ = 1:N, _ = 1:N),
    initializecoefficients(T,N),
    )
end

function initializeorder!(
    Z::ΔSumColProduct{T},
    Δ::Matrix{Vector{T}},
    x::Vector{Vector{T}},
    ) where T
    
    return increaseorder!(Z, Δ, x)
end

function increaseorder!(
    Z::ΔSumColProduct{T},
    Δ::Matrix{Vector{T}},
    x::Vector{Vector{T}},
    ) where T

    c = Z.c
    B = Z.buf_product_mat

    # product.
    for i in axes(Δ,2)
        for k in axes(Δ,1)
            #@assert (length(c[i]) + 1 == length(x[i]))
            push!(B[k,i], conv(Δ[k,i], x[k]))
        end
    end

    # sum.
    for i in eachindex(c)
        #@assert length(c[i]) == length(B[d1,i]) - 1

        tmp = zero(T)
        for k in Iterators.take(axes(B,1), i-1)
            tmp += B[k,i][end]
        end
        
        # sum from index i+1 to end
        for k in Iterators.drop(axes(B,1), i)
            tmp += B[k,i][end]
        end

        # slow version.
        # tmp = sum( B[k,i][end] for k in axes(B,1) if k != i)
        
        push!(c[i], tmp)
    end

    return nothing
end

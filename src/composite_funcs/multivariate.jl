## 
# given at least the 0-th order term, compute and append the L-th coefficient of each composition transformation.



struct SelfSquaredLoss{T}
    # not using Array{T,3} since we might resize the order on-the-fly.
    Δ::Matrix{Vector{T}} # i x j, order.
    Δ_sq::Matrix{Vector{T}} # i x j, order.
end

function SelfSquaredLoss(::Type{T}, N::Integer, L::Integer) where T
    return SelfSquaredLoss(
        collect( zeros(T, L) for _ = 1:N, _ = 1:N),
        collect( zeros(T, L) for _ = 1:N, _ = 1:N),
    )
end

function initializeorder!(
    A::SelfSquaredLoss{T},
    c::Vector{Vector{T}},
    ) where T

    return increaseorder!(A, c)
end

function increaseorder!(
    A::SelfSquaredLoss{T},
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

            tmp_sq = conv(A.Δ[i,j])
            push!(Δ_sq[i,j], tmp_sq)
        end
    end

    return nothing
end

struct SumCols{T} # W4.
    c::Vector{Vector{T}} # [variable index][order index].
end

function SumCols(::Type{T}, N::Integer)::SumCols{T} where T
    return SumCols(initializecoefficients(T,N))
end

function initializeorder!(A::SumCols{T}, S::Matrix{Vector{T}}) where T

    return increaseorder!(A, S)
end

function increaseorder!(A::SumCols{T}, S::Matrix{Vector{T}}) where T

    c = A.c

    for d1 in axes(S,1)
        @assert length(c[d1]) == length(S[d1,d2]) - 1

        tmp = sum( S[d1,d2][end] for d2 in axes(S,2) )
        push!(c[d1], tmp)
    end

    return nothing
end





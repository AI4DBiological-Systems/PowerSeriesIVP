
# struct Coefficients{T}
#     c::Vector{Vector{T}} # [variable index][order index].
# end

struct SelfSquareLoss{T}
    
    # not using Array{T,3} since we might resize the order on-the-fly.
    Δ::Matrix{Vector{T}} # i x j, order.
    Δ_sq::Matrix{Vector{T}} # i x j, order.
end

function SelfSquareLoss(::Type{T}, N::Integer, L::Integer) where T
    return SelfSquareLoss(
        collect( zeros(T, L) for _ = 1:N, _ = 1:N),
        collect( zeros(T, L) for _ = 1:N, _ = 1:N),
    )
end

function increaseorder!(
    A::SelfSquareLoss{T},
    c::Vector{Vector{T}},
    ) where T

    # traverse lower triangle.
    for j in eachindex(c)
        for i = j+1:length(c)
            #@show (i,j)

            tmp = c[i][end] - c[j][end]
            
            push!(A.Δ[i,j], tmp)
            push!(A.Δ[j,i], -tmp)

            # product: conv.
            tmp_sq = conv(A.Δ[i,j])
            push!(A.Δ_sq[i,j], tmp_sq)
        end
    end

    return nothing
end

struct SummedCols{T} # W4.
    c::Vector{Vector{T}} # [variable index][order index].
end

function increaseorder!(
    A::SummedCols{T},
    S::Matrix{Vector{T}},
    ) where T

    c = A.c

    for d1 in axes(S,1)
        tmp = sum( S[d1,d2][end] for d2 in axes(S,2) )
        push!(c[d1], tmp)
    end

    return nothing
end

struct AddedConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::AddedConstant{T},
    args...
    ) where T

    return nothing
end

struct MultipliedConstant{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::MultipliedConstant{T},
    x::Vector{Vector{T}},
    a::T,
    ) where T

    c = A.c

    for i in eachindex(c)
        push!(c[i], x[i][end]*a)
    end

    return nothing
end


struct Reciprocal{T}
     # [variable index][order index].
    c::Vector{Vector{T}} # power series coefficients.
    b_tilde::Vector{Vector{T}} # buffer, b_tilde in my notes.
    α::Vector{Vector{T}} # buffer, index d is the α(k,0) from (Rodriguez-Bermúdez, 2022) for variable d.
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::Reciprocal{T},
    x::Vector{Vector{T}},
    ) where T

    c = A.c

    for i in eachindex(c)
        push!(c[i], x[i][end]*a)
    end

    return nothing
end
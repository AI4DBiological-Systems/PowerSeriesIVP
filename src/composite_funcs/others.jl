# uses the fundamental theorem of calculus to get a first-derivative's power series without numerical integration.
struct IntegralSequence{T}
    c::Vector{Vector{T}} # [variable index][order index].
end

function IntegralSequence(::Type{T}, N::Integer)::IntegralSequence{T} where T
    return IntegralSequence(initializecoefficients(T, N))
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function initializeorder!(
    A::IntegralSequence{T},
    x0::Vector{T},
    ) where T

    c = A.c
    for i in eachindex(c)
        #@assert (length(c[i]) + 1 == length(x[i]))
        push!(c[i], x0[i])
    end

    return nothing
end

# adding a constant only affects order zero, which wouldn't show up with increaseorder!(), but rather initializeorder0.
function increaseorder!(
    A::IntegralSequence{T},
    x::Vector{Vector{T}}, # [variable index][order index].
    ) where T

    # power series order starts at 0, but Vector{T} is 1-indexing.
    # - therefore, we need to -1 to get the power series order.
    # power_series_order_of_x = length(x[begin]) - 1
    # next_order = power_series_order_of_x + 1
    next_order = length(x[begin])
    
    c = A.c
    for i in eachindex(c) #Iterators.drop(eachindex(c), 1)
        
        push!(c[i], x[i][end]/next_order)
    end

    return nothing
end

# # This is for updating from an x that has just been updated, but we really want to simultaneously update A.c along with x.c.
# # we'd need to look at the second last entry instead of the last entry.
# function increaseorderfromderivative!(
#     A::IntegralSequence{T},
#     x::Vector{Vector{T}}, # [variable index][order index].
#     ) where T

#     next_order = length(x[begin])

#     c = A.c
#     for i in eachindex(c) #Iterators.drop(eachindex(c), 1)
#         #@assert (length(c[i]) + 1 == length(x[i]))
#         #push!(c[i], x[i][end]/(length(x[i])+1))
#         push!(c[i], x[i][end-1]/next_order)
#     end

#     return nothing
# end
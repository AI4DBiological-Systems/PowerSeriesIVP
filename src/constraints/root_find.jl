# root finding routine, based on Decartes rule of sign change.

# unused.
function findrootinterval(
    x::Vector{Vector{T}},
    #as::Vector{Vector{T}},
    #bs::Vector{T},
    a::Vector{T},
    b::T,
    ) where T

    @assert length(x) == length(as) == length(bs)

    # out.position[d] = evaltaylor(c.x[d], t, a)
    L_p1 = length(x[d][begin])
    c = Vector{T}(undef, L_p1)

    for l in eachindex(x[d][begin])
        for d in eachindex(x)
            c[l] = x[d][l]*a[d]
        end
    end
    c[begin] -= b

    # need Taylor shift.
    # https://math.stackexchange.com/questions/694565/polynomial-shift

    findfirstroot()
end


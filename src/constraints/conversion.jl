######## prepare the intersection polynomials from the given Taylor polynomial x and the affine cosntraints as, bs.

# mutates cs.
# the intersection polynomial operate on the variable t-t0, t0 is the expansion point, a constant.
function updateintersectionpolynomials!(
    cs::Vector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    as::Vector{Vector{T}}, # [constraints][variable]
    bs::Vector{T}, # [constraints]
    ) where T

    N_constraints = length(as)
    @assert N_constraints == length(bs)
    resize!(cs, N_constraints)

    L_p1 = length(x[begin])

    for m in eachindex(cs)
        
        #cs[m] = zeros(T, L_p1)
        resize!(cs[m], L_p1)
        fill!(cs[m], zero(T))

        a = as[m]
        b = bs[m]

        for d in eachindex(x)
            for l in eachindex(x[d])
                cs[m][l] += x[d][l]*a[d]
            end
        end
        cs[m][begin] -= b
    end

    return nothing
end

# c[i] is the coefficient of the i-th degree monomial. i ∈ {0,1,2,...,L}.
# modifies c such that the last entry, i.e. the coefficient for the L-th degree, is 1.
function standardizecoefficients!(c::Vector{T}) where T

    Z = c[end]
    for i in eachindex(c)
        c[i] = c[i]/Z
    end

    return nothing
end
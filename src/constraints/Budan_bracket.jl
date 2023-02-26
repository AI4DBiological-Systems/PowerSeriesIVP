

### Budan's upper bound on the number of real roots of a polynomial equation on an interval.




# no error checking.
function copybuffer!(x::Vector{Vector{T}}, src::Vector{Vector{T}}) where T
    
    resize!(x, length(src))

    for m in eachindex(x)
        resize!(x[m], length(src[m]))
        for l in eachindex(x)
            x[m][l] = src[m][l]
        end
    end

    return nothing
end


################### Budan's theorem for smallest positive real root isolation.

function setupbinomialcoefficients(order::Integer)::Matrix{Int}
    @assert order > 0
    
    A = Matrix{Int}(undef, order, order)
    for m in axes(A, 2)
        for n in axes(A, 1)
            A[n,m] = binomial(n,m)
        end
    end
    
    return A
end

function updateshiftedpolynomial!(
    out::Vector{T},
    c::Vector{T},
    b::T,
    a::T,
    binomial_mat::Matrix{Int}, # output of setupbinomialcoefficients()
    ) where T

    resize!(out, length(c))
    L = length(c) - 1

    out[begin] = evaltaylor(c, b, a)

    x = b-a
    for m = 1:L
        # TODO: use Horner's method to speed up the evaluation here.
        out[m+1] = sum( c[n+1] * binomial_mat[n,m] * x^(n-m) for n = m:L )
    end
    return nothing
end

function updateshiftedpolynomial!(
    cs_shifted::Vector{Vector{T}},
    cs::Vector{Vector{T}},
    b::T,
    a::T,
    binomial_mat::Matrix{Int}, # output of setupbinomialcoefficients()
    ) where T

    @assert length(cs) == length(cs_shifted)

    for m in eachindex(cs)
        updateshiftedpolynomial!(cs_shifted[m], cs[m], b, a, binomial_mat)
    end

    return nothing
end

# mutates A.
function upperboundintersections!(
    A::BudanIntersectionBuffers{T},
    sol::PiecewiseTaylorPolynomial{T},
    constraints,
    bino_mat::Matrix{Int},
    ) where T
    
    N_pieces = length(sol.coefficients)
    @assert length(sol.steps) == N_pieces

    ubs = Vector{T}(undef, N_pieces)
    inds = Vector{Int}(undef, N_pieces)

    for n in eachindex(sol.coefficients)

        ubs[n], inds[n] = upperboundNroots!(
            A,
            sol.coefficients[n].x,
            sol.steps[n],
            constraints,
            bino_mat,
        )
    end

    return ubs, inds
end


# based on findfirstroot!()
function upperboundNroots!(
    A::BudanIntersectionBuffers{T},
    x::Vector{Vector{T}},
    h::T,
    constraints::AffineConstraints{T},
    bino_mat::Matrix{Int},
    ) where T

    updateintersectionpolynomials!(
        A.cs, x, 
        constraints.normals, constraints.offsets,
    )

    cs = A.cs
    cs_right = A.cs_right

    # initialize.
    #copybuffer!(A.cs_left, A.cs)
    #cs_left = A.cs_left
    cs_left = cs
    
    x_right = h

    # get shifted polynomial.
    updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

    max_val, max_ind = findmax( upperboundroots(cs_left[m], cs_right[m]) for m in eachindex(cs_left))
    
    return max_val, max_ind
end



function upperboundroots(c_left::Vector{T}, c_right::Vector{T})::Int where T
    return countsignflips(c_left) - countsignflips(c_right)
end

function countsignflips(c::Vector{T})::Int where T
    prev_sign = signbit(c[begin]) # true for negative sign.

    N_flips = 0
    for i in Iterators.drop(eachindex(c),1)
        
        current_sign = signbit(c[i])
        
        if abs(c[i]) > eps() # ignore coefficients with very small magnitude.
            if prev_sign != current_sign
                N_flips += 1
                
                prev_sign = current_sign
            end
        end
    end

    return N_flips
end


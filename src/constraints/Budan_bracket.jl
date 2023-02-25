

### Budan's upper bound on the number of real roots of a polynomial equation on an interval.


function allocatebuffer!(
    A::BudanIntersectionBuffers{T},
    cs::Vector{Vector{T}},
    ) where T

    N_constraints = length(cs)
    @assert N_constraints > 0

    resize!(A.cs_left, N_constraints)
    resize!(A.cs_right, N_constraints)
    resize!(A.cs_center, N_constraints)

    L = length(cs[begin]) - 1
    
    for m in eachindex(cs)
    
        L_p1 = length(cs[m])
        @assert L == L_p1-1
    
        resize!(A.cs_left[m], L_p1)
        resize!(A.cs_right[m], L_p1)
        resize!(A.cs_center[m], L_p1)
    end

    return nothing
end

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

    updateintersectionpolynomial!(
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


# in development. not used for now.
function findfirstroot!(
    A::BudanIntersectionBuffers{T},
    h::T,
    bino_mat::Matrix{Int};
    min_bracket_len = 1e-4,
    max_iters = 100000,
    ) where T

    # initialize.
    #allocatebuffer!(A, cs)
    copybuffer!(A.cs_left, A.cs)

    cs = A.cs
    cs_left = A.cs_left
    cs_right = A.cs_right
    cs_center = A.cs_center

    x_left = zero(T)
    x_right = h

    # get shifted polynomial.
    updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

    max_val, max_ind = findmax( upperboundroots(cs_left[m], cs_right[m]) for m in eachindex(cs_left))
    if max_val < 0
        return h, :redo_quartic # error! contracticts Budan's theorem. Redo this solution piece using quartic polynomials.
    end

    if max_val == 0
        # Guaranteed no roots in this interval.
        return h, :continue_simulation # safe to declare no intersection on current piece's interval.
        
    elseif max_val == 1
        # there is a root in this interval.
        x_intersection = runITP(x_left, x_right)

        if isinite(x_intersection)
            return x_intersection, :end_simulation
        end

    end
    
    # more than one root in this interval, or we failed to find the root in the interval.
    # we shall subdivide the interval further, and attempt to find the root, or get a new
    # starting time that is less than h for the next solution polynomial for the IVP simulation.
    
    new_h, instruction = decideintersectionbracket!(
        cs_center,
        cs_left,
        cs_right,
        x_left,
        x_right,
        h,
        cs,
        bino_mat;
        min_bracket_len = min_bracket_len,
        max_iters = max_iters,
    )
    
    return new_h, instruction
end

# next simulation is either use quartic degree, or we end simulation on the current polynomial.
function decideintersectionbracket!(
    cs_center::Vector{Vector{T}},
    cs_left::Vector{Vector{T}},
    cs_right::Vector{Vector{T}},
    x_left::T,
    x_right::T,
    x_ub::T,
    cs::Vector{Vector{T}},
    bino_mat::Matrix{Int};
    min_bracket_len = min_bracket_len,
    max_iters = max_iters,
    ) where T

    # when we return with an indication to use quartic polynomials in the next solution piece:
    # set the starting time as x_left for the next polynomial piece, because
    # we hope to use explicit quartic root solve formula to get the intersection in the next solution piece.

    for _ = 1:max_iters

        x_center = (x_left+x_right)/2
        ub_roots, m_ind, region = splitquery!(
            cs_center,
            cs_left,
            cs_right,
            x_left,
            x_right,
            x_center,
            cs,
            bino_mat;
            min_bracket_len = min_bracket_len,
        )

        # see if we can attempt to solve for the root.
        if ub_roots == 1
            
            x_intersection = runITP()
            if isinite(x_intersection)
                return x_intersection, :end_simulation
            end

            # in case we couldn't find a root using ITP, for whatever reason.
            if region == :right

                return x_center, :use_quartic_next
            else
                return x_left, :use_quartic_next
            end
        

        elseif ub_roots == 0
            # this only happens if instructions == :right: see splitquery!()
            
            # no roots in the [x_left, x_right] interval.
            if (x_ub - x_right) < min_bracket_len
                return x_right, :continue_simulation
            end

            # next interval is: [x-right, x_ub].
            x_left = x_right
            x_right = x_ub
        
        # we shouldn't solve for the root, so we subdivide the interval further.
        elseif region == :tol_reached
            # the region between [x_left, x_right] does not have Budan's root upperbound of 1
            # but, the length of this region is smaller than min_bracket_len.
            return x_right, :continue_simulation

        elseif region == :right
            x_left = x_center
            #cs_left, cs_center = cs_center, cs_left
        else
            x_right = x_center
            #cs_right, cs_center = cs_center, cs_right
        end

    end
    
    return x_left, :use_quartic_next
end


# the assumption is that Budan's test for the interval [x_left, x_right] is non-zero.
function splitquery!(
    cs_center,
    cs_left,
    cs_right,
    x_left, x_right, x_center, cs, bino_mat;
    min_bracket_len = 1e-4)
    
    if abs(x_left-x_right) < min_bracket_len
        return tol_reached
    end

    #x_center = (x_left + x_right)/2

    # # query left.
    
    updateshiftedpolynomial!(cs_center, cs, x_center, zero(T), bino_mat)
    updateshiftedpolynomial!(cs_left, cs, x_left, zero(T), bino_mat)

    max_val1, max_ind1 = findmax( upperboundroots(cs_left[m], cs_center[m]) for m in eachindex(cs_left))
    
    if max_val1 == 0
        # # query right.

        updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

        max_val12, max_ind2 = findmax( upperboundroots(cs_center[m], cs_right[m]) for m in eachindex(cs_center))

        return max_val12, max_ind2, :right
    end

    return max_val1, max_ind1, :left
end

function upperboundroots(c_left::Vector{T}, c_right::Vector{T})::Int where T
    return countsignflips(c_left) - countsignflips(c_right)
end

function countsignflips(c::Vector{T})::Int where T
    prev_sign = signbit(c[begin]) # true for negative sign.

    N_flips = 0
    for i in Iterators.drop(eachindex(c),1)
        
        current_sign = signbit(c[i])
        
        if prev_sign != current_sign
            N_flips += 1
            
            prev_sign = current_sign
        end
    end

    return N_flips
end


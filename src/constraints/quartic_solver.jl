
######### routines that solve univariate cubic and quartic polynomial equations.

# http://numerical.recipes/book/book.html.
"""
solvecubicequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}} where T

solves the cubic equation z^3 + sum( a[begin+i]*z^i for i = 1:2 ) + a[begin]
i.e., a[begin+i]  is the i-th coefficient for monomial z^i, a[begin] is the constant coefficient.

Denote by ri the i-th root.
This function returns three complex numbers.
"""
function solvecubicequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}} where T
    @assert length(a) >= 3
    return solvecubicequation(a[begin], a[begin+1], a[begin+2])
end

function solvecubicequation(a0::T, a1::T, a2::T)::Tuple{Complex{T}, Complex{T}, Complex{T}} where T

    q = a1/3 - a2^2/9
    r = (a1*a2 - 3*a0)/6 - a2^3/27

    q_cubed = q^3
    condition_eval = r^2 + q_cubed
    if condition_eval > zero(T)
        
        # case: one real solution.
        A = cbrt( abs(r) + sqrt(condition_eval) )
        q_div_A = q/A

        t1 = A - q_div_A
        if r < zero(T)
            t1 = -t1
        end

        z1 = t1 - a2/3
        
        x2 = -t1/2 - a2/3 # real part of z2.
        y2 = (A+q_div_A)*sqrt(3)/2 # imaginary part of z2

        return Complex(z1), Complex(x2, y2), Complex(x2, -y2)
    end

    # case all are real solutions.
    θ = zero(T)
    if q_cubed < zero(T)
        θ = acos(r/(sqrt(-q_cubed)))
    end

    ϕ1 = θ/3
    ϕ2 = ϕ1 - 2*π/3
    ϕ3 = ϕ1 + 2*π/3
    
    sqrt_mq = zero(T)
    z1 = -a2/3
    z2 = z1
    z3 = z1
    if q < zero(T)
        sqrt_mq = sqrt(-q)
        z1 += 2*sqrt_mq*cos(ϕ1)
        z2 += 2*sqrt_mq*cos(ϕ2)
        z3 += 2*sqrt_mq*cos(ϕ3)
    end

    return Complex(z1), Complex(z2), Complex(z3) # roots are in decending order.
end

"""
solvequarticequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}, Complex{T}} where T

Returns 4 complex numbers.
"""
# David Wolter's modified Euler algorithm: https://quarticequations.com/
function solvequarticequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}, Complex{T}} where T
    
    @assert length(a) >= 4

    A0, A1, A2, A3 = a # ingnores the entries after 4.
    C = A3/4
    
    # # prepare resolvent cubic equation intermediates.
    C_sq = C^2
    C_cube = C_sq*C
    C_quartic = C_cube*C

    b2 = A2 - 6*C_sq
    b1 = A1 - 2*A2*C + 8*C_cube
    b0 = A0 - A1*C + A2*C_sq - 3*C_quartic
    σ = -1
    if b1 > zero(T) # could use flipsign() here, but -0.0 and 0.0 worries me. however, -0.0 == 0.0 is true on my machine.
        σ = 1
    end

    # # solve resolvent cubic equation.
    a0 = -b1^2/64
    a1 = (b2^2 - 4*b0)/16
    a2 = b2/2
    z1, z2, z3 = solvecubicequation(a0, a1, a2)
    x1 = real(z1)
    x2 = real(z2)
    x3 = real(z3)
    y2 = imag(z2)

    # # solve for quartic roots.
    # y1 should be zero.
    # x1 should b the largest real root, and non-negative.
    # Either y2 and y3 are both zeros, or ( y3 = -y2 and  x2 = x3 )
    # x2*x3 + y2^2 should be non-negative.

    # # no need to do explicit check here. Do explicit check on constraint qualification later.
    inner_term = x2*x3 + y2^2
    # if inner_term < zero(T) || x1 < zero(T)
    #     return Complex(NaN), Complex(NaN), Complex(NaN), Complex(NaN)
    # end
    # σ_term = 2*σ*sqrt(inner_term)
    # term1 = sqrt(x1)
    
    term1 = sqrt(Complex(x1))

    σ_term = 2*σ*sqrt(Complex(inner_term))
    x2_plus_x3 = x2+x3
    
    term12 = sqrt( (x2_plus_x3 - σ_term) )
    term34 = sqrt( (x2_plus_x3 + σ_term) )

    r12 = term1 - C
    T1 = r12 + term12
    T2 = r12 - term12
    
    r34 = -term1 - C
    T3 = r34 + term34
    T4 = r34 - term34

    return T1, T2, T3, T4
end

function solvequarticequation!(roots::Vector{Complex{T}}, a::Vector{T}) where T
    return roots[1], roots[2], roots[3], roots[4] = solvequarticequation(a)
end

function isapproxreal(z::Complex{T}; atol = 1e-14) where T
    if abs(imag(z)) < atol
        return true
    end
    
    return false
end

# returns one of the minimum positive (approximately) real entries in zs.
# returns Inf if there are no approximately real entries in zs.
function findminpositivereal(
    zs::Vector{Complex{T}};
    atol = 1e-14,
    )::T where T <: AbstractFloat
    
    r = convert(T, Inf)

    for i in eachindex(zs)
        if isapproxreal(zs[i]; atol = atol)
            zr = real(zs[i])
            if zr > 0
                r = min(r, zr)
            end
        end
    end

    return r
end

# second output is 0 if no roots. An index between [0, N_constraints] if a real root was found.
function refinestep!(
    A::RootsBuffer{T},
    h::T,
    x::Vector{Vector{T}},
    constraints::ConstraintType,
    )::Tuple{T,Int} where T <: AbstractFloat
    
    cs = A.intersection_coefficients
    all_roots = A.all_roots
    smallest_positive_roots = A.smallest_positive_roots
    atol = A.zero_tol

    updateintersectionpolynomials!(cs, x, constraints)

    for m in eachindex(cs)
        standardizecoefficients!(cs[m]) # put in standard form.

        solvequarticequation!(all_roots, cs[m])
        #@show all_roots, cs[m]
        #println()
        smallest_positive_roots[m] = findminpositivereal(all_roots; atol = atol)
    end
    t, constraint_ind = findmin(smallest_positive_roots)
    #@show smallest_positive_roots
    if zero(T) < t <= h
        return t, constraint_ind
    end

    return h, 0
end

# front end.
function refinestep!(
    C::ConstraintsContainer{T},
    h::T,
    x::Vector{Vector{T}},
    )::Tuple{T,Int} where T <: AbstractFloat

    return refinestep!(
        C.explicit_roots_buffer,
        h,
        x,
        C.constraints,
    )
end

# front end.
function refinestep!(
    ::NoConstraints,
    h::T,
    args...
    )::Tuple{T,Int} where T <: AbstractFloat
    
    return h, 0
end

# mutates intersection_buf.
function searchintersection!(
    intersection_buf::RootsBuffer{T},
    sol::PiecewiseTaylorPolynomial{T},
    constraints,
    ) where T
    
    N_pieces = length(sol.coefficients)
    @assert length(sol.steps) == N_pieces

    ts = Vector{T}(undef, N_pieces)
    inds = Vector{Int}(undef, N_pieces)

    for n in eachindex(sol.coefficients)

        ts[n], inds[n] = refinestep!(
            intersection_buf,
            sol.steps[n],
            sol.coefficients[n].x,
            constraints,
        )
    end

    return ts, inds
end
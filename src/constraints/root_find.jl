# root finding routine, based on Decartes rule of sign change.

# I am here. utilize the endpoints of the intervals of the solution's coeff, instead of Taylor exp of the solution.
# - the idea is to see if the IVP solution have roots, if it were atruncated Taylor polynomial. If there might be roots, then we use quartic order over this interval.
# then write routine that loops this check over all constraints.
# then write routine that resolves the current piece if there is a root.
# - solve using order 4, then Budan interval val check. if might have root again, then solve using quartic.
# Review Budan's theorem again.
function Budanintervalcheck(
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

# http://numerical.recipes/book/book.html.
"""
solvecubicequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}} where T

solves the cubic equation z^3 + sum( a[begin+i]*z^i for i = 1:2 ) + a[begin]
i.e., a[begin+i]  is the i-th coefficient for monomial z^i, a[begin] is the constant coefficient.

Denote by ri the i-th root.
This function returns three complex numbers.
"""
function solvecubicequation(a::Vector{T})::Tuple{Complex{T}, Complex{T}, Complex{T}} where T
    @assert length(a) == 3
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
    
    @assert length(a) == 4

    A0, A1, A2, A3 = a
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

function isapproxreal(z::Complex{T}; atol = 1e-14) where T
    if abs(imag(z)) < atol
        return true
    end
    
    return false
end
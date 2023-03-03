# run intersect.jl first.

########## test root solving code. Cubic and quartic euqations.

c_cube = randn(3)
z1, z2, z3 = PowerSeriesIVP.solvecubicequation(c_cube)
cubefunc = tt->(c_cube[begin]+ c_cube[begin+1]*tt + c_cube[begin+2]*tt^2 + tt^3)

@show cubefunc(z1), cubefunc(z2), cubefunc(z3)
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3)
println()

"""
graph of x^4 + 3*x^3 - 2*x^2 +4.3*x -3
x == -3.85917 || x == 0.657164
{x == 0.101001 +/- 1.08292 im}
"""

c_test = [-3; 4.3; -2; 3]
z1, z2, z3, z4 = PowerSeriesIVP.solvequarticequation(c_test)
@show z1, z2, z3, z4
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3), PowerSeriesIVP.isapproxreal(z4)
println()

"""
graph of x^4 - 0.1*x^3 - 2*x^2 +0.3*x +0.1
x ≈ -1.4211 || x ≈ -0.16203 || x ≈ 0.31819 || x ≈ 1.36493
"""

c_test = [0.1; 0.3; -2; -0.1]
z1, z2, z3, z4 = PowerSeriesIVP.solvequarticequation(c_test)
@show z1, z2, z3, z4
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3), PowerSeriesIVP.isapproxreal(z4)
println()



############### intersection code.

###### shifted polynomial algorithm: set up test scenario.
piece_select = 17
c_p = sol.coefficients[piece_select].x
c_u = sol.coefficients[piece_select].u
h = sol.steps[piece_select]
t0 = sol.expansion_points[piece_select]

c_p_next = sol.coefficients[piece_select+1].x
c_u_next = sol.coefficients[piece_select+1].u

# N_constraints = 5
# as = collect( randn(N_vars) for _ = 1:N_constraints )
# bs = randn(N_constraints)

constraint_select = 2
a = as[constraint_select]
b = bs[constraint_select]


########## general root isolation. investigate in the future.

L = Lm1_fixed + 1
T = Float64


all_roots = Vector{Complex{T}}(undef, L)
smallest_positive_roots = Vector{T}(undef, N_constraints) # real roots.



complex_zero_tol = 1e-8
intersection_buf = PowerSeriesIVP.RootsBuffer(
    complex_zero_tol,
    L,
    N_constraints,
)


t_intersect, constraint_ind = PowerSeriesIVP.refinestep!(
    intersection_buf,
    h,
    c_p,
    constraints,
)
@show h, t_intersect, constraint_ind, intersection_buf.smallest_positive_roots

# assumes cs is already updated.

t_ITP, ITP_status = PowerSeriesIVP.runITP(
    intersection_buf.intersection_coefficients,
    t0,
    h,
    PowerSeriesIVP.ITPConfig(Float64),
)
@show t_ITP, ITP_status
println()

# println("timing: runITP()")
# @btime PowerSeriesIVP.runITP(
#     intersection_buf.intersection_coefficients,
#     t0,
#     h,
#     PowerSeriesIVP.ITPConfig(Float64),
# )
# println()
# 3.658 μs (1 allocation: 48 bytes)

cs = intersection_buf.intersection_coefficients

PowerSeriesIVP.evaltaylor(cs[2], t_intersect-1e-3, 0.0)
PowerSeriesIVP.evaltaylor(cs[2], t_intersect, 0.0)
PowerSeriesIVP.evaltaylor(cs[2], t_intersect+1e-3, 0.0)

x_lp = tt->collect( PowerSeriesIVP.evaltaylor(sol.coefficients[end].x[d], tt, sol.expansion_points[end]) for d = 1:N_vars )
@show norm(x_evals[end]-x_lp(t_viz[end]))

t_root = t_intersect
t_root = 0.16852762022517737
x_19 = tt->collect( PowerSeriesIVP.evaltaylor(sol.coefficients[19].x[d], tt, sol.expansion_points[19]) for d = 1:N_vars )
p = x_19(t_root+t0)
dot(as[1],p) - bs[1]
dot(as[2],p) - bs[2] # very small!!



#### show all 4-th order intersection analysis for all peices.

ts, constraint_inds = PowerSeriesIVP.searchintersection!(intersection_buf, sol, constraints)

# the two intersections.
pieces = findall(xx->xx>0, constraint_inds)
#@show pieces

t1 = ts[pieces[1]]
a1 = as[constraint_inds[pieces[1]]]

t2 = ts[pieces[2]]
a2 = as[constraint_inds[pieces[2]]]



### Budan bound.
root_ub_buf = PowerSeriesIVP.RootsUpperBoundBuffer(Float64, N_constraints, L)
bino_mat = PowerSeriesIVP.setupbinomialcoefficients(L_max)

ubs, ubs_inds = PowerSeriesIVP.upperboundintersections!(root_ub_buf, sol, constraints, bino_mat)

budan_inds = findall(xx->xx>0, ubs)
println("[budan_inds ubs[budan_inds]]")
display([budan_inds ubs[budan_inds]])
println()

@show issubset(pieces, budan_inds) # should be true if Budan's upperbound works.

Q = sol.coefficients[budan_inds]
Q[1].x

######## find root via ITP.

#PolynomialEvalParams()
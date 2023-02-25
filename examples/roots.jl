#using LinearAlgebra
#import PowerSeriesIVP
#import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1


###### get x and u.

a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2)

N_vars = 3
t_start = 0.0 # starting time.
t_fin = 25.0 # finish time.

x0 = randn(N_vars) # starting position.
u0 = randn(N_vars) # starting velocity.

# transport 2 vector fields.
N_parallel_vector_fields = 1

## set to same as u.
v0_set = collect( rand(N_vars) for _ = 1:N_parallel_vector_fields)

L_test_max = 10
N_analysis_terms = 2
L_max = L_test_max + N_analysis_terms

adaptive_order_config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L_test_max = L_test_max, # increase this for maximum higher-order power series.
    r_order = 0.1,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = N_analysis_terms,
)
Lm1_fixed = 3
fixed_order_config = PowerSeriesIVP.FixedOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L = Lm1_fixed, # actual order is L+1.
    h_zero_error = Inf,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
)

config = adaptive_order_config
#config = fixed_order_config
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, u0, v0_set)
sol = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    h_initial = 1.0
)

# @btime PowerSeriesIVP.solveIVP(
#     prob_params,
#     PowerSeriesIVP.EnableParallelTransport(),
#     # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
#     # h_initial,
#     t_start,
#     t_fin,
#     config;
#     h_initial = 1.0
# );
# @assert 1==2
# # 1.995 ms (21716 allocations: 3.17 MiB) using adaptive config.
# # 24.307 ms (419958 allocations: 29.75 MiB) using 4th order.

config2 = adaptive_order_config
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, 2 .* u0, v0_set)
sol2 = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    h_initial = 1.0
)

orders = collect( length(sol.coefficients[i].x[begin]) for i in eachindex(sol.coefficients) )
println("The min and max orders of the solution pieces:")
@show minimum(orders), maximum(orders)

println("The number of solution pieces: ", length(sol.coefficients))

#### visualize trajectory.

# set up eval points.
N_viz = 1000
t_viz = LinRange(t_start, t_fin, N_viz)

# evals.
x_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals = collect( ones(N_vars) for _ = 1:N_viz )

N_transports = PowerSeriesIVP.getNtransports(sol)
vs_evals = collect( collect( ones(N_vars) for _ = 1:N_transports ) for _ = 1:N_viz )
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, vs_evals, sol, t_viz)



y_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals2 = collect( ones(N_vars) for _ = 1:N_viz )
vs_evals2 = collect( collect( ones(N_vars) for _ = 1:N_transports ) for _ = 1:N_viz )
status_flags2 = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags2, y_evals, u_evals2, vs_evals2, sol2, t_viz)

@show norm(x_evals[1] - x0)
@show norm(y_evals[1] - x0)

# plot trajectory vs. time.

d_select = 2
x_psm_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )
y_psm_viz = collect( y_evals[n][d_select] for n in eachindex(y_evals) )

slope2 = u0[d_select]
intercept2 = x_psm_viz[begin] - slope2 * t_viz[begin]
line_geodesic = tt->( slope2*tt + intercept2)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_viz, x_psm_viz, "purple", label = "sol", linewidth = 3)
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")
PyPlot.plot(t_viz .* 2, y_psm_viz, "--", label = "sol, twice speed", linewidth = 3)

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("numerical solution, dim $d_select")






# make two hyperplanes.
N_constraints = 2
as = collect( randn(N_vars) for _ = 1:N_constraints )
bs = randn(N_constraints) .- 25

if as[1][end] > 0
    as[1] = -as[1]
    bs[1] = -bs[1]
end
if as[2][end] > 0
    as[2] = -as[2]
    bs[2] = -bs[2]
end
constraints = PowerSeriesIVP.AffineConstraints(as, bs)

x1_range = LinRange(-20, 20, 20)
x2_range = LinRange(-20, 20, 20)

m = 1
x31 = collect( (bs[m] - dot(as[m][1:2], [x1_range[i]; x2_range[j]]))/as[m][3] for i in eachindex(x1_range), j in eachindex(x2_range) )

m = 2
x32 = collect( (bs[m] - dot(as[m][1:2], [x1_range[i]; x2_range[j]]))/as[m][3] for i in eachindex(x1_range), j in eachindex(x2_range) )


# plot trajectory.
x1_viz = collect( x_evals[n][1] for n in eachindex(x_evals) )
x2_viz = collect( x_evals[n][2] for n in eachindex(x_evals) )
x3_viz = collect( x_evals[n][3] for n in eachindex(x_evals) )

line_evals = collect( x0 + t .* u0 for t in t_viz )
line_x1_viz = collect( line_evals[n][1] for n in eachindex(line_evals) )
line_x2_viz = collect( line_evals[n][2] for n in eachindex(line_evals) )
line_x3_viz = collect( line_evals[n][3] for n in eachindex(line_evals) )

y1_viz = collect( y_evals[n][1] for n in eachindex(y_evals) )
y2_viz = collect( y_evals[n][2] for n in eachindex(y_evals) )
y3_viz = collect( y_evals[n][3] for n in eachindex(y_evals) )

PyPlot.figure(fig_num)
PyPlot.subplot(111, projection="3d")
fig_num += 1


PyPlot.plot_surface(x1_range, x2_range, x31, rstride=2, edgecolors="k", cstride=2,
   cmap=PyPlot.ColorMap("gray"), alpha=0.5, linewidth=0.25,
)
PyPlot.plot_surface(x1_range, x2_range, x32, rstride=2, edgecolors="k", cstride=2,
   cmap=PyPlot.ColorMap("copper"), alpha=0.5, linewidth=0.25,
)
#
PyPlot.plot(x1_viz, x2_viz, x3_viz, label = "sol")
PyPlot.plot(line_x1_viz, line_x2_viz, line_x3_viz, label = "line")

PyPlot.plot(y1_viz, y2_viz, y3_viz, label = "sol, twice speed")

PyPlot.legend()
PyPlot.xlabel("x1")
PyPlot.ylabel("x2")
PyPlot.zlabel("x3")
PyPlot.title("trajectory and constraints")

# x1 vs x2
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(x1_viz, x2_viz, "purple", label = "sol", linewidth = 3)
PyPlot.plot(line_x1_viz, line_x2_viz, label = "line")

PyPlot.plot(y1_viz, y2_viz, "--", label = "sol, twice speed", linewidth = 3)

PyPlot.legend()
PyPlot.xlabel("x1")
PyPlot.ylabel("x2")
PyPlot.title("trajectory x1, x2")


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
piece_select = 20
c_p = sol.coefficients[piece_select].x
c_u = sol.coefficients[piece_select].u
h = sol.steps[piece_select]

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
intersection_buf = PowerSeriesIVP.IntersectionBuffer(complex_zero_tol, L, N_constraints)


t_intersect, constraint_ind = PowerSeriesIVP.refinestep!(
    intersection_buf,
    c_p,
    h,
    constraints,
)
@show h, t_intersect, constraint_ind, intersection_buf.smallest_positive_roots

ts, constraint_inds = PowerSeriesIVP.searchintersection!(intersection_buf, sol, constraints)

# the two intersections.
pieces = findall(xx->xx>0, constraint_inds)
@show pieces

t1 = ts[pieces[1]]
a1 = as[constraint_inds[pieces[1]]]

t2 = ts[pieces[2]]
a2 = as[constraint_inds[pieces[2]]]



### Budan bound.
root_ub_buf = PowerSeriesIVP.BudanIntersectionBuffers(Float64, N_constraints, L)
bino_mat = PowerSeriesIVP.setupbinomialcoefficients(L_max)

ubs, ubs_inds = PowerSeriesIVP.upperboundintersections!(root_ub_buf, sol, constraints, bino_mat)

budan_inds = findall(xx->xx>0, ubs)
println("[budan_inds ubs[budan_inds]]")
display([budan_inds ubs[budan_inds]])
println()

@show issubset(pieces, budan_inds) # should be true if Budan's upperbound works.

Q = sol.coefficients[budan_inds]
Q[1].x

# I am here. 
@assert 5==4

########## general root isolation. investigate in the future.
# implementation in Budan_bracket.jl

T = Float64
cs = Vector{Vector{T}}(undef, length(as))
root_bracket = PowerSeriesIVP.BudanIntersectionBuffers(cs)

bino_mat = PowerSeriesIVP.setupbinomialcoefficients(L_max)

# mutates cs, root_bracket
function postprocesspiece!(
    cs::Vector{Vector{T}},
    root_bracket::BudanIntersectionBuffers{T},
    c_lb::Vector{Vector{T}},
    h::T,
    as::Vector{Vector{T}},
    bs::Vector{T},
    ) where T

    # get the polynomial equation.
    PowerSeriesIVP.updateintersectionpolynomial!(cs, c_lb, as, bs)
    #root_bracket = PowerSeriesIVP.BudanIntersectionBuffers(cs)
    PowerSeriesIVP.allocatebuffer!(root_bracket, cs)

    # figure out if there is a root from x ∈ [0, h] for the current polynomials cs.
    new_h, instruction = PowerSeriesIVP.findfirstroot!(
        root_bracket,
        h,
        bino_mat;
        min_bracket_len = 1e-4,
        max_iters = 100000,
    )
end




# then write routine that resolves the current piece if there is a root.
# - solve using order 4, then Budan interval val check. if might have root again, then solve using quartic.
# Review Budan's theorem again.
# mutates c_right.
# x := t-t0, t0 is the expansion point of the Taylor polynomial c.
# c[i] is the i-th order coefficient of a Taylor polynomial, i = 0, 1, 2, ..., L.
# function getshiftedpolynomial!(
#     c_shifted::Vector{T},
#     x_shift::T,
#     c::Vector{T},
#     bino_mat::Matrix{Int},
#     ) where T
    
#     L = length(c) - 1
#     @assert size(bino_mat,1) >= L

#     reisze!(c_shifted, length(c))
#     updateshiftedpolynomial!(c_shifted, c, x_shift, zero(T), bino_buffer)

#     return nothing
# end
# then write routine that loops this check over all constraints.


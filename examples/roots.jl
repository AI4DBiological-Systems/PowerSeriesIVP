#using LinearAlgebra
#import PowerSeriesIVP
#import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

include("helpers/constraint_helpers.jl")

###### get x and u.

# pick metric parameters.
# a = 2.0
# b = 23.0
# graph of (2^2 + x^2)/(23^2 + x^2) # notch.

a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2) # unimodal.

N_vars = 3

# bound constraints.
lbs = -10.5 .* ones(N_vars)
ubs = 13.3 .* ones(N_vars)

# ## IVP simulation interval and initial conditions.
t_start = 0.0 # starting time.
t_fin = 25.0 # finish time.

x0 = randn(N_vars) # starting position.
u0 = randn(N_vars) # starting velocity.

# transport 2 vector fields.
N_parallel_vector_fields = 1

## set to same as u.
v0_set = collect( rand(N_vars) for _ = 1:N_parallel_vector_fields)

## constraints.


### configs.
L_min = 4
L_max = 21

ITP_config = PowerSeriesIVP.ITPConfig(
    Float64;
    f_tol = 1e-8,
    x_tol = 1e-15,
    k1 = 0.1,
    k2 = 0.98*(1+MathConstants.golden), # see ITP paper, equation 24.
    n0 = 0,
)

# # constraints configs & buffers.

# # ### test case: 2 hyperplane constraints.
# N_constraints = 2
# #as, bs = generateHyperplaneConstraintscase1(N_constraints)
# as, bs = generateHyperplaneConstraintscase1(Float64)
# constraints = PowerSeriesIVP.HyperplaneConstraints(as, bs)
# constraints_info = PowerSeriesIVP.ConstraintsContainer(
#     as,
#     bs;
#     complex_zero_tol = 1e-8,
#     L_min = 4,
#     L_max = 10,
#     max_divisions = 0,
#     solver_config = ITP_config,
# )

# # ### test case: 6 bound constraints (box constrains for all vars)
# as, bs = PowerSeriesIVP.boxtohyperplanes(lbs, ubs)
# constraints_info = PowerSeriesIVP.ConstraintsContainer(
#     lbs,
#     ubs;
#     complex_zero_tol = 1e-8,
#     L_min = 4,
#     L_max = 10,
#     max_divisions = 0,
#     solver_config = ITP_config,
# )
# end_time = 12.539794907619783

# ### test case: hyperplane and bound constraints.
ordering_constraints = [(3,1); (1,2)]
# ordering_constraints = Vector{Tuple{Int,Int}}(undef, 0) # should get the same result as the above.
ordering_shifts = zeros(length(ordering_constraints))
constraints_info = PowerSeriesIVP.ConstraintsContainer(
    lbs,
    ubs,
    ordering_constraints, # if the m-th entry is (i,j), then it means x[i] <= x[j] + b[m]. for the m-th cosntraint.
    ordering_shifts; # the m-th entry is bs[m].
    complex_zero_tol = 1e-8,
    L_min = 4,
    L_max = 10,
    max_divisions = 0,
    solver_config = ITP_config,
)
as = constraints_info.constraints.affine.normals
bs = constraints_info.constraints.affine.offsets

# # ### No constraints.
# constraints_info = PowerSeriesIVP.NoConstraints()
# t_fin = 9.93219966069891 # the intersection time of the 2 hyperplane test case.
# t_fin = 4.0495153834453275 # the intersection time of the 6 bound, 2 hyperplane test case.

# I am here. post processing routine to detect which constraint was hit.
# I am here. Line version of IVP. data structure, engine, etc. not using IVP, or special case of RQgeodeisc with fixed order 1.
#   The latter keeps with the theme of this framework.
# clean up the I am here's.

## step selection configs.

strategy = PowerSeriesIVP.ContinuitySecondDerivative(
    PowerSeriesIVP.GuentherWolfStep(),
)
# strategy = PowerSeriesIVP.ContinuitySecondDerivative(
#     PowerSeriesIVP.VelocityContinuityStep(),
# )
# strategy = PowerSeriesIVP.GuentherWolfStep()
# strategy = PowerSeriesIVP.VelocityContinuityStep()

step_config = PowerSeriesIVP.StepConfig(
    Float64,
    strategy;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    h_max = Inf,
    reduction_factor = 1,
    discount_factor = 0.9,        
)

adaptive_order_config = PowerSeriesIVP.AdaptOrderConfig(
    step_config;
    L_min = L_min,
    L_max = L_max, # increase this for maximum higher-order power series.
    order_increase_factor = 1.35,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
)
L_fixed = 4
fixed_order_config = PowerSeriesIVP.FixedOrderConfig(
    step_config;
    L = L_fixed,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
)

config = adaptive_order_config
#config = fixed_order_config

metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, u0, v0_set)
sol, exit_flag = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    constraints_info = constraints_info,
)
@show length(sol.coefficients)
orders = PowerSeriesIVP.getpieceorders(sol)
end_time = PowerSeriesIVP.getendtime(sol)

# (c_x, c_u, next_x, next_u, h) = ([[-0.40243248293794137, -0.9754736537655263, -0.028283214576385357, 0.008649894756430094, 0.005392191906369124], [0.8540414903329187, -0.341379056315598, -0.1582639812836908, -0.032351280660213325, 0.003273286430183644], [-0.6651248667822778, -1.0410555755312705, 0.003524793084258744, 0.014106457829175854, 0.004967531777339607]], [[-0.9754736537655263, -0.056566429152770714, 0.02594968426929028, 0.021568767625476496, 0.0009422561798294823], [-0.341379056315598, -0.3165279625673816, -0.09705384198063997, 0.013093145720734577, 0.011047314839612375], [-1.0410555755312705, 0.007049586168517488, 0.04231937348752756, 0.019870127109358426, -0.0013052033827502707]], [[-0.40707354090780734, -0.9757421559656958, -0.02815903716542786], [0.852413933441089, -0.34288700416659107, -0.15872522860720406], [-0.6700771836505584, -1.0410210801714856, 0.003726784484356198]], [[-0.9757421559656958, -0.05631807433085572, 0.026257624356730923], [-0.34288700416659107, -0.3174504572144081, -0.09686548487608504], [-1.0410210801714856, 0.007453568968712396, 0.042602766460731106]], 0.004757092965693477)
# xd = tt->PowerSeriesIVP.evaltaylorAD(c_x[1], tt, 0.0)
# xnd = tt->PowerSeriesIVP.evaltaylorAD(next_x[1], tt, 0.0)
# und = tt->PowerSeriesIVP.evaltaylorAD(next_u[1], tt, 0.0)
# ud = tt->PowerSeriesIVP.evaltaylorAD(c_u[1], tt, 0.0)
#
# # same.
# xnd_0 = ForwardDiff.derivative(xnd, 0.0)
# und(0.0)
#
# # almost same.
# ud(h)
# und(0.0)
#
# # discrepancy.
# xd_h = ForwardDiff.derivative(xd, h)
# xnd_0 = ForwardDiff.derivative(xnd, 0.0)
#@assert 1==23

# @btime sol, exit_flag = PowerSeriesIVP.solveIVP(
#     prob_params,
#     PowerSeriesIVP.EnableParallelTransport(),
#     t_start,
#     t_fin,
#     config;
#     constraints_info = constraints_info,
# );
# @assert 1==2
# # ContinuitySecondDerivative, GuentherWolfStep:
# # 3.231 ms (19791 allocations: 2.23 MiB) # using adaptive config.
# # 14.415 ms (211695 allocations: 14.98 MiB #  using 4th order.
# # with 2 hyper constraints.
# 6.237 ms (89927 allocations: 6.38 MiB) # 2 hyper constraints, 4th order fixed..
# 6.069 ms (89925 allocations: 6.38 MiB) # no constraints.
# 1.514 ms (9320 allocations: 917.77 KiB) # adaptive order, constraints. 44 pieces.
# 1.315 ms (8142 allocations: 909.23 KiB) # adaptive order, constraints. 18 pieces.

# 694.140 μs (4097 allocations: 431.34 KiB) # 9 pieces, adaptive, 6 bound & 2 hyper constraints.
# 573.573 μs (3741 allocations: 408.56 KiB) # 8 pieces. adaptive. no constraints.
# @assert 1==23





config2 = adaptive_order_config
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, 2 .* u0, v0_set)
sol2, exit_flag2 = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    constraints_info = constraints_info,
)

orders = collect( PowerSeriesIVP.getorder(sol.coefficients[i]) for i in eachindex(sol.coefficients) )
println("The min and max orders of the solution pieces:")
@show minimum(orders), maximum(orders)

println("The number of solution pieces: ", length(sol.coefficients))
#@assert 1==2

#### visualize trajectory.

# set up eval points.
N_viz = 1000
t_viz = LinRange(t_start, end_time, N_viz)

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
PyPlot.plot(t_viz .* 2, y_psm_viz, "--", label = "sol, twice speed (time *2)", linewidth = 3)

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("numerical solution, dim $d_select")


#@assert 1==23




# plot constraint hyperplanes.
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

#@assert 1==2


### next, box constraints. and mix of box and affine.
# then interface of how one caneasily specify box and affine.
# move onto RCG after that.

# I am here.
# implement combo Constraint container where
# box, general affine, and no constraints can be assigned to any combination of variables.
# implement backtrack Budan,
# write intersection tests.
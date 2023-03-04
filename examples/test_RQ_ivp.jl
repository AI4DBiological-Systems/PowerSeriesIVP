
Random.seed!(25)

include("./helpers/example_taylor_series.jl")
include("./helpers/test_helpers.jl")

PyPlot.close("all")
fig_num = 1

######## set up oracle..

N_vars = 3
θ_sin = 1.23
θ_cos = 3.42
xs = Vector{Function}(undef, N_vars)
xs[1] = tt->sin(θ_sin*tt)
xs[2] = tt->(log(tt)^2+1)
xs[3] = tt->cos(θ_cos*tt)

dxs = collect( tt->ForwardDiff.derivative(xs[d], tt) for d in eachindex(xs) )

# initial conditions.
t_start = rand()*3.2
t_fin = 100.0
#t_fin = 1000.0

# pick metric parameters.
# a = 2.0
# b = 23.0
# graph of (2^2 + x^2)/(23^2 + x^2)

a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2)

x0 = collect( xs[i](t_start) for i in eachindex(xs) )
u0 = collect( dxs[i](t_start) for i in eachindex(dxs) )

# ## transport vector fields.

# transport 2 vector fields.
N_parallel_vector_fields = 2

## set to same as u.
v0_set = collect( copy(u0) for _ = 1:N_parallel_vector_fields)

# ## set to random.
# v0_set = collect( randn(N_vars) for _ = 1:N_parallel_vector_fields )

########### test DE expansion.
T = Float64
#h_initial = convert(T, NaN) # use default strategy to get h_intial for solving each solution piece.


############ solve for piece-wise solution, then show it is C^{1} differentiable.

# println("solveIVP!")
# ϵ = 1e-6
# config = PowerSeriesIVP.IVPConfig(
#     ϵ;
#     L_test_max = 10,
#     #L_test_max = 30,
#     r_order = 0.3,
#     h_max = one(T),
#     step_reduction_factor = 2.0,
#     max_pieces = 100000000,
# )
# metric_params = PowerSeriesIVP.RQ22Metric(a,b)
# prob_params = PowerSeriesIVP.GeodesicIVPStatement(metric_params, x0, u0)
# sol = PowerSeriesIVP.solveIVP!(
#     prob_params,
#     PowerSeriesIVP.DisableParallelTransport(),
#     h_initial,
#     t_start,
#     t_fin,
#     config,
# )


##### solve IVP with parallel transport.
println("solveIVP, with parallel transport")
config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    ϵ = 1e-6,
    L_test_max = 10, # increase this for maximum higher-order power series.
    r_order = 0.3,
    h_max = one(T),
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = 2,
)
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPStatement(metric_params, x0, u0, v0_set)
sol = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    #h_initial,
    t_start,
    t_fin,
    config;
    #checkconstraintsfunc = xx->false
)

orders = collect( length(sol.coefficients[i].x[begin]) for i in eachindex(sol.coefficients) )
println("The min and max orders of the solution pieces:")
@show minimum(orders), maximum(orders)

println("The number of solution pieces: ", length(sol.coefficients))

# package up for analysis.
orders = collect( length(sol.coefficients[d].x[begin]) for d in eachindex(sol.coefficients) )
sol_fun = tt->evalsolutionvectorout(sol, tt)


###### check against DiffEq.




###### diff eqns.

# Simple Pendulum Problem
#import OrdinaryDiffEq
import DifferentialEquations

#Initial Conditions
y_initial = [x0; u0]
tspan = (t_start, t_fin)

function RQproblem(dy, y, p, t)
    
    a, b, N_vars = p
    x_range = 1:N_vars
    u_range = (N_vars+1):(N_vars*2)
    
    x = y[x_range]
    u = y[u_range]

    theta = evalθ(a, b, x, u)[begin]

    dy[x_range] .= u
    dy[u_range] .= theta

    return nothing
end

prob_toolbox = DifferentialEquations.ODEProblem(
    RQproblem,
    y_initial,
    tspan,
    (a,b,length(x0))
)

#toolbox_algorithm = DifferentialEquations.Tsit5()
toolbox_algorithm = DifferentialEquations.TsitPap8()
#toolbox_algorithm = DifferentialEquations.VCABM5()
#toolbox_algorithm = DifferentialEquations.Feagin10() # issues.
#toolbox_algorithm = DifferentialEquations.Feagin14() # issues.
sol_toolbox = DifferentialEquations.solve(
    prob_toolbox,
    toolbox_algorithm,
)

# set up eval points.
N_viz = 1000
t_viz = LinRange(t_start, t_fin, N_viz)

# evals.
y_toolbox = sol_toolbox.(t_viz)

x_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals = collect( ones(N_vars) for _ = 1:N_viz )

N_transports = PowerSeriesIVP.getNtransports(sol)
vs_evals = collect( collect( ones(N_vars) for _ = 1:N_transports ) for _ = 1:N_viz )
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, vs_evals, sol, t_viz)


# in this example script, the evaluated vector fields should be the same as u.

function sumvectorfielddiff(vs_evals, u_evals)
    discrepancy = 0.0
    for i in eachindex(vs_evals)
        for m in eachindex(vs_evals[i])
            discrepancy += norm(vs_evals[i][m] - u_evals[i])
        end
    end
    return discrepancy
end

# This should be zero since our initial conditions for all of the vector fields are the same as the initial conditions for u.
println("transport vector field evals vs. u:")
@show sumvectorfielddiff(vs_evals, u_evals)
println()

## prepare for plot.
d_select = 3
y_tb_viz = collect( y_toolbox[n][d_select] for n in eachindex(y_toolbox) )
y_psm_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )

# straight line from (t_start, x0) to (t_fin, x_sol(t_fin)).
# - use the following systems of equations.
# y_psm_viz[end] = slope * t_viz[end] + intercept
# y_psm_viz[begin] = slope * t_viz[begin] + intercept

mat = [ t_viz[end] 1;
        t_viz[begin] 1;]

obs = [y_psm_viz[end]; y_psm_viz[begin]]
tmp = mat\obs
slope, intercept = tmp
line_func = tt->( slope*tt + intercept)

# straight line from (t_start, x0) with slope u0.
# use the following to get intercept:
# y_psm_viz[begin] = slope * t_viz[begin] + intercept
slope2 = u0[d_select]
intercept2 = y_psm_viz[begin] - slope2 * t_viz[begin]
line_geodesic = tt->( slope2*tt + intercept2)

println("variable $(d_select): discrepancy over t_viz between $(toolbox_algorithm) and implemented method: ")
@show norm(y_tb_viz-y_psm_viz)
println()

# plot.
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_viz, y_tb_viz, label = "toolbox")
PyPlot.plot(t_viz, y_psm_viz, label = "sol")
PyPlot.plot(t_viz, line_func.(t_viz), "--", label = "line: solution endpoints")
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("numerical solution, dim $d_select")




##### straight line.

# the relationship between the index-normalized Taylor coefficients of 
# u := dx/dt  and x: they are the derivatives of x.
piece_select = 1
dim_select = 2
display("""
index-normalized Taylor coefficients for u := dx/dt (column 1) and x (column 2)
""")
display_mat = [collect( sol.coefficients[piece_select].u[dim_select][i] ./ i for i in eachindex(sol.coefficients[piece_select].u[dim_select]) )  sol.coefficients[piece_select].x[dim_select] ]
display( display_mat )


line = PowerSeriesIVP.createline(x0, u0, t_start, t_fin)

# evals.
x_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals = collect( ones(N_vars) for _ = 1:N_viz )
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, line, t_viz)

## prepare for plot.
d_select = 1
y_line_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )


PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_viz, y_tb_viz, label = "toolbox")
PyPlot.plot(t_viz, y_line_viz, label = "line")
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("line, dim $d_select")

println("check line against line_geodesic:")
@show norm(line_geodesic.(t_viz) - y_line_viz) # should be zero if createline() is correctly implemented.
println()



########### show C^{1} differentiable at piece boundaries.

struct EvalXTrait end
struct EvalUTrait end

# based on PowerSeriesIVP.evalsolution!().
function evalpiecewisesolAD(
    ::Type{FT},
    A::PowerSeriesIVP.PiecewiseTaylorPolynomial,
    t,
    ) where FT

    expansion_points = A.expansion_points
    t_start = PowerSeriesIVP.getstarttime(A)
    t_fin = PowerSeriesIVP.getendtime(A)

    if !(t_start <= t <= t_fin)

        return NaN
    end

    for k in Iterators.drop(eachindex(expansion_points), 1)
    
        if t < expansion_points[k]
            
            return evalpiecewisesolAD(FT, A.coefficients[k-1], t, expansion_points[k-1])
        end
    end

    if t < expansion_points[end] + A.steps[end]
        return evalpiecewisesolAD(FT, A.coefficients[end], t, expansion_points[end])
    end

    # case: our solver algorithm did not reach t_fin, and t is beyond the last solution piece's estimated interval of validity.
    return NaN
end

function evalpiecewisesolAD(::Type{EvalXTrait}, c::PowerSeriesIVP.GeodesicPiece, t, a)

    @assert length(c.x) == length(c.u)

    return collect( evaltaylorAD(c.x[d], t, a) for d in eachindex(c.x) )
end

function evalpiecewisesolAD(::Type{EvalUTrait}, c::PowerSeriesIVP.GeodesicPiece, t, a)

    @assert length(c.x) == length(c.u)

    return collect( evaltaylorAD(c.u[d], t, a) for d in eachindex(c.u) )
end

# based on horner's method.
function evaltaylorAD(c, x, a)
    
    x0 = x-a

    b = c[end]
    for n = length(c):-1:2
        b = c[n-1] + b*x0
    end

    return b
end


x_AD = tt->evalpiecewisesolAD(EvalXTrait, sol, tt)
u_AD = tt->evalpiecewisesolAD(EvalUTrait, sol, tt)

PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, vs_evals, sol, t_viz)
x_AD_viz = x_AD.(t_viz)
discrepancies_x = collect( norm(x_AD_viz[i] - x_evals[i]) for i in eachindex(x_evals) )
@assert sum(discrepancies_x) < eps(Float64)*10

d_select = 3
y_psm_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )
y_psm_AD_viz = collect( x_AD_viz[n][d_select] for n in eachindex(x_AD_viz) )

# PyPlot.figure(fig_num)
# fig_num += 1

# PyPlot.plot(t_viz, y_psm_AD_viz, label = "AD")
# PyPlot.plot(t_viz, y_psm_viz, label = "sol")

# PyPlot.legend()
# PyPlot.xlabel("time")
# PyPlot.ylabel("discrepancy")
# PyPlot.title("dim = $d_select")



dx_AD = tt->collect( ForwardDiff.derivative(xx->(x_AD(xx)[d]), tt) for d = 1:N_vars )

# sanity check. should be zero.
println("sanity check on the AD implementation: should yield the same as initial conditions")
@show norm( x_AD(t_start) - x0 )
@show norm( u_AD(t_start) - u0 )
println()

t = t_start + 0.1
println("derivative of position curve at t = $t:")
@show norm( dx_AD(t) - u_AD(t) )

t_tests = LinRange(t_start, PowerSeriesIVP.getendtime(sol) - 1e-12, 2000)
discrepancies = collect( norm( dx_AD(t_tests[i]) - u_AD(t_tests[i]) ) for i in eachindex(t_tests) )
@show sum(discrepancies)
@show findmax(discrepancies)

dx_AD(t_tests[849]) - u_AD(t_tests[849])



PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(t_tests, discrepancies)
#PyPlot.plot(t_viz, discrepancies_x)
PyPlot.plot(t_viz, y_psm_AD_viz, label = "AD")
PyPlot.plot(t_viz, y_psm_viz, label = "sol")

PyPlot.xlabel("time")
PyPlot.ylabel("discrepancy")
PyPlot.title("dx vs. u, norm discrepancy")
# It is strange that there is a spike in discrepancy near the middle.

@assert 1==2

println("t at an expansion point:")
t = sol.expansion_points[3]
@show norm( dx_AD(t) - u_AD(t) )
println()

println("t at just before an expansion point:")
t = sol.expansion_points[3]-1e-14
@show norm( dx_AD(t) - u_AD(t) )
println()

sol_c_x, sol_c_u = PowerSeriesIVP.extractcoefficients(sol)

piece_select = max(length(sol_c_x) - 3, 1)
c_x = sol_c_x[piece_select]
c_u = sol_c_u[piece_select]


var_select = 1
L = length(c_u[var_select])-1
factorial_factors = collect( factorial(l) for l = 0:L )
#scaled_u = c_u[var_select] ./ factorial_factors

cx1 = c_x[var_select]
cu1 = c_u[var_select]
scaled_u = cu1 .* factorial_factors
table = Table(x = cx1, u = cu1, scaled_u = scaled_u )
display(table)

@show cu1[1:end-1] ./ cx1[2:end]
# here. use the recursive formula for x, and the derivaitve of polynomials to show that
#   u is (probaly?) the derivative of x (see if true always?), even at piece boundary.
#  then, verify retraction and vector transport.
# for completeness, then verify against other DE solvers.
println()
println()

# select piece and test.
piece_select = 4
N_tests = 100

t = sol.expansion_points[piece_select]

discrepancies = ones(N_tests) .* Inf
t_records = ones(N_tests) .* Inf
for n = 1:N_tests
    t = convertcompactdomain(rand(), 0.0, 1.0, t0, t0+h)
    
    discrepancies[n] = norm( dx_AD(t) - u_func(t) )
    t_records[n] = t
end
@show maximum(discrepancies)
_, max_ind = findmax(discrepancies)

t_records[max_ind]

t = t0 + h
@show norm(dx_AD(t) - u_func(t))
@show dx_AD(t)
@show u_func(t)

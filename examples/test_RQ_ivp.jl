
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

########### test DE expansion.
T = Float64
#h_initial = convert(T, NaN) # use default strategy to get h_intial for solving each solution piece.
h_initial = one(T)

############ solve for piece-wise solution, then show it is C^{1} differentiable.

println("solveIVP!")
ϵ = 1e-6
config = PowerSeriesIVP.IVPConfig(
    ϵ;
    L_test_max = 10,
    #L_test_max = 30,
    r_order = 0.3,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000000,
)
prob_params = PowerSeriesIVP.RQGeodesicIVP(a, b, x0, u0)
sol = PowerSeriesIVP.solveIVP!(prob_params, h_initial, t_start, t_fin, config)

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
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, sol, t_viz)

## prepare for plot.
d_select = 1
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
[collect( curve.u[1][i] ./ i for i in eachindex(curve.u[1]) )  curve.x[1] ]


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


@assert 2==44

########### show C^{1} differentiable at piece boundaries.

struct EvalXTrait end
struct EvalUTrait end

# based on PowerSeriesIVP.evalsolution!().
function evalpiecewisesolAD(::Type{FT}, A::PowerSeriesIVP.PiecewiseTaylorPolynomial, t) where FT

    expansion_points = A.expansion_points

    if !(A.t_start <= t <= A.t_fin)

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

function evalpiecewisesolAD(::Type{EvalXTrait}, c::PowerSeriesIVP.RQGeodesicPiece, t, a)

    @assert length(c.x) == length(c.u)

    return collect( evaltaylorAD(c.x[d], t, a) for d in eachindex(c.x) )
end

function evalpiecewisesolAD(::Type{EvalUTrait}, c::PowerSeriesIVP.RQGeodesicPiece, t, a)

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

@assert 1==2

dx_AD = tt->collect( ForwardDiff.derivative(xx->(x_AD(xx)[d]), tt) for d = 1:N_vars )

# sanity check. should be zero.
@show norm( x_AD(t_start) - x0 )
@show norm( u_AD(t_start) - u0 )
println()

t = t_start + 0.1
@show norm( dx_AD(t) - u_AD(t) )

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

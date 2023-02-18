# PowerSeriesIVP
Numerical solver for a specific initial value problem ordinary differential equation problem, via piece-wise power series approximation. Based on a recent advances of the Parker-Sochacki method (Guenther 2019). We added an automated order selection strategy in addition to their adaptive step-size selection strategy.

The problem we solve is the geodesic for a particular diagonal metric: (WIP)

# Usage

Set up example problem:
```julia
# pick metric parameters.
a = 2.0
b = 23.0
# graph of (2^2 + x^2)/(23^2 + x^2)
```
or 
```julia
a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2)
```


Define initial value problem (IVP) for `N_vars` number of variables:
```julia
N_vars = 3

t_start = rand()*3.2 # starting time.
t_fin = 100.0 # finish time.

x0 = randn(N_vars) # starting position.
u0 = randn(N_vars) # starting velocity.
```

Solve:
```julia
import PowerSeriesIVP

ϵ = 1e-6
h_initial = 1.0
config = PowerSeriesIVP.IVPConfig(
    ϵ;
    L_test_max = 10, # increase this for maximum higher-order power series.
    r_order = 0.3,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
)
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, u0)
sol = PowerSeriesIVP.solveIVP!(
    prob_params,
    PowerSeriesIVP.DisableParallelTransport(),
    h_initial,
    t_start,
    t_fin,
    config,
)
println("The estimated interval of validity of the piece-wise numerical solution is from time ", t_start, " to ", PowerSeriesIVP.getendtime(sol))
println()

```
If we used `max_pieces` in our piece-wise solution of the IVP, then we terminate. Therefore, the caller of `PowerSeriesIVP.solveIVP!` should check if it simulated up to the finishing time `t_fin`.

To query:
```julia
# set up eval points.
N_viz = 1000
t_viz = LinRange(t_start, t_fin, N_viz)

x_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals = collect( ones(N_vars) for _ = 1:N_viz )
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, sol, t_viz)

```

Visualize dimension `d_select`, using `PyPlot.jl` from the Julia REPL in a terminal (i.e., not tested in a notebook environment):
```julia
import PyPlot
PyPlot.close("all") # comment this out if you don't want to close previously opened PyPlot windows.
PyPlot.matplotlib["rcParams"][:update](["font.size" => 16, "font.family" => "serif"])
fig_num = 1

## evaluate solution for the chosen dimension.
d_select = 2
y_tb_viz = collect( y_toolbox[n][d_select] for n in eachindex(y_toolbox) )
y_psm_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )

## evluate the line discribed by the initial derivative condition.
slope2 = u0[d_select]
intercept2 = y_psm_viz[begin] - slope2 * t_viz[begin]
line_geodesic = tt->( slope2*tt + intercept2)

## visualize.
PyPlot.figure(fig_num)
fig_num += 1

#PyPlot.plot(t_viz, y_tb_viz, label = "toolbox")
PyPlot.plot(t_viz, y_psm_viz, label = "sol")
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("numerical solution, dim $d_select")

```

To run a single query, do the following:
```julia
T = eltype(x0)
t = clamp(t_start + rand(), t_start, t_fin)
sol_eval = PowerSeriesIVP.RQGeodesicEvaluation(T, PowerSeriesIVP.getNvars(sol))
status_flag = PowerSeriesIVP.evalsolution!(sol_eval, sol, t)
@show sol_eval.position # the queried position.
@show sol_eval.velocity # the queried velocity.
```

To create a line without solving an IVP, then evaluate the line, do the following:
```julia
# create line
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

PyPlot.plot(t_viz, y_line_viz, label = "line")
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("line, dim $d_select")

println("check line against line_geodesic:")
@show norm(line_geodesic.(t_viz) - y_line_viz) # should be zero.
println()

```

# TODO:
- Make the test-oriented scripts in `\examples` into test sets.
- Add API documentation strings, and tutorials with visualization.
- Write latex or render an image to show the metric equation on the README.md and tutorial examples.
- Separately export the simulate IVP and simulate single piece functionalities.

# Reference
1. Guenther, Jenna, and Morgan Wolf. "An adaptive, highly accurate and efficient, parker-sochacki algorithm for numerical solutions to initial value ordinary differential equation systems." Online 12 (2019): 257-281.
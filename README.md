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

# transport 2 vector fields.
N_parallel_vector_fields = 2

## set to same as u.
v0_set = collect( copy(u0) for _ = 1:N_parallel_vector_fields)

```

Solve:
```julia
import PowerSeriesIVP

config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    ϵ = 1e-6,
    L_test_max = 10, # increase this for maximum higher-order power series.
    r_order = 0.1,
    h_max = 1.0,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = 2,
)
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPStatement(metric_params, x0, u0, v0_set)
sol = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config,
)
println("The estimated interval of validity of the piece-wise numerical solution is from time ", t_start, " to ", t_start + PowerSeriesIVP.getendtime(sol))
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


```

Visualize dimension `d_select`, using `PyPlot.jl` from the Julia REPL in a terminal (i.e., not tested in a notebook environment):
```julia
import PyPlot
PyPlot.close("all") # comment this out if you don't want to close previously opened PyPlot windows.
PyPlot.matplotlib["rcParams"][:update](["font.size" => 16, "font.family" => "serif"])
fig_num = 1

## evaluate solution for the chosen dimension.
d_select = 2
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
T = Float64
t_select = 3 # index selected for this test, from the eval positions, t_viz.
t = clamp(t_viz[t_select], t_start, t_fin)
sol_eval = PowerSeriesIVP.GeodesicEvaluation(
    T,
    PowerSeriesIVP.getdim(sol),
    PowerSeriesIVP.getNtransports(sol),
)
status_flag = PowerSeriesIVP.evalsolution!(
    sol_eval,
    PowerSeriesIVP.EnableParallelTransport(),
    sol,
    t,
)
@show norm(sol_eval.position - x_evals[t_select]) # the queried position.
@show norm(sol_eval.velocity - u_evals[t_select]) # the queried velocity.
@show norm(sol_eval.vector_fields - vs_evals[t_select]) # the queried velocity.
println()

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

For now, always assign `L_min` to `4` if there are constraints.

# TODO:
- constrained, fixed order.


# Future:
- implement backtrack initial step for current order when Budan's upper bound is inconclusive, Future release.
- impelemnt simultaneous root solve for all constraint polynomials, and eventually hone in on the interval that is closest to zero, using the ANewDsc algorithm.
- Allow `L_min` to be less than 4 for constrained simulation.

# Release check list:
- Make the test-oriented scripts in `\examples` into test sets.
- Add API documentation strings, and tutorials with visualization.
- Write latex or render an image to show the metric equation on the README.md and tutorial examples.
- Separately export the simulate IVP and simulate single piece functionalities.
- intersection tests, numerical, graphical / tutorial.

# Reference
1. Guenther, Jenna, and Morgan Wolf. "An adaptive, highly accurate and efficient, parker-sochacki algorithm for numerical solutions to initial value ordinary differential equation systems." Online 12 (2019): 257-281.
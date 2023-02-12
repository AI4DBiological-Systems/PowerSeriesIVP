# PowerSeriesIVP
Numerical solver for a specific initial value problem ordinary differential equation problem, via piece-wise power series approximation. Based on a recent advances of the Parker-Sochacki method, but with an automated order selection heurestic.

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
prob_params = PowerSeriesIVP.RQGeodesicIVP(a, b, x0, u0)
sol = PowerSeriesIVP.solveIVP!(prob_params, h_initial, t_start, t_fin, config);
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
PowerSeriesIVP.batchevalsolution!(x_evals, u_evals, sol, t_viz)

```

Visualize dimension `d_select`, using `PyPlot.jl`:
```julia
import PyPlot

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

# TODO:
- Add a gain parameter to the quadratic rational metric.
- Make the test-oriented scripts in `\examples` into test sets.
- Add API documentation strings, and tutorials with visualization.
- Write latex or render an image to show the metric equation on the README.md and tutorial examples.
- Separately export the simulate IVP and simulate single piece functionalities.
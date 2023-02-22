using LinearAlgebra
#using BenchmarkTools

import Random
Random.seed!(25)

import ForwardDiff
#import FiniteDifferences





###### get x and u.

a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2)

N_vars = 3
t_start = rand()*3.2 # starting time.
t_fin = 100.0 # finish time.

x0 = randn(N_vars) # starting position.
u0 = randn(N_vars) # starting velocity.

# transport 2 vector fields.
N_parallel_vector_fields = 1

## set to same as u.
v0_set = collect( rand(N_vars) for _ = 1:N_parallel_vector_fields)

import PowerSeriesIVP

adaptive_order_config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L_test_max = 10, # increase this for maximum higher-order power series.
    r_order = 0.1,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = 2,
)
fixed_order_config = PowerSeriesIVP.FixedOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L = 6, # increase this for maximum higher-order power series.
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

###### test up test scenario.
piece_select = 29
c_p = sol.coefficients[piece_select].x
c_u = sol.coefficients[piece_select].u

c_p_next = sol.coefficients[piece_select+1].x
c_u_next = sol.coefficients[piece_select+1].u

N_constraints = 5
as = collect( randn(N_vars) for _ = 1:N_constraints )
bs = randn(N_constraints)

constraint_select = 3
a = as[constraint_select]
b = bs[constraint_select]

##########

# get the polynomial equation.
c = sum( a[d] .* c_p[d] for d in eachindex(a) )
c[begin] += b

c_next = sum( a[d] .* c_p_next[d] for d in eachindex(a) )
c_next[begin] += b

## taylor polynomial for shifted version.

function evaltaylor(c::Vector, x, a)
    τ = x-a
    return sum( c[i]*τ^(i-1) for i in eachindex(c) )
end

L_p1 = length(c_p[begin])

t0 = sol.expansion_points[piece_select]
h = sol.steps[piece_select]

t0_next = sol.expansion_points[piece_select+1]

#evaltaylor(c, t, t0)
f = tt->evaltaylor(c, tt, t0)
f_next = tt->evaltaylor(c_next, tt, t0_next)

df = xx->ForwardDiff.derivative(f, xx)

orders = collect( length(sol.coefficients[i].x[begin]) for i in eachindex(sol.coefficients) )
println("The min and max orders of the solution pieces:")
@show minimum(orders), maximum(orders)

println("The number of solution pieces: ", length(sol.coefficients))

########## 

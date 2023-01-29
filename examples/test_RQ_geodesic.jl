
Random.seed!(25)

include("./helpers/example_taylor_series.jl")
include("./helpers/test_compositions.jl")

L = 13
max_discrepancy_tol = 1e-10

N_approximations = 10
N_tests_per_approx = 10
getafunc = xx->10*rand()

######## multivariate cases.

N_vars = 3
θ_sin = 1.23
θ_cos = 3.42
xs = Vector{Function}(undef, N_vars)
xs[1] = tt->sin(θ_sin*tt)
xs[2] = tt->(log(tt)^2+1)
xs[3] = tt->cos(θ_cos*tt)

getseqfuncs_x = Vector{Function}(undef, N_vars)
getseqfuncs_x[1] = (aa,LL)->taylorsin(aa, LL, θ_sin)
getseqfuncs_x[2] = (aa,LL)->generateseqlogexample1(aa, LL)
getseqfuncs_x[3] = (aa,LL)->taylorcos(aa, LL, θ_cos)

# the derivative of xs.
us = Vector{Function}(undef, N_vars)

# I am here. get the derivative function's coeffs by shifting xs' coeff?
getseqfuncs_u = Vector{Function}(undef, N_vars)
getseqfuncs_u[1] = (aa,LL)->taylorsin(aa, LL, θ_sin)
getseqfuncs_u[2] = (aa,LL)->generateseqlogexample1(aa, LL)
getseqfuncs_u[3] = (aa,LL)->taylorcos(aa, LL, θ_cos)





a = rand()*2.31
p = a + 0.01

###### test individual routines.

### Δ

delta = PowerSeriesIVP.InterVariableDifference(Float64, N_vars) # test target.

for l = 0:L
    x = collect( getseqfuncs[i](a, l) for i = 1:N_vars )
    PowerSeriesIVP.initializeorder!(delta, x)
end

Δfunc = collect( xx->(xs[i](xx)-xs[j](xx)) for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect(Δfunc[i,j](p) for i = 1:N_vars, j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(delta.Δ[i,j], p, a)
    for i = 1:N_vars, j = 1:N_vars )

println("Δ:")
@show norm(eval_oracle-eval_taylor)

Δsqfunc = collect( xx->(xs[i](xx)-xs[j](xx))^2 for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect(Δsqfunc[i,j](p) for i = 1:N_vars, j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(delta.Δ_sq[i,j], p, a)
    for i = 1:N_vars, j = 1:N_vars )

println("Δ squared:")
@show norm(eval_oracle-eval_taylor)


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
fs = Vector{Function}(undef, N_vars)
fs[1] = tt->sin(θ_sin*tt)
fs[2] = tt->(log(tt)^2+1)
fs[3] = tt->cos(θ_cos*tt)

gs = Vector{Function}(undef, N_vars)

getseqfuncs = Vector{Function}(undef, N_vars)
getseqfuncs[1] = (aa,LL)->taylorsin(aa, LL, θ_sin)
getseqfuncs[2] = (aa,LL)->generateseqlogexample1(aa, LL)
getseqfuncs[3] = (aa,LL)->taylorcos(aa, LL, θ_cos)

a = rand()*2.31
p = a + 0.01

###### test individual routines.

### Δ

delta = PowerSeriesIVP.InterVariableDifference(Float64, N_vars) # test target.

for l = 0:L
    x = collect( getseqfuncs[i](a, l) for i = 1:N_vars )
    PowerSeriesIVP.initializeorder!(delta, x)
end

Δfunc = collect( xx->(fs[i](xx)-fs[j](xx)) for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect(Δfunc[i,j](p) for i = 1:N_vars, j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(delta.Δ[i,j], p, a)
    for i = 1:N_vars, j = 1:N_vars )

println("Δ:")
@show norm(eval_oracle-eval_taylor)

Δsqfunc = collect( xx->(fs[i](xx)-fs[j](xx))^2 for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect(Δsqfunc[i,j](p) for i = 1:N_vars, j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(delta.Δ_sq[i,j], p, a)
    for i = 1:N_vars, j = 1:N_vars )

println("Δ squared:")
@show norm(eval_oracle-eval_taylor)

### ΔSumCol

buf_delta = PowerSeriesIVP.InterVariableDifference(Float64, N_vars)

sc = PowerSeriesIVP.ΔSumCol(Float64, N_vars) # test target.

for l = 0:L
    x = collect( getseqfuncs[i](a, l) for i = 1:N_vars )
    PowerSeriesIVP.initializeorder!(buf_delta, x)
    PowerSeriesIVP.initializeorder!(sc, buf_delta.Δ)
end

Δfunc = collect( xx->(fs[i](xx)-fs[j](xx)) for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect( sum( Δfunc[k,j](p) for k = 1:N_vars ) for j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(sc.c[i], p, a)
    for i = 1:N_vars )

println("ΔSumCol:")
@show norm(eval_oracle-eval_taylor)

### ΔSumColProduct

buf_delta = PowerSeriesIVP.InterVariableDifference(Float64, N_vars)

scp = PowerSeriesIVP.ΔSumColProduct(Float64, N_vars) # test target.

for l = 0:L
    x = collect( getseqfuncs[i](a, l) for i = 1:N_vars )
    y = x
    PowerSeriesIVP.initializeorder!(buf_delta, x)
    PowerSeriesIVP.initializeorder!(scp, buf_delta.Δ, y)
end

Δfunc = collect( xx->(fs[i](xx)-fs[j](xx)) for i = 1:N_vars, j = 1:N_vars )
eval_oracle = collect( sum( Δfunc[k,j](p)*fs[k](p) for k = 1:N_vars ) for j = 1:N_vars )
eval_taylor = collect(
    PowerSeriesIVP.evaltaylorguarded(scp.c[i], p, a)
    for i = 1:N_vars )

println("ΔSumColProduct:")
@show norm(eval_oracle-eval_taylor)


######### univariate cases.
θ_sin = 1.23
s = 3.4

f = tt->sin(θ_sin*tt)
g = tt->(log(tt)^2+1)
h = tt->s*f(tt)/g(tt)
#c_sine = taylorsin(a, L, θ_sin)

# ## manual check.
# a = rand()*2.31
# pr = PowerSeriesIVP.Product(Float64, 1)

# b1 = collect( taylorln(a, 0) for _ = 1:1 )
# b2 = collect( taylorsin(a, 0, θ_sin) for _ = 1:1 )
# PowerSeriesIVP.initializeorder!(pr, b1, b2)

# b1 = collect( taylorln(a, 1) for _ = 1:1 )
# b2 = collect( taylorsin(a, 1, θ_sin) for _ = 1:1 )
# PowerSeriesIVP.increaseorder!(pr, b1, b2)

# b1 = collect( taylorln(a, 2) for _ = 1:1 )
# b2 = collect( taylorsin(a, 2, θ_sin) for _ = 1:1 )
# PowerSeriesIVP.increaseorder!(pr, b1, b2)

# p = a + 0.005
# sol_b12 = PowerSeriesIVP.evaltaylor(b2[1], p, a)*PowerSeriesIVP.evaltaylor(b1[1], p, a)
# sol_pr = PowerSeriesIVP.evaltaylor(pr.c[1], p, a)
# @show sol_b12, sol_pr
# @show abs(sol_b12 - sol_pr)
# @assert 1==2

s_x = -4.3
ds = runtestcompositiontwoinput(
    PowerSeriesIVP.SubtractFromScaled,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(s_x*f(tt)-g(tt)),
    getafunc,
    s_x;
    ϵ_test_radius = max_discrepancy_tol*0.1/(2*abs(s_x)), # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("SubtractFromScaled:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

s_x = 4.3
s_y = -8.4
ds = runtestcompositiontwoinput(
    PowerSeriesIVP.LinearOperation,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(s_x*f(tt)+s_y*g(tt)),
    getafunc,
    s_x,
    s_y;
    ϵ_test_radius = max_discrepancy_tol*0.1/(2*max(abs(s_x),abs(s_y))), # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("LinearOperation:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol


ds = runtestcompositiontwoinput(
    PowerSeriesIVP.Addition,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(f(tt)+g(tt)),
    getafunc;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("Addition:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol


ds = runtestcompositiontwoinput(
    PowerSeriesIVP.Subtraction,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(f(tt)-g(tt)),
    getafunc;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("Subtraction:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositiontwoinput(
    PowerSeriesIVP.ScaledProduct,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(s*f(tt)*g(tt)),
    getafunc,
    s;
    ϵ_test_radius = max_discrepancy_tol*0.1/(2*s), # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("ScaledProduct:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositiontwoinput(
    PowerSeriesIVP.Product,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(f(tt)*g(tt)),
    getafunc;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, θ_sin),
    getseqfunc2 = generateseqlogexample1,
)
println("Product:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositionsingleinput(
    PowerSeriesIVP.Squared,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(g(tt))^2,
    getafunc;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc = generateseqlogexample1,
)
println("Squared:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol


ds = runtestcompositionsingleinput(
    PowerSeriesIVP.AddConstant,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(s+g(tt)),
    getafunc,
    s;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc = generateseqlogexample1,
)
println("AddConstant:")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositionsingleinput(
    PowerSeriesIVP.Reciprocal,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(1/g(tt)),
    getafunc;
    ϵ_test_radius = max_discrepancy_tol*0.1/2, # heurestic.
    getseqfunc = generateseqlogexample1,
)
println("Reciprocal")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositionsingleinput(
    PowerSeriesIVP.ScaledReciprocal,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    tt->(s/g(tt)),
    getafunc,
    s;
    ϵ_test_radius = max_discrepancy_tol*0.1/(2*s), # heurestic.
    getseqfunc = generateseqlogexample1,
)
println("ScaledReciprocal")
@show maximum(ds), maximum(ds) < max_discrepancy_tol

ds = runtestcompositiontwoinput(
    PowerSeriesIVP.ScaledQuotient,
    Float64,
    N_approximations,
    N_tests_per_approx,
    L,
    h,
    getafunc,
    s;
    ϵ_test_radius = max_discrepancy_tol*0.1/(2*s), # heurestic.
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, 1.23),
    getseqfunc2 = generateseqlogexample1,
)
println("ScaledQuotient")
@show maximum(ds), maximum(ds) < max_discrepancy_tol


#@assert 1==2

## test horner.
a = randn()
x = a + 0.1

sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
@show abs(sol-sol_horner)

@btime sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
@btime sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
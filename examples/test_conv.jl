
Random.seed!(25)

include("./helpers/example_taylor_series.jl")
include("./helpers/test_compositions.jl")

L = 13
max_discrepancy_tol = 1e-10

θ_sin = 1.23
a = rand()*2.31
#a = 0.0043433638998377775 # problem.
s = 3.4

f = tt->sin(θ_sin*tt)
g = tt->(log(tt)^2+1)
h = tt->s*f(tt)/g(tt)
#h = tt->1/g(tt)

#c_sine = taylorsin(a, L, θ_sin)


N_approximations = 300
N_tests_per_approx = 50

getafunc = xx->10*rand()
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
@show maximum(ds), maximum(ds) < max_discrepancy_tol



getafunc = xx->10*rand()
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
@show maximum(ds), maximum(ds) < max_discrepancy_tol



getafunc = xx->10*rand()
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
@show maximum(ds), maximum(ds) < max_discrepancy_tol


@assert 1==2

## test horner.
a = randn()
x = a + 0.1

sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
@show abs(sol-sol_horner)

@btime sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
@btime sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
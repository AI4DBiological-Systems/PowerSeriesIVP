
Random.seed!(25)

include("./helpers/example_taylor_series.jl")
include("./helpers/test_helpers.jl")

L = 3
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
us = collect( tt->ForwardDiff.derivative(xs[i], tt) for i in eachindex(xs) )

# the sequences for us.
getseqfuncs_u = collect(
    (aa,LL)->getderivativesequence(getseqfuncs_x[i](aa,LL+1))
    for i = 1:N_vars)


# # Test θ

# pick exapansion time t0, and query time t.
t0 = rand()*2.31
t = t0 + 0.01

# pick metric parameters.
a = 2.0
b = 1.0
x0 = collect( getseqfuncs_x[i](t0, 0)[begin] for i in eachindex(getseqfuncs_x) )
u0 = collect( getseqfuncs_u[i](t0, 0)[begin] for i in eachindex(getseqfuncs_u) )

# sanity check:
u0_oracle = collect( us[i](t0) for i in eachindex(us) )
x0_oracle = collect( xs[i](t0) for i in eachindex(xs) )
@show norm(x0-x0_oracle)
@show norm(u0-u0_oracle)

###### test θ

prob = PowerSeriesIVP.RQGeodesicBuffer(a, b, x0, u0) # test target.
theta = prob.θ

l = 0
x = collect( getseqfuncs_x[i](t0, l) for i = 1:N_vars )
u = collect( getseqfuncs_u[i](t0, l) for i = 1:N_vars )
PowerSeriesIVP.initializeorder!(theta, x, u)

for l = 1:L
    x = collect( getseqfuncs_x[i](t0, l) for i = 1:N_vars )
    u = collect( getseqfuncs_u[i](t0, l) for i = 1:N_vars )
    PowerSeriesIVP.increaseorder!(theta, x, u)
end


#### explicit eval of θ.


x_t = evalvecfunc(xs, t)
u_t = evalvecfunc(us, t)

# sanity check
l = L
x = collect( getseqfuncs_x[i](t0, l) for i = 1:N_vars )
u = collect( getseqfuncs_u[i](t0, l) for i = 1:N_vars )
x_t_taylor = evalvectaylor(x, t, t0)
u_t_taylor = evalvectaylor(u, t, t0)
@show norm(x_t - x_t_taylor)
@show norm(u_t - u_t_taylor)



# test Δ
Δ_t = evalθ(a, b, xs, us, t)[2]
Δ_t_taylor = evalvectaylor(vec(theta.delta.Δ), t, t0)
@show norm(Δ_t - Δ_t_taylor)
delta_mat = reshape(Δ_t, N_vars, N_vars)

# test Δ
Δ_sq_t = evalθ(a, b, xs, us, t)[3]
Δ_sq_t_taylor = evalvectaylor(vec(theta.delta.Δ_sq), t, t0)
@show norm(Δ_sq_t - Δ_sq_t_taylor)

# test W4
W4_t = evalθ(a, b, xs, us, t)[4]
W4_t_taylor = evalvectaylor(theta.W4.c, t, t0)
@show norm(W4_t - W4_t_taylor)

# test W8
W8_t = evalθ(a, b, xs, us, t)[5]
W8_t_taylor = evalvectaylor(theta.W8.c, t, t0)
@show norm(W8_t - W8_t_taylor)

# test A
A_t = evalθ(a, b, xs, us, t)[6]
A_t_taylor = evalvectaylor(theta.A.c, t, t0)
@show norm(A_t - A_t_taylor)

# test B
B_t = evalθ(a, b, xs, us, t)[7]
B_t_taylor = evalvectaylor(theta.B.c, t, t0)
@show norm(B_t - B_t_taylor)

# test R
R_t = evalθ(a, b, xs, us, t)[8]
R_t_taylor = evalvectaylor(theta.R.c, t, t0)
@show norm(R_t - R_t_taylor)

# test η
η_t = evalθ(a, b, xs, us, t)[9]
η_t_taylor = evalvectaylor(theta.η.c, t, t0)
@show norm(η_t - η_t_taylor)

# test C
C_t = evalθ(a, b, xs, us, t)[10]
C_t_taylor = evalvectaylor(theta.C.c, t, t0)
@show norm(C_t - C_t_taylor)

# test W9
W9_t = evalθ(a, b, xs, us, t)[11]
W9_t_taylor = evalvectaylor(theta.W9.c, t, t0)
@show norm(W9_t - W9_t_taylor)


# stage 3
W3_t = evalθ(a, b, xs, us, t)[12]
W3_t_taylor = evalvectaylor(theta.W3.c, t, t0)
@show norm(W3_t - W3_t_taylor)

W5_t = evalθ(a, b, xs, us, t)[13]
W5_t_taylor = evalvectaylor(theta.W5.c, t, t0)
@show norm(W5_t - W5_t_taylor)

W1_t = evalθ(a, b, xs, us, t)[14]
W1_t_taylor = evalvectaylor(theta.W1.c, t, t0)
@show norm(W1_t - W1_t_taylor)


#@assert 1==232

# test θ
θ_t = evalθ(a, b, xs, us, t)[1]

θ_t_taylor = evalvectaylor(theta.c, t, t0)
@show norm(θ_t - θ_t_taylor)






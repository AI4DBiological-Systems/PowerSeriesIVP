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

config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    ϵ = 1e-6,
    L_test_max = 10, # increase this for maximum higher-order power series.
    r_order = 0.1,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = 2,
)
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, u0, v0_set)
sol = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config,
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

#df(x)

# struct ShiftedPolynomialBuffer{T}
#     c::Vector{T}
#     times::Vector{T} # Taylor expansion time, shift time.

#     product_buffer::Matrix{T}
# end

function setupbinomialcoefficients(order::Integer)::Matrix{Int}
    @assert order > 0
    
    A = Matrix{Int}(undef, order, order)
    for m in axes(A, 2)
        for n in axes(A, 1)
            A[n,m] = binomial(n,m)
        end
    end
    
    return A
end


# TODO cache prod( n-i for i = 0:m-1), use Horner's method.
# x := t-t1, t1 is lower endpoint on interval.
function higherderivativespolynomial(c::Vector{T}, t, a, m) where T
    L = length(c) - 1
    x = t-a
    return sum( c[n+1] * prod( n-i for i = 0:m-1) * x^(n-m) for n = m:L )
end

# TODO pass prod() cache to higherderivativespolynomial()
function makeshiftedcoefficients(c::Vector{T}, b, a) where T
    out = similar(c)
    out[begin] = evaltaylor(c, b, a)
    for l = 1:length(out)-1
        out[l+1] = higherderivativespolynomial(c, b, a, l)/factorial(l)
    end
    return out
end

# this is a version of makeshiftedcoefficients() that uses binomial coefficients.
function updateshiftedpolynomial!(
    out::Vector{T},
    c::Vector{T},
    b::T,
    a::T,
    binomial_mat::Matrix{Int},
    ) where T

    resize!(out, length(c))
    L = length(c) - 1

    out[begin] = evaltaylor(c, b, a)

    x = b-a
    for m = 1:L
        #out[l+1] = higherderivativespolynomial(c, b, a, l)/factorial(l)

        out[m+1] = sum( c[n+1] * binomial_mat[n,m] * x^(n-m) for n = m:L )
    end
    return nothing
end

#z = t0+1.5
#z = t0+h
z = h # take value from 0 to h.
c_z0 = makeshiftedcoefficients(c, z, t0)

bino_buffer = setupbinomialcoefficients(length(c)-1)
c_z = similar(c)
updateshiftedpolynomial!(c_z, c, z, t0, bino_buffer)
@assert norm(c_z0 - c_z) < 1e-12

#t = t0
t = randn()
@show f(t + z)
@show f_next(t0_next)

# the second slot of evaltaylor is subtracted by the third slot to yield the monomital base value.
# therefore, it should be just t and zero offset for the shifted version.
g1 = tt->evaltaylor(c_z, tt, 0.0)
@show g1(t)

# show numerically, that g1(t) is h(t+z).

c_z2 = similar(c)
updateshiftedpolynomial!(c_z2, c, z, 0.0, bino_buffer)
g2 = tt->evaltaylor(c_z2, tt-t0, 0.0)

h = tt->evaltaylor(c, tt-t0, 0.0)
@show t, norm(g2(t)-h(t+z))
@show -t, norm(g2(-t)-h(-t+z))
println()

println("[c c_z c_z2 c_next]")
display([c c_z c_z2 c_next])

@assert 1==43


############# zero test.

@show PowerSeriesIVP.countsignflips(c_z)


###### continuous x, u check.
@show t0, f(t0+h)
@show t0_next, f_next(t0_next)
println("[c c_next]")
display([c c_next])

#evaltaylor(sol.coefficients[piece_select])
orders = collect( length(sol.coefficients[i].x[begin]) for i in eachindex(sol.coefficients) )
@show minimum(orders), maximum(orders)

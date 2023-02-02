
Random.seed!(25)

include("./helpers/example_taylor_series.jl")
include("./helpers/test_helpers.jl")

L = 3
max_discrepancy_tol = 1e-10

N_approximations = 10
N_tests_per_approx = 10
getafunc = xx->10*rand()

######## set up oracle..

N_vars = 3
θ_sin = 1.23
θ_cos = 3.42
xs = Vector{Function}(undef, N_vars)
xs[1] = tt->sin(θ_sin*tt)
xs[2] = tt->(log(tt)^2+1)
xs[3] = tt->cos(θ_cos*tt)

dxs = collect( tt->ForwardDiff.derivative(xs[d], tt) for d in eachindex(xs) )


# set up initial-value problem parameters.
t_start = rand()*2.31
t_fin = 10.0

# pick metric parameters.
a = 2.0
b = 23.0
x0 = collect( xs[i](t_start) for i in eachindex(xs) )
u0 = collect( dxs[i](t_start) for i in eachindex(dxs) )

########### test DE expansion.
T = Float64
#h_initial = convert(T, NaN) # use default strategy to get h_intial for solving each solution piece.
h_initial = one(T)

############ solve for piece-wise solution, then show it is C^{1} differentiable.

println("solveIVP!")
t_start = rand()*3.2
t_fin = 10.0
ϵ = 1e-6
config = PowerSeriesIVP.IVPConfig(
    ϵ;
    L_test_max = 10,
    r_order = 0.3,
    h_zero_error = Inf,
    step_reduction_factor = 2.0,
    max_pieces = 100000000,
)
prob_params = PowerSeriesIVP.RQGeodesicIVP(a, b, x0, u0)
sol = PowerSeriesIVP.solveIVP!(prob_params, h_initial, t_start, t_fin, config)

orders = collect( length(sol.coefficients[d].x[begin]) for d in eachindex(sol.coefficients) )


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

    if t < expansion_point[end] + A.steps[end]
        return evalpiecewisesolAD(FT, A.coefficients[end], t, expansion_points[end])
    end

    # case: our solver algorithm did not reach t_fin, and t is beyond the last solution piece's estimated interval of validity.
    return NaN
end

function evalpiecewisesolAD(::Type{EvalXTrait}, c::PowerSeriesIVP.RQGeodesicPiece, t, a)

    @assert length(c.x) == length(c.u)

    return collect( evaltaylorAD(c.x[d], t, a) for d in eachindex(c.x) )
end

function evalpiecewiseuAD(::Type{EvalUTrait}, c::PowerSeriesIVP.RQGeodesicPiece, t, a)

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
dx_AD = tt->collect( ForwardDiff.derivative(xx->(x_AD(xx)[d]), tt) for d in eachindex(x_AD) )

# sanity check. should be zero.
@show norm( x_AD(t_start) - x0 )
@show norm( u_AD(t_start) - u0 )
println()

t = t_start + 0.1
@show norm( dx_AD(t) - u_AD(t) )


@assert 1==232

# select piece and test.
piece_select = 4
N_tests = 100

t = t_expansion[piece_select]

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

#using LinearAlgebra
#import PowerSeriesIVP
#import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

#include("helpers/test_helpers.jl")

###### get x and u.

a = 12.24
b = 3.34
# graph of (12.24^2 + x^2)/(3.34^2 + x^2)

N_vars = 3
t_start = 0.0 # starting time.
t_fin = 25.0 # finish time.

x0 = randn(N_vars) # starting position.
u0 = randn(N_vars) # starting velocity.

# transport 2 vector fields.
N_parallel_vector_fields = 1

## set to same as u.
v0_set = collect( rand(N_vars) for _ = 1:N_parallel_vector_fields)

L_min = 4
L_max = 21
N_analysis_terms = 3
#L_max = L_max + N_analysis_terms

adaptive_order_config = PowerSeriesIVP.AdaptOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L_min = L_min,
    L_max = L_max, # increase this for maximum higher-order power series.
    r_order = 0.3, # lower means higher-order, fewer pieces. very sensitive around 0.1 to 0.3.
    h_max = 1.0,
    step_reduction_factor = 2.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    N_analysis_terms = N_analysis_terms,
)
Lm1_fixed = 4
fixed_order_config = PowerSeriesIVP.FixedOrderConfig(
    Float64;
    #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
    ϵ = 1e-6,
    L = Lm1_fixed, # actual order is L+1.
    h_max = 1.0,
    max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
)

config = adaptive_order_config
#config = fixed_order_config
metric_params = PowerSeriesIVP.RQ22Metric(a,b)
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, u0, v0_set)
constraints_info = PowerSeriesIVP.NoConstraints()
sol = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    constraints_info = constraints_info,
)

# @btime PowerSeriesIVP.solveIVP(
#     prob_params,
#     PowerSeriesIVP.EnableParallelTransport(),
#     # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
#     # h_initial,
#     t_start,
#     t_fin,
#     config;
#     h_initial = 1.0
# );
# @assert 1==2
# # 1.995 ms (21716 allocations: 3.17 MiB) using adaptive config.
# # 24.307 ms (419958 allocations: 29.75 MiB) using 4th order.

config2 = adaptive_order_config
prob_params = PowerSeriesIVP.GeodesicIVPProblem(metric_params, x0, 2 .* u0, v0_set)
sol2 = PowerSeriesIVP.solveIVP(
    prob_params,
    PowerSeriesIVP.EnableParallelTransport(),
    # PowerSeriesIVP.DisableParallelTransport(), # use this line isntead for faster computation, if don't want to parallel transport the vector fields in v0_set.
    # h_initial,
    t_start,
    t_fin,
    config;
    constraints_info = constraints_info,
)

orders = collect( PowerSeriesIVP.getorder(sol.coefficients[i]) for i in eachindex(sol.coefficients) )
println("The min and max orders of the solution pieces:")
@show minimum(orders), maximum(orders)

println("The number of solution pieces: ", length(sol.coefficients))

#### visualize trajectory.

# set up eval points.
N_viz = 1000
t_viz = LinRange(t_start, t_fin, N_viz)

# evals.
x_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals = collect( ones(N_vars) for _ = 1:N_viz )

N_transports = PowerSeriesIVP.getNtransports(sol)
vs_evals = collect( collect( ones(N_vars) for _ = 1:N_transports ) for _ = 1:N_viz )
status_flags = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags, x_evals, u_evals, vs_evals, sol, t_viz)



y_evals = collect( ones(N_vars) for _ = 1:N_viz )
u_evals2 = collect( ones(N_vars) for _ = 1:N_viz )
vs_evals2 = collect( collect( ones(N_vars) for _ = 1:N_transports ) for _ = 1:N_viz )
status_flags2 = falses(N_viz) # true for good eval by sol.
PowerSeriesIVP.batchevalsolution!(status_flags2, y_evals, u_evals2, vs_evals2, sol2, t_viz)

@show norm(x_evals[1] - x0)
@show norm(y_evals[1] - x0)

# plot trajectory vs. time.

d_select = 2
x_psm_viz = collect( x_evals[n][d_select] for n in eachindex(x_evals) )
y_psm_viz = collect( y_evals[n][d_select] for n in eachindex(y_evals) )

slope2 = u0[d_select]
intercept2 = x_psm_viz[begin] - slope2 * t_viz[begin]
line_geodesic = tt->( slope2*tt + intercept2)

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(t_viz, x_psm_viz, "purple", label = "sol", linewidth = 3)
PyPlot.plot(t_viz, line_geodesic.(t_viz), "--", label = "line geodesic")
PyPlot.plot(t_viz .* 2, y_psm_viz, "--", label = "sol, twice speed", linewidth = 3)

PyPlot.legend()
PyPlot.xlabel("time")
PyPlot.ylabel("trajectory")
PyPlot.title("numerical solution, dim $d_select")






# make two hyperplanes.
N_constraints = 2
as = collect( randn(N_vars) for _ = 1:N_constraints )
bs = randn(N_constraints) .- 25

if as[1][end] > 0
    as[1] = -as[1]
    bs[1] = -bs[1]
end
if as[2][end] > 0
    as[2] = -as[2]
    bs[2] = -bs[2]
end
constraints = PowerSeriesIVP.AffineConstraints(as, bs)

x1_range = LinRange(-20, 20, 20)
x2_range = LinRange(-20, 20, 20)

m = 1
x31 = collect( (bs[m] - dot(as[m][1:2], [x1_range[i]; x2_range[j]]))/as[m][3] for i in eachindex(x1_range), j in eachindex(x2_range) )

m = 2
x32 = collect( (bs[m] - dot(as[m][1:2], [x1_range[i]; x2_range[j]]))/as[m][3] for i in eachindex(x1_range), j in eachindex(x2_range) )


# plot trajectory.
x1_viz = collect( x_evals[n][1] for n in eachindex(x_evals) )
x2_viz = collect( x_evals[n][2] for n in eachindex(x_evals) )
x3_viz = collect( x_evals[n][3] for n in eachindex(x_evals) )

line_evals = collect( x0 + t .* u0 for t in t_viz )
line_x1_viz = collect( line_evals[n][1] for n in eachindex(line_evals) )
line_x2_viz = collect( line_evals[n][2] for n in eachindex(line_evals) )
line_x3_viz = collect( line_evals[n][3] for n in eachindex(line_evals) )

y1_viz = collect( y_evals[n][1] for n in eachindex(y_evals) )
y2_viz = collect( y_evals[n][2] for n in eachindex(y_evals) )
y3_viz = collect( y_evals[n][3] for n in eachindex(y_evals) )

PyPlot.figure(fig_num)
PyPlot.subplot(111, projection="3d")
fig_num += 1


PyPlot.plot_surface(x1_range, x2_range, x31, rstride=2, edgecolors="k", cstride=2,
   cmap=PyPlot.ColorMap("gray"), alpha=0.5, linewidth=0.25,
)
PyPlot.plot_surface(x1_range, x2_range, x32, rstride=2, edgecolors="k", cstride=2,
   cmap=PyPlot.ColorMap("copper"), alpha=0.5, linewidth=0.25,
)
#
PyPlot.plot(x1_viz, x2_viz, x3_viz, label = "sol")
PyPlot.plot(line_x1_viz, line_x2_viz, line_x3_viz, label = "line")

PyPlot.plot(y1_viz, y2_viz, y3_viz, label = "sol, twice speed")

PyPlot.legend()
PyPlot.xlabel("x1")
PyPlot.ylabel("x2")
PyPlot.zlabel("x3")
PyPlot.title("trajectory and constraints")

# x1 vs x2
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(x1_viz, x2_viz, "purple", label = "sol", linewidth = 3)
PyPlot.plot(line_x1_viz, line_x2_viz, label = "line")

PyPlot.plot(y1_viz, y2_viz, "--", label = "sol, twice speed", linewidth = 3)

PyPlot.legend()
PyPlot.xlabel("x1")
PyPlot.ylabel("x2")
PyPlot.title("trajectory x1, x2")


########### continuity check.

continuity_ϵ = eps(Float64)*10
continuity_pass_flag = PowerSeriesIVP.continuitycheck(sol; atol = continuity_ϵ)
@show continuity_pass_flag

continuity_analysis_order = 4
xu0_err, xuh_err, endpoint_err, endpoint_err_x, endpoint_err_u, interval_err, orders = PowerSeriesIVP.getderivativediscrepancies(sol, continuity_analysis_order)
@show minimum(endpoint_err), maximum(endpoint_err)

println("[orders'; endpoint_err; maximum(endpoint_err, dims=1)]:")
display( [orders'; endpoint_err; maximum(endpoint_err, dims=1)] )
# we can see our step size selection algorithm gives higher derivative continuity error for higher-order approximations.
println()

max_endpoint_err = vec(maximum(endpoint_err, dims=1)) # for each piece.
max_interval_err = vec(maximum(interval_err, dims=1)) # for each piece.
table = Table(
    piece = 1:length(orders),
    orders = orders,
    endpoint_error = max_endpoint_err,
    interval_error = max_interval_err,
    step_size = sol.steps,
)
display(table)

@show eps(Float64)
@show maximum(interval_err)
@show maximum(endpoint_err)
@show maximum(xu0_err) # zero, since PSM is based on this relation.
@show maximum(xuh_err-endpoint_err) # should be 0.
# use either xuh_err or endpoint_err. they are the same due to initial condition carry-over. xuh_err does not require starting a new piece.

max_err_x_order_set = collect( maximum(endpoint_err_x[i]) for i in eachindex(endpoint_err_x) )
min_err_x_order_set = collect( minimum(endpoint_err_x[i]) for i in eachindex(endpoint_err_x) )
@show min_err_x_order_set
@show max_err_x_order_set

max_err_u_order_set = collect( maximum(endpoint_err_u[i]) for i in eachindex(endpoint_err_u) )
min_err_u_order_set = collect( minimum(endpoint_err_u[i]) for i in eachindex(endpoint_err_u) )
@show min_err_u_order_set
@show max_err_u_order_set
println()

# I am here. devise strategy based on the derivaitve continuity to control step size
# and order increases.
# Need Lipschitz continuous on d( f ∘ R ), so this is kind of important!!
# check where in the proof we need this condition, and see if we can get away with it...

PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(orders, max_endpoint_err, "x")

PyPlot.xlabel("order")
PyPlot.ylabel("abs(dx_n - u_{n+1}) error")
PyPlot.title("error vs. order")

## verify my equation for the step inequality from notes.

d = d_select
piece_select = 3
c_x = sol.coefficients[piece_select].x
c_u = sol.coefficients[piece_select].u

c_x_next = sol.coefficients[piece_select+1].x
c_u_next = sol.coefficients[piece_select+1].u

t0 = sol.expansion_points[piece_select]
h = sol.steps[piece_select]
t_next = t0 + h

dx_coefficients = Vector{Float64}(undef, 0)
PowerSeriesIVP.differentiatepolynomial!(dx_coefficients, c_x[d])

#endpoint_err1 = abs(PowerSeriesIVP.evaltaylor(dx_coefficients, t_next, t0) - c_u_next[d][begin])
# shold be the same as endpoint_err1
xuh_err1 = abs(PowerSeriesIVP.evaltaylor(dx_coefficients, t_next, t0) - PowerSeriesIVP.evaltaylor(c_u[d], t_next, t0))
#xuh_err1 = abs(PowerSeriesIVP.evaltaylordirect(dx_coefficients, t_next, t0) - PowerSeriesIVP.evaltaylordirect(c_u[d], t_next, t0))



b = c_u[d_select]
c = c_x[d_select]

dx_h = PowerSeriesIVP.evaltaylordirect(dx_coefficients, t_next, t0)
dx_h_AN = sum( n*c[begin+n]*h^(n-1) for n = 1:length(c)-1 ) # same as above.
@show abs(dx_h - dx_h_AN)

u_h = PowerSeriesIVP.evaltaylordirect(c_u[d], t_next, t0)
u_h_AN = sum( b[begin+n]*h^n for n = 1:length(b)-1 )+b[begin] # same as above.
@show abs(u_h - u_h_AN)
@show abs( xuh_err1 - abs(dx_h - u_h) )


function geterrlastline(h, b::Vector{T}) where T
    B = zero(T)
    L = length(b)-1

    for n = 2:L
        B += b[begin+n-1]*h^(n-1) - b[begin+n]*h^n
    end
    B = B-b[begin+1]*h

    return B
end

function geterrfirstline(h, b::Vector{T}, c::Vector{T}) where T
    
    L = length(b)-1

    B = sum( ( n*c[begin+n]*h^(n-1) - b[begin+n]*h^n ) for n = 1:L ) -b[begin]
    # dx_h_AN = sum( n*c[begin+n]*h^(n-1) for n = 1:length(c)-1 )
    # u_h_AN = sum( b[begin+n]*h^n for n = 1:length(b)-1 )+b[begin]
    # B = dx_h_AN- u_h_AN

    return B
end

# don't forget that xuh_err1 has abs, so need to take abs() when comparing with geterr-first /last lines.
@show abs(geterrlastline(h, b) - geterrfirstline(h, b, c))
@show abs(abs(geterrfirstline(h, b, c)) - xuh_err1)
@show abs(abs(geterrlastline(h, b)) - xuh_err1)
L = length(b)-1
@show abs( abs(b[end]*h^L) - xuh_err1 )

println("This is a key simplification")
@show abs( -b[end]*h^L -geterrfirstline(h, b, c) )

ϵ = 1e-6
h_old = PowerSeriesIVP.stepsizeformula(ϵ, c[end], length(c)-1)

function newstepfunc(ϵ, x, order)
    return (ϵ/abs(x))^(1/order)
end
h_new = newstepfunc(ϵ, b[end], length(b)-1)

h_news = collect( newstepfunc(ϵ, b[end-m], length(b)-1-m) for m = 0:L-1)
h_orders = collect( length(b)-1-m for m = 0:L-1 )

@show h_new, geterrfirstline(h_new, b, c)
@show h_old, geterrfirstline(h_old, b, c)

f = hh->geterrfirstline(hh, b, c)
display( [h_orders f.(h_news)] )

## we arrive at the old step rule.
# julia> PowerSeriesIVP.stepsizeformula(1e-6,b[end],length(b)-1)
# 0.5162609891145944

# julia> newstepfunc(1e-6/2, b[end], length(b)-2)
# 0.5162609891145944

@assert 555==4

########## test root solving code. Cubic and quartic euqations.

c_cube = randn(3)
z1, z2, z3 = PowerSeriesIVP.solvecubicequation(c_cube)
cubefunc = tt->(c_cube[begin]+ c_cube[begin+1]*tt + c_cube[begin+2]*tt^2 + tt^3)

@show cubefunc(z1), cubefunc(z2), cubefunc(z3)
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3)
println()

"""
graph of x^4 + 3*x^3 - 2*x^2 +4.3*x -3
x == -3.85917 || x == 0.657164
{x == 0.101001 +/- 1.08292 im}
"""

c_test = [-3; 4.3; -2; 3]
z1, z2, z3, z4 = PowerSeriesIVP.solvequarticequation(c_test)
@show z1, z2, z3, z4
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3), PowerSeriesIVP.isapproxreal(z4)
println()

"""
graph of x^4 - 0.1*x^3 - 2*x^2 +0.3*x +0.1
x ≈ -1.4211 || x ≈ -0.16203 || x ≈ 0.31819 || x ≈ 1.36493
"""

c_test = [0.1; 0.3; -2; -0.1]
z1, z2, z3, z4 = PowerSeriesIVP.solvequarticequation(c_test)
@show z1, z2, z3, z4
@show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3), PowerSeriesIVP.isapproxreal(z4)
println()



############### intersection code.

###### shifted polynomial algorithm: set up test scenario.
piece_select = 17
c_p = sol.coefficients[piece_select].x
c_u = sol.coefficients[piece_select].u
h = sol.steps[piece_select]
t0 = sol.expansion_points[piece_select]

c_p_next = sol.coefficients[piece_select+1].x
c_u_next = sol.coefficients[piece_select+1].u

# N_constraints = 5
# as = collect( randn(N_vars) for _ = 1:N_constraints )
# bs = randn(N_constraints)

constraint_select = 2
a = as[constraint_select]
b = bs[constraint_select]


########## general root isolation. investigate in the future.

L = Lm1_fixed + 1
T = Float64


all_roots = Vector{Complex{T}}(undef, L)
smallest_positive_roots = Vector{T}(undef, N_constraints) # real roots.



complex_zero_tol = 1e-8
intersection_buf = PowerSeriesIVP.IntersectionBuffer(complex_zero_tol, L, N_constraints)


t_intersect, constraint_ind = PowerSeriesIVP.refinestep!(
    intersection_buf,
    h,
    c_p,
    constraints,
)
@show h, t_intersect, constraint_ind, intersection_buf.smallest_positive_roots

# assumes cs is already updated.

t_ITP, ITP_status = PowerSeriesIVP.runITP(
    intersection_buf.intersection_coefficients,
    t0,
    h,
    PowerSeriesIVP.ITPConfig(Float64),
)
@show t_ITP, ITP_status
println()

# println("timing: runITP()")
# @btime PowerSeriesIVP.runITP(
#     intersection_buf.intersection_coefficients,
#     t0,
#     h,
#     PowerSeriesIVP.ITPConfig(Float64),
# )
# println()
# 3.658 μs (1 allocation: 48 bytes)

cs = intersection_buf.intersection_coefficients

PowerSeriesIVP.evaltaylor(cs[2], t_intersect-1e-3, 0.0)
PowerSeriesIVP.evaltaylor(cs[2], t_intersect, 0.0)
PowerSeriesIVP.evaltaylor(cs[2], t_intersect+1e-3, 0.0)

x_lp = tt->collect( PowerSeriesIVP.evaltaylor(sol.coefficients[end].x[d], tt, sol.expansion_points[end]) for d = 1:N_vars )
@show norm(x_evals[end]-x_lp(t_viz[end]))

t_root = t_intersect
t_root = 0.16852762022517737
x_19 = tt->collect( PowerSeriesIVP.evaltaylor(sol.coefficients[19].x[d], tt, sol.expansion_points[19]) for d = 1:N_vars )
p = x_19(t_root+t0)
dot(as[1],p) - bs[1]
dot(as[2],p) - bs[2] # very small!!



#### show all 4-th order intersection analysis for all peices.

ts, constraint_inds = PowerSeriesIVP.searchintersection!(intersection_buf, sol, constraints)

# the two intersections.
pieces = findall(xx->xx>0, constraint_inds)
#@show pieces

t1 = ts[pieces[1]]
a1 = as[constraint_inds[pieces[1]]]

t2 = ts[pieces[2]]
a2 = as[constraint_inds[pieces[2]]]



### Budan bound.
root_ub_buf = PowerSeriesIVP.BudanIntersectionBuffers(Float64, N_constraints, L)
bino_mat = PowerSeriesIVP.setupbinomialcoefficients(L_max)

ubs, ubs_inds = PowerSeriesIVP.upperboundintersections!(root_ub_buf, sol, constraints, bino_mat)

budan_inds = findall(xx->xx>0, ubs)
println("[budan_inds ubs[budan_inds]]")
display([budan_inds ubs[budan_inds]])
println()

@show issubset(pieces, budan_inds) # should be true if Budan's upperbound works.

Q = sol.coefficients[budan_inds]
Q[1].x

######## find root via ITP.

#PolynomialEvalParams()

@assert 5==4



# mutates cs, root_bracket
function postprocessstep!(
    A::BudanIntersectionBuffers{T},
    sol_piece::GeodesicPiece{T},
    x_right::T,
    constraints::AffineConstraints{T},
    bino_mat::Matrix{Int};
    allow_reduce_step = false,
    ) where T
    
    # the left boundary is the piece's expansion point.
    cs = A.cs
    cs_left = cs
    #x_left = zero(T)

    if x_right <= zero(T)
        return :error_received_non_positive_step
    end

    if getorder(sol_piece) <= 4
        out = getexplicitstep(sol_piece)
        return out
    end

    # store the intersection polynomial equations in cs.
    updateintersectionpolynomials!(
        cs, sol_piece.x, 
        constraints.normals, constraints.offsets,
    )

    # get shifted polynomial.
    cs_right = A.cs_right
    updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

    # upper bound.
    min_x = x_right
    success_exit_tag = :continue_simulation

    for m in eachindex(cs_left)
        ub_m = upperboundroots(cs_left[m], cs_right[m])

        # when Budan's theorem doesn't work.
        if ub_m < 0
            return -one(T), :reduce_order # this case shouldn't happen. Fallback to a lower order.
        end

        # when Budan's theorem works.
        if ub_m == 1
            # there is single root for this constraint.
            x, status_flag = getroot()

            if status_flag
                if 0 < x <= x_right
                    min_x = min(x, min_x)
                    success_exit_tag = :stop_simulation
                end

                # if status_flag is true, but x is invalid, it means the root solve is sure there are no real roots.
                # this means Budan's upperbound of 0 is false, for whatever reason.
                # We request reduce order, since this shouldn't happen.
            end
            
            # It could be that the solver tolerance is too high.
            # we failed to find the root using ITP. Try another solver, or fallback to a lower order.
            return -one(T), :check_tol
        
        elseif ub_m == 0
            # for sure no intersection with this constraint on the interval.
            # don't need to do anything.

        else
            # inconclusive upper bound.
            # TODO: use a root isolation algorithm here.
            # for now, we just check the left region to see if we can keep the current order, but reduce the current simulation interval h, such that Budan's upper bound is 0.
            
            if allow_reduce_step
                return -one(T), :reduce_step
            end
            
            return -one(T), :reduce_order
        end
    end

    return min_t, success_exit_tag
end




# then write routine that resolves the current piece if there is a root.
# - solve using order 4, then Budan interval val check. if might have root again, then solve using quartic.
# Review Budan's theorem again.
# mutates c_right.
# x := t-t0, t0 is the expansion point of the Taylor polynomial c.
# c[i] is the i-th order coefficient of a Taylor polynomial, i = 0, 1, 2, ..., L.
# function getshiftedpolynomial!(
#     c_shifted::Vector{T},
#     x_shift::T,
#     c::Vector{T},
#     bino_mat::Matrix{Int},
#     ) where T
    
#     L = length(c) - 1
#     @assert size(bino_mat,1) >= L

#     reisze!(c_shifted, length(c))
#     updateshiftedpolynomial!(c_shifted, c, x_shift, zero(T), bino_buffer)

#     return nothing
# end
# then write routine that loops this check over all constraints.


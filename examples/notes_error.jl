
# run intersect.jl first.

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
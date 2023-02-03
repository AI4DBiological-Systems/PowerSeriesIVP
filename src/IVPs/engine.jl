
struct IVPConfig{T}
    ϵ::T
    L_test_max::Int
    r_order::T
    h_zero_error::T
    step_reduction_factor::T
    max_pieces::Int
end

function IVPConfig(
    ϵ::T;
    L_test_max::Int = convert(Int, 10),
    r_order = convert(T, 0.3),
    h_zero_error = convert(T, Inf),
    step_reduction_factor = convert(T, 2),
    max_pieces::Int = typemax(Int),
    ) where T

    return IVPConfig(ϵ, L_test_max, r_order, h_zero_error, step_reduction_factor, max_pieces)
end

struct RQGeodesicPiece{T}
    
    # [variable index][order index]
    x::Vector{Vector{T}} # coefficients for the solution state.
    u::Vector{Vector{T}} # coefficients for first-derivative of state.
end

function getNvars(c::RQGeodesicPiece)::Int
    return length(c.x)
end

struct PiecewiseTaylorPolynomial{T,PT}

    # [piece index]
    coefficients::Vector{PT}

    # [piece index]
    expansion_points::Vector{T}
    steps::Vector{T} # diagnostic information. Use this to determine if t_fin was actually reached by our solver algorithm for a given ODE solution.

    t_start::T
    t_fin::T
end

function getNvars(A::PiecewiseTaylorPolynomial)::Int
    return length(A.coefficients[begin].x)
end

struct RQGeodesicEvaluation{T}
    position::Vector{T}
    velocity::Vector{T}
end

function RQGeodesicEvaluation(::Type{T}, N::Integer)::RQGeodesicEvaluation{T} where T
    return RQGeodesicEvaluation(ones(T,N) .* NaN,  ones(T,N) .* NaN )
end

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::RQGeodesicEvaluation{T},
    c::RQGeodesicPiece{T},
    t,
    a,
    ) where T

    @assert length(c.x) == length(c.u) == length(out.position) == length(out.velocity)

    for d in eachindex(c.x)
        out.position[d] = evaltaylor(c.x[d], t, a)
        out.velocity[d] = evaltaylor(c.u[d], t, a)
    end

    return nothing
end

# handles the selection of a solution piece from the piece-wise solution.
function evalsolution!(
    out,
    A::PiecewiseTaylorPolynomial,
    t::T,
    ) where T

    expansion_points = A.expansion_points

    if !(A.t_start <= t <= A.t_fin)

        return false
    end

    for k in Iterators.drop(eachindex(expansion_points), 1)
    
        if t < expansion_points[k]
            
            evalsolution!(out, A.coefficients[k-1], t, expansion_points[k-1])
            return true
        end
    end

    if t < expansion_points[end] + A.steps[end]
        evalsolution!(out, A.coefficients[end], t, expansion_points[end])
        return true
    end

    # case: our solver algorithm did not reach t_fin, and t is beyond the last solution piece's estimated interval of validity.
    return false
end


function evalsolution(
    A::PiecewiseTaylorPolynomial{T,RQGeodesicPiece{T}},
    t,
    )::RQGeodesicEvaluation{T} where T

    out = RQGeodesicEvaluation(T, getNvars(A))
    evalsolution!(out, A, t)

    return out
end

function batchevalsolution!(
    positions_buffer::Vector{Vector{T}},
    velocities_buffer::Vector{Vector{T}},
    A::PiecewiseTaylorPolynomial{T,RQGeodesicPiece{T}},
    ts,
    ) where T

    @assert length(positions_buffer) == length(ts) == length(velocities_buffer)
    out = RQGeodesicEvaluation(T, getNvars(A))

    for n in eachindex(ts)
        evalsolution!(out, A, ts[n])
        positions_buffer[n][:] = out.position
        velocities_buffer[n][:] = out.velocity
    end

    return nothing
end

###################


struct RQGeodesicIVP{T}
    a::T
    b::T
    x0::Vector{T}
    u0::Vector{T}
end

# generate a piece-wise polynomials (seperately for each variable) that approximately solve an IVP of the form:
# - starts at t = 0, stop at t = t_fin,
# - h_initial used for adaption of the first polynomial.
# Subsequent h_initials are based on the solved step size for the previous polynomial segment.
function solveIVP!(
    prob_params::RQGeodesicIVP{T},
    h_initial::T,
    t_start::T,
    t_fin::T,
    config::IVPConfig{T},
    ) where T

    # set up.
    t_expansion = t_start

    a = prob_params.a
    b = prob_params.b
    prob = RQGeodesicBuffer(a, b, prob_params.x0, prob_params.u0)

    sol = PiecewiseTaylorPolynomial(
        Vector{RQGeodesicPiece{T}}(undef,0),
        Vector{T}(undef,0),
        Vector{T}(undef,0),
        t_start,
        t_fin,
    )

    N_vars = length(prob.x0)
    next_conditions = RQGeodesicEvaluation(T, N_vars)

    # solve for the first solution piece.
    t_next = solvesegmentIVP!(
        sol,
        next_conditions,
        prob,
        t_expansion,
        h_initial,
        config,
    )
    
    # check stopping condition
    if t_next > t_fin
        return sol
    end

    # each iteration is a piece.
    for _ = 2:config.max_pieces 

        # set up new IVP problem for the current expansion time.
        x0_current = next_conditions.position
        u0_current = next_conditions.velocity
        prob_current = RQGeodesicBuffer(a, b, x0_current, u0_current)
        t_expansion = t_next

        # solve for the current solution piece.
        t_next = solvesegmentIVP!(
            sol,
            next_conditions,
            prob_current,
            t_expansion,
            #sol.steps[end], # use the step from the last solution piece as the initial segment. actually, might not make sense for the geodesic setting, where extreme accuracy might not be needed.
            h_initial, # more direct control over how large the step sizes could be. set this to be large to encourage taking a big step and high-order.
            config,
        )
    
        # check stopping condition.
        if t_next > t_fin
            return sol
        end
    
    end

    return sol
end


# mutates sol and next_conditions.
function storesolutionpiece!(
    sol::PiecewiseTaylorPolynomial,
    next_conditions::RQGeodesicEvaluation{T},
    prob::RQGeodesicBuffer{T},
    t_expansion::T,
    h::T,
    ) where T

    # add the coefficients for the solution piece.
    new_coefficients = RQGeodesicPiece(prob.x.c, prob.u.c)
    push!(sol.coefficients, new_coefficients)
    push!(sol.expansion_points, t_expansion)
    push!(sol.steps, h)

    # get the initial conditions for the next IVP that the next solution piece solves.
    t_next = t_expansion + h
    evalsolution!(next_conditions, new_coefficients, t_next, t_expansion)

    return nothing
end

function solvesegmentIVP!(
    sol::PiecewiseTaylorPolynomial,
    next_conditions::RQGeodesicEvaluation{T},
    prob::RQGeodesicBuffer{T},
    t_expansion::T,
    h_initial,
    config,
    ) where T

    # solve first one.
    h = computetaylorsolution!(
        prob;
        ϵ = config.ϵ,
        h_initial = h_initial,
        L_test_max = config.L_test_max,
        r_order = config.r_order,
        h_zero_error = config.h_zero_error,
    )
    h = h/config.step_reduction_factor

    # update solution, and prepare for the next IVP.
    t_next = t_expansion + h
    storesolutionpiece!(sol, next_conditions, prob, t_expansion, h)
    
    return t_next
end

# move to geodesic.
function vectortransport()
    #
end
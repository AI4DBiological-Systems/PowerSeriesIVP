
######### evaluation of line geodesic.

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::GeodesicEvaluation,
    ::DisableParallelTransport,
    c::GeodesicLine,
    t,
    a,
    )

    @assert length(c.x) == length(c.u) == length(out.position) == length(out.velocity)

    for d in eachindex(c.x)
        out.position[d] = evalline(c.x[d], t-a, c.u[d])
    end
    out.velocity[:] = c.u

    return nothing
end

# no checking against interval of validity here. That responsibility is on the calling routine.
function evalsolution!(
    out::GeodesicEvaluation,
    ::EnableParallelTransport,
    c::GeodesicLine,
    t,
    a,
    )

    evalsolution!(out, DisableParallelTransport(), c, t, a)

    for m in eachindex(c.vs)
        out.vector_fields[m][:] = c.vs[m]
    end

    return nothing
end

function evalsolution!(
    out::GeodesicEvaluation,
    pt_trait::ParallelTransportTrait,
    sol::PiecewiseTaylorPolynomial{T, GeodesicLine{T}},
    t,
    a,
    ) where T

    return evalsolution!(out, pt_trait, sol.coefficients[begin], t, a)
end

function evalline(p::T, t::T, u::T)::T where T
    return p + t*u
end

function intersectline(
    C::HyperplaneConstraints,
    line::GeodesicLine{T},
    )::T where T

    as = C.normals
    bs = C.offsets

    p = line.x
    u = line.u

    min_t = convert(T, Inf)
    for m in eachindex(as)
        t = (bs[m]-dot(as[m], p))/dot(as[m], u)
        
        if t >= zero(T)
            # @show t, as[m], bs[m], t .* u + p
            min_t = min(t, min_t)
        end
    end

    return min_t
end

########### intersection with constraints.

function intersectline(
    C::BoundConstraints,
    line::GeodesicLine{T},
    )::T where T

    lbs = C.lbs
    lb_dims = C.lb_dims
    ubs = C.ubs
    ub_dims = C.ub_dims

    p = line.x
    u = line.u

    min_t = convert(T, Inf)

    # lower bounds.
    for m in eachindex(lbs)
        
        d = lb_dims[m]
        L = lbs[m]

        t = (L-p[d])/u[d]
        
        if t >= zero(T)
            #@show t, L, t .* u + p
            min_t = min(t, min_t)
        end
    end

    # upper bounds.
    for m in eachindex(ubs)
        
        d = ub_dims[m]
        U = ubs[m]

        t = (U-p[d])/u[d]
        
        if t >= zero(T)
            #@show t, U, t .* u + p
            min_t = min(t, min_t)
        end
    end

    return min_t
end

function intersectline(
    C::AffineConstraints,
    line::GeodesicLine{T},
    )::T where T

    t1 = intersectline(C.hyperplane, line)
    t2 = intersectline(C.bound, line)

    return min(t1,t2)
end

function intersectline(
    C::ConstraintsContainer,
    line::GeodesicLine{T},
    )::T where T

    return intersectline(C.constraints, line)
end

function intersectline(
    ::NoConstraints,
    line::GeodesicLine{T},
    )::T where T

    return convert(T, Inf)
end

############# solve Euclidean metric IVP, which gives a line geodesic.

function assembleline(
    ::DisableParallelTransport,
    p::GeodesicIVPStatement{EuclideanMetric,T},
    )::GeodesicLine{T} where T
    
    return GeodesicLine(p.x0, p.u0, empty(p.v0_set))
end

function assembleline(
    ::EnableParallelTransport,
    p::GeodesicIVPStatement{EuclideanMetric,T},
    )::GeodesicLine{T} where T
    
    return GeodesicLine(p.x0, p.u0, p.v0_set)
end

function solveIVP!(
    sol::PiecewiseTaylorPolynomial{T,GeodesicLine{T}},
    ::VariableContainer,
    ::LineIVPTrait,
    prob_params::GeodesicIVPStatement{EuclideanMetric,T},
    pt_trait::ParallelTransportTrait,
    t_start::T,
    t_fin::T,
    config::IVPConfig,
    constraints_info::ConstraintsTrait,
    )::Symbol where T

    #sol_line = PowerSeriesIVP.PiecewiseTaylorPolynomial(T, getsolpiecetype(prob_params))
    
    line = assembleline(pt_trait, prob_params)
    push!(sol.coefficients, line)
    push!(sol.expansion_points, t_start)

    # check intersections.
    s_intersect = intersectline(constraints_info, line)

    # determine exit flag.
    if isfinite(s_intersect)
        if zero(T) < s_intersect < (t_fin-t_start) # s is normalized to start at 0, t starts at t_start.
            push!(sol.steps, t_start + s_intersect)

            return :intersection_found
        end
    end

    # case: no itnersection within the interval [t_start, t_fin]
    push!(sol.steps, t_fin)

    return :success
end
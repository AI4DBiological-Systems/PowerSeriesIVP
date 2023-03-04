####  front end methods for use with other parts of this library.

# ::Type{BT} is for use with trait function getbuffertype(), which is defined in Levi-Civita metric scripts like RQ22.jl
function getIVPbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    )::BT where BT <: GeodesicIVPBuffer
    
    return getIVPbuffer(
        prob_params.metric_params,
        prob_params.x0,
        prob_params.u0,
        prob_params.v0_set,
    )
end

function getIVPbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    initial_conditions::GeodesicEvaluation,
    )::BT where BT <: GeodesicIVPBuffer

    return getIVPbuffer(
        prob_params.metric_params,
        # copy(initial_conditions.position),
        # copy(initial_conditions.velocity),
        # collect( copy(initial_conditions.vector_fields[i]) for i in eachindex(initial_conditions.vector_fields) )
        initial_conditions.position,
        initial_conditions.velocity,
        initial_conditions.vector_fields,
    )
end

function copyIVPtestbuffer(
    ::Type{BT},
    prob_params::GeodesicIVPStatement,
    )::BT where BT <: GeodesicIVPBuffer
    
    return getIVPbuffer(
        prob_params.metric_params,
        copy(prob_params.x0),
        copy(prob_params.u0),
        empty(prob_params.v0_set), # no parallel transport. test is on x, which only needs x and u to be updated.
    )
end

# used only for continuity test on x, which doesn't need the computation involving the parallel transported vector field solutions.
# therefore, we ignore resetting the parallel transport field to reduce computation.
function resetbuffer!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    x0::Vector{T},
    u0::Vector{T},
    ) where T <: AbstractFloat

    resetbuffer!(p.θ)
    resetbuffer!(p.u)
    resetbuffer!(p.x)

    p.x0[:] = x0
    p.u0[:] = u0

    return nothing
end

############### GeodesicIVPBuffer methods.

#####  main routine for computing the power series method iterations for the generic geodesic IVP problem.

function getorder(p::GeodesicIVPBuffer)::Int
    return length(p.x.c[begin])-1
end

function getfirstorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    )

    # variables
    initializeorder!(p.x, p.x0)
    initializeorder!(p.u, p.u0)

    # du/dt
    initializeorder!(p.θ, p.x.c, p.u.c)

    # variables: first update.
    increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
    increaseorder!(p.u, p.θ.c)
    
    return nothing
end


function increaseorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    )

    # du/dt
    increaseorder!(p.θ, p.x.c, p.u.c)

    # variables
    increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
    increaseorder!(p.u, p.θ.c)

    return nothing
end

# only decrease x and u. Does not decrease intermediates variables and their buffers, such as  θ.
function decreaseorder!(
    p::GeodesicIVPBuffer,
    ::DisableParallelTransport,
    min_order::Integer,
    )

    # if length(p.x.c) < 1
    #     println("Warning, decreaseorder!() received an polynomial less than order 1. Did not decrease order further.")
    #     return nothing
    # end

    # variables
    decreaseorder!(p.x, 1, min_order)
    decreaseorder!(p.u, 1, min_order)

    return nothing
end

function getfirstorder!(
    p::GeodesicIVPBuffer,
    ::EnableParallelTransport,
    )

    pt = p.parallel_transport

    # set up the geodesic variables and related quantities.
    getfirstorder!(p, DisableParallelTransport())

    # set up the transport vector variables and quantities.
    for m in eachindex(pt)
        
        # variables
        initializeorder!(pt[m].v, pt[m].v0)

        # du/dt
        initializeorder!(pt[m].ζ, p.θ, pt[m].v.c)

        # variables: first update.
        increaseorder!(pt[m].v, pt[m].ζ.c)

        # increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
        # increaseorder!(p.u, p.θ.c)
    end

    return nothing
end


function increaseorder!(
    p::GeodesicIVPBuffer,
    ::EnableParallelTransport,
    )

    pt = p.parallel_transport

    # set up the geodesic variables and related quantities.
    increaseorder!(p, DisableParallelTransport())

    # update the transport vector variables and quantities.
    for m in eachindex(pt)
    
        # dv/dt
        increaseorder!(pt[m].ζ, p.θ, pt[m].v.c)

        # variables: first update.
        increaseorder!(pt[m].v, pt[m].ζ.c)

        # increaseorder!(p.x, p.u.c) # x must be updated before u, since we hardcoded updates to always use the last element.
        # increaseorder!(p.u, p.θ.c)
    end
    
    return nothing
end

# only decrease x and u and {pt[m].v}_m. Does not decrease intermediates variables and their buffers, such as  θ and ζ.
function decreaseorder!(
    p::GeodesicIVPBuffer, 
    ::EnableParallelTransport,
    min_order::Integer,
    )

    # if length(p.x.c) < 1
    #     println("Warning, decreaseorder!() received an polynomial less than order 1. Did not decrease order further.")
    #     return nothing
    # end

    pt = p.parallel_transport

    # decrease the position p.x and geodesic velocity/vector field p.u.
    decreaseorder!(p, DisableParallelTransport(), min_order)

    # update the transport vector variables and quantities.
    for m in eachindex(pt)
        decreaseorder!(pt[m].v, 1, min_order)
    end

    return nothing
end
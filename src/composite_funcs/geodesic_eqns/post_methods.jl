# "post" means include this file as the last file in this directory, in the library module.

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
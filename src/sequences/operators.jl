## the work horse behind the elemntary operations in elementary.jl

# does s/f, where b is the sequence for function f, s is a constant.
# this is getquotientseries() with c being the array with zeros everywhere except the first entry, with a value of 1.
function getreciprocalseries(b::Vector{T}, L::Integer) where T
    @assert length(b) >= L + 1

    α = zeros(T, 1)
    α[begin] = one(T)

    #b0_powers = collect( b[begin]^l for l = 0:L+1 )

    b_tilde = Vector{T}(undef, 0)
    b0 = b[begin]

    q = Vector{T}(undef, L+1) # output.
    q[begin] = one(T)/b0

    for k = 1:L
        # update.
        
        push!(b_tilde, b0^(k-1) * b[begin+k])
        #push!(b_tilde, b0_powers[begin+k-1] * b[begin+k])

        pushfirst!(α, -dot(α, b_tilde))

        # compute coefficient.
        q[begin+k] = α[begin]/b0^(k+1)
        #q[begin+k] = α[begin]/b0_powers[begin+k+1]
    end

    return q
end

function appendreciprocalcoeff!(q::Vector{T}, b::Vector{T}, b_tilde::Vector{T}, α::Vector{T}) where T
    
    @assert length(q) == length(b) - 1
    @assert length(b_tilde) == length(b) - 2
    @assert length(α) == length(b_tilde) + 1

    k = length(b) -1 # the order of the cofficient to be appened to q.

    b0 = b[begin]
    push!(b_tilde, b0^(k-1) * b[end]) # see equation for α(3,0) for an example from notes.
    pushfirst!(α, -dot(α, b_tilde))
    
    tmp = α[begin]/b0^(k+1)
    push!(q, tmp)

    return nothing
end

# use when length(b) > 400 (machine dependent).
function conv2(b::Vector{T}, c::Vector{T})::T where T

    reverse!(b)
    out = dot(b, c)
    reverse!(b)

    return out
end

function conv(b::Vector{T}, c::Vector{T})::T where T
    #@assert N <= length(b)
    #@assert N <= length(c)

    # sum( b[begin+i-1]*c[end-i+1] for i = 1:length(b) )

    # explicit loop instead of sum because b might be empty. In which case, we want to return zero.
    out = zero(T)
    for i in eachindex(b)
        out += b[begin+i-1]*c[end-i+1]
    end
    
    return out
end
## benchmark:
# N = 40
# b = randn(N)
# c = randn(N)
# @btime conv($b,$c);
# @btime conv2($b, $c);
# when N = 40
# 40.125 ns (0 allocations: 0 bytes)
# 65.454 ns (0 allocations: 0 bytes)
#
# when N = 400
# 378.598 ns (0 allocations: 0 bytes)
# 348.093 ns (0 allocations: 0 bytes)

function conv(b::Vector{T})::T where T
    #@assert N <= length(b)
    #@assert N <= length(c)

    #out = sum( b[begin+i-1]*b[end-i+1] for i = 1:length(b) )

    # explicit loop instead of sum because b might be empty. In which case, we want to return zero.
    out = zero(T)
    for i in eachindex(b)
        out += b[begin+i-1]*b[end-i+1]
    end

    return out
end
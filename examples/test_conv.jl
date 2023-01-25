
Random.seed!(25)

function conv2(b::Vector{T}, c::Vector{T})::T where T

    reverse!(b)
    out = dot(b, c)
    reverse!(b)

    return out
end

function conv(b::Vector{T}, c::Vector{T})::T where T
    #@assert N <= length(b)
    #@assert N <= length(c)
    return sum( b[begin+i-1]*c[end-i+1] for i = 1:length(b) )
end

function conv(b::Vector{T})::T where T
    #@assert N <= length(b)
    #@assert N <= length(c)
    return sum( b[begin+i-1]*b[end-i+1] for i = 1:length(b) )
end

N = 400
b = randn(N)
c = randn(N)

# @btime conv($b,$c);
# @btime conv2($b, $c);
# @btime conv($b);


N = 40
b = randn(N)
c = randn(N)

# @btime conv($b,$c);
# @btime conv2($b, $c);
# @btime conv($b);

# 378.402 ns (0 allocations: 0 bytes)
# 346.553 ns (0 allocations: 0 bytes)
# 376.980 ns (0 allocations: 0 bytes)
# 40.125 ns (0 allocations: 0 bytes)
# 56.080 ns (0 allocations: 0 bytes)
# 39.811 ns (0 allocations: 0 bytes)

##### correctness of appendreciprocalcoeff!() as an incremental version of getreciprocalseries().
# add test for reciprocal vs. actual reciprocal. then package as test.

c = randn(9)
q = PowerSeriesIVP.getreciprocalseries(c, length(c)-1)

q2 = ones(1)
q2[begin] = 1/c[begin]
b_tilde = Vector{Float64}(undef, 0)
α = ones(Float64,1)

for k = 1:length(c)-1
    PowerSeriesIVP.appendreciprocalcoeff!(q2, c[1:k+1], b_tilde, α)
end
@show norm(q-q2)

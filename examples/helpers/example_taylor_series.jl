# truncated taylor series for some functionss, expanded at x = a, order L.

function taylorexp(a::T, L::Integer, λ)::Vector{T} where T

    # (λ^0) * exp(λ*a)/factorial(0)
    # (λ^1) * exp(λ*a)/factorial(1)

    return collect( (λ^l) * exp(λ*a)/factorial(l) for l = 0:L )
end

# https://www.wolframalpha.com/input?i=n-th+derivative+of+cos%28cx%29
function taylorcos(a::T, L::Integer, c)::Vector{T} where T


    return collect( (c^l) *cos(l/2 * π + c*a)/factorial(l) for l = 0:L )
end

function taylorsin(a::T, L::Integer, c)::Vector{T} where T


    return collect( (c^l) *sin(l/2 * π + c*a)/factorial(l) for l = 0:L )
end

# https://www.wolframalpha.com/input?i=n-th+derivative+of+log%28x%29
function taylorln(a::T, L::Integer)::Vector{T} where T

    out = Vector{T}(undef, L+1)
    
    out[begin] = log(a)
    
    for l = 1:L    
        #out[begin+l] = ( (-1)^(1+l) * factorial(l-1)/a^l )/factorial(l)
        out[begin+l] = ( ((-1)^(1+l))/a^l )/l
    end

    return out
end




############ compositions.

# generates L-th order Taylor polynomial of g = tt->(log(tt)^2+1) at t = a.
function generateseqlogexample1(a::T, L) where T
    b = taylorln(a, L)
    b = collect( dot(b[begin:i],reverse(b[begin:i])) for i in eachindex(b) ) # conv(b,b)
    b[begin] += 1

    return b
end


# use Polynomial.jl (possible large indirect dependencies) or implement horner's method to eval polynomials.
function evaltaylordirect(c::Vector{T}, x::T, a)::T where T
    τ = x-a
    return sum( c[i]*τ^(i-1) for i in eachindex(c) )
end

# Horner's method for evaluating a Taylor polynomial at expansion point `a` and coefficient `c`.
function evaltaylor(c::Vector{T}, x::T, a::T)::T where T
    
    x0 = x-a

    b = c[end]
    for n = length(c):-1:2
        b = c[n-1] + b*x0
    end

    return b
end

#### adaptive

# Univariate search to find a large step `α` from expansion center `a`, such that the discrepancy between the Taylor series `q` and `q_analysis`  is less than tolerance `ϵ`.
# The stopping condition is checked only in the positive time direction.
# Returns (evaluation of Taylor series at p, α), where p = a + α. 
# Returns `NaN` if no step satisfying tolerance `ϵ` after `max_iters` iterations.
function getacceptanceradius(
    q::Vector{T},
    q_analysis::Vector{T},
    a::T,
    ϵ;
    α::T = one(T),
    max_iters::Integer = 100,
    discount_rate = convert(T,0.9),
    )::Tuple{T,T} where T
    
    @assert α > zero(T)
    @assert length(q) < length(q_analysis)
    @assert zero(T) < discount_rate < one(T)

    for _ = 1:max_iters

        p = a + α # positive time direction.

        eval_p_analysis = evaltaylor(q_analysis, p, a)
        eval_p = evaltaylor(q, p, a)

        discrepancy = abs(eval_p_analysis-eval_p)
        if discrepancy < ϵ
            return eval_p, α
        end

        α = α*discount_rate
    end

    return convert(T, NaN), convert(T, NaN)
end

function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: Real

    return (x-a)*(d-c)/(b-a)+c
end

#### generic

function generatecasesetup(
    a::T,
    L::Integer,
    N_tests::Integer;
    getseqfunc = generateseqlogexample1,
    ϵ_test_radius = 1e-10,
    α = one(T),
    max_iters = 100,
    discount_rate = convert(T, 0.9),
    )::Tuple{Vector{T},Vector{T}} where T

    @assert L > 0

    c = getseqfunc(a, L)
    c_analysis = getseqfunc(a, L+1)
    _, max_step_len = PowerSeriesIVP.getacceptanceradius(
        c,
        c_analysis,
        a,
        ϵ_test_radius;
        α = α,
        max_iters = max_iters,
        discount_rate = discount_rate,
    )
    if isfinite(max_step_len)

        return c, collect( a + max_step_len*rand() for _ = 1:N_tests )
    end

    return c, Vector{T}(undef, 0)
end


### for composition routines that take one sequence as input.

function runtestcompositionsingleinput(
    RT,
    ::Type{T},
    N_approximations,
    N_tests_per_approx,
    L::Integer,
    g::Function,
    getafunc,
    args...;
    ϵ_test_radius = 1e-10,
    getseqfunc = generateseqlogexample1,
    #testcompositionfunc = testReciprocal,
    ) where T

    discrepancy = Vector{T}(undef, N_approximations)
    fill!(discrepancy, Inf)

    for n = 1:N_approximations
        
        a::T = getafunc(0)
        b, Xs = generatecasesetup(a, L, N_tests_per_approx;
            getseqfunc = getseqfunc, ϵ_test_radius = ϵ_test_radius)

        if !isempty(Xs)
            #discrepancy[n] = testcompositionfunc(Xs, b, a, g)
            discrepancy[n] = testcompositionroutine(RT, Xs, b, a, g, args...)
            
        end
    end

    return discrepancy # Inf entries if no test points were bound to satisfy the ϵ tolerance.
end

function testcompositionroutine(RT, Xs::Vector{T}, b::Vector{T}, a::T, f, args...) where T

    N_tests = length(Xs)
    L = length(b) - 1 # max order of the Taylor polynomial stored in b.

    ### test Reciprocal.
    N = 1 # single DE case.
    r = RT(T, N)
    #r = PowerSeriesIVP.Reciprocal(T, N)

    b_pkg = collect( Vector{T}(undef,0) for _ = 1:N )

    push!(b_pkg[begin], b[begin])
    PowerSeriesIVP.initializeorder!(r, b_pkg, args...)

    for l = 1:L
        push!(b_pkg[begin], b[begin+l])
        PowerSeriesIVP.increaseorder!(r, b_pkg, args...)
    end

    discrepancy = Vector{T}(undef, N_tests)
    for n = 1:N_tests
        p = Xs[n]
        discrepancy[n] = abs( f(p) - PowerSeriesIVP.evaltaylor(r.c[begin], p, a) )
    end

    return maximum(discrepancy)
end


##### for composition routines that take two input sequences.

function runtestcompositiontwoinput(
    RT,
    ::Type{T},
    N_approximations,
    N_tests_per_approx,
    L::Integer,
    g::Function,
    getafunc,
    args...;
    ϵ_test_radius = 1e-10,
    getseqfunc1 = (aa,LL)->taylorsin(aa, LL, 1.23),
    getseqfunc2 = generateseqlogexample1,
    ) where T

    discrepancy = Vector{T}(undef, N_approximations)
    fill!(discrepancy, Inf)

    for n = 1:N_approximations
        
        a::T = getafunc(0)
        b1, Xs1 = generatecasesetup(a, L, N_tests_per_approx;
            getseqfunc = getseqfunc1, ϵ_test_radius = ϵ_test_radius)
        #
        b2, Xs2 = generatecasesetup(a, L, N_tests_per_approx;
            getseqfunc = getseqfunc2, ϵ_test_radius = ϵ_test_radius)
        
        
        if !isempty(Xs1) && !isempty(Xs2)
            Xs = collect( min(Xs1[i], Xs2[i]) for i in eachindex(Xs1) )

            discrepancy[n] = testcompositionroutine2(RT, Xs, b1, b2, a, g, args...)
        end
    end

    return discrepancy # Inf entries if no test points were bound to satisfy the ϵ tolerance.
end

function testcompositionroutine2(RT, Xs::Vector{T}, b1::Vector{T}, b2::Vector{T}, a::T, f, args...) where T

    N_tests = length(Xs)
    L = length(b2) - 1 # max order of the Taylor polynomial stored in b.
    @assert length(b1) == length(b2)

    # set up composition routine data structure.
    N = 1 # single DE case.
    r = RT(T, N)

    # push as we increment the order, like how the DE solver would utilize these compositions.
    b1_pkg = collect( Vector{T}(undef,0) for _ = 1:N )
    b2_pkg = collect( Vector{T}(undef,0) for _ = 1:N )

    push!(b1_pkg[begin], b1[begin])
    push!(b2_pkg[begin], b2[begin])
    PowerSeriesIVP.initializeorder!(r, b1_pkg, b2_pkg, args...)

    for l = 1:L
        push!(b1_pkg[begin], b1[begin+l])
        push!(b2_pkg[begin], b2[begin+l])
        PowerSeriesIVP.increaseorder!(r, b1_pkg, b2_pkg, args...)
    end

    discrepancy = Vector{T}(undef, N_tests)
    for n = 1:N_tests
        p = Xs[n]
        discrepancy[n] = abs( f(p) - PowerSeriesIVP.evaltaylor(r.c[begin], p, a) )
    end

    return maximum(discrepancy)
end



############## test RQGeodesicθ

function evalθ(a::T, b::T, x::Vector{T}, u::Vector{T}) where T

    N = length(x)
    @assert length(u) == N

    # stage 1
    Δ = collect( x[i] - x[j] for i in eachindex(x), j in eachindex(x) )
    Δ_sq = Δ .^2

    W4 = collect( sum( Δ_sq[i,k] for k in axes(Δ_sq,2) ) for i in axes(Δ_sq, 1) )
    W8 = collect( sum( Δ[i,k] for k in axes(Δ,2) ) for i in axes(Δ, 1) )

    W4_flip = collect( sum( Δ_sq[i,k] for i in axes(Δ_sq, 1) ) for k in axes(Δ_sq,2) )
    @assert norm(W4-W4_flip) < 1e-15 # should be practically zero.

    # RQGeodesicθ uses W8_flip instead of W8.
    W8_flip = collect( sum( Δ[i,k] for i in axes(Δ, 1) ) for k in axes(Δ,2) )
    @assert norm(W8+W8_flip) < 1e-15 # since W8 is anti-symmetric.

    # stage 1
    a_sq = a^2
    b_sq = b^2
    A = collect( a_sq + W4[i] for i in eachindex(W4) )
    B = collect( b_sq + W4[i] for i in eachindex(W4) )
    R = collect( (a_sq - b_sq)/A[i] for i in eachindex(A) )

    # stage 2
    η = collect( u[i]/B[i] for i in eachindex(B) )
    C = collect( η[i]^2 for i in eachindex(η) )
    #W9 = collect( sum(Δ[i,j]*(u[j]/B[j])^2 for j in eachindex(B) ) for i in axes(Δ,1) )
    #@assert norm(W9+W9_flip) < 1e-15 # since W8 is anti-symmetric.
    W9 = collect( sum(Δ[j,i]*(u[j]/B[j])^2 for j in axes(Δ,1) ) for i in eachindex(B) )

    # stage 3
    W3 = collect( u[i]*sum(Δ[k,i] for k in axes(Δ,1)) for i in eachindex(u) )
    W5 = collect( sum( Δ[j,i]*u[j] for j in axes(Δ,1) ) for i in eachindex(u) )
    W1 = 2 .* W5 .- W3
    
    # stage 4
    W6 = η .* W1
    W7 = B .* W9
    W2 = W6 + W7

    θ_factors = R .* W2
    
    # assemble θ without using stages.
    θ = Vector{T}(undef, N)
    for i in eachindex(θ)

        term1 = -2*sum( Δ[i,j]*u[i]*u[j]/B[i] for j in axes(Δ,2) )
        
        term2 = -sum( B[i]*Δ[i,j]*(u[j]/B[j])^2 for j in axes(Δ,2) )

        term3 = sum( Δ[i,k] for k in axes(Δ,2) )*(u[i]^2)/B[i]

        θ[i] = R[i]*(term1 + term2 + term3)
    end
    @assert norm(θ_factors-θ) < 1e-14

    return θ,
        vec(Δ),  vec(Δ_sq), W4, W8_flip, A, B, R,
        η, C, W9,
        W3, W5, W1,
        W6, W7, W2
    #return θ
end

function evalvecfunc(xs, t::T)::Vector{T} where T
    return collect( xs[i](t) for i in eachindex(xs) )
end

function evalθ(a::T, b::T, xs, us, t::T) where T
    
    return evalθ(a, b, evalvecfunc(xs, t), evalvecfunc(us, t))
end

function evalvectaylor(c::Vector{Vector{T}}, t::T, t0::T)::Vector{T} where T
    return collect( PowerSeriesIVP.evaltaylorwithguard(c[i], t, t0) for i in eachindex(c) )
end
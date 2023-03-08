using Test
using LinearAlgebra

import PowerSeriesIVP

include("../examples/helpers/constraint_helpers.jl")

import Random
Random.seed!(25)

@test true # trivial test for now.

# test the following power series from compositions: reciprocal, multiplication, square, affine transform.

# test horner's algorithm for evaluating polynomials.
@testset "Horner tests" begin

    L = 13 # order of Taylor polynomial.
    N_tests = 1000
    zero_tol = 1e-12

    for _ in N_tests    

        c = randn(L)

        a = randn()
        x = a + 0.1

        sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
        sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
        @test abs(sol-sol_horner) < zero_tol

        #@btime sol_horner = PowerSeriesIVP.evaltaylor(c, x, a)
        #@btime sol = PowerSeriesIVP.evaltaylordirect(c, x, a)
    end

end


# test intersection.
@testset "Intersection tests: verify no intersections in simulated interval" begin

        
    N_tests = 3000
    max_D = 10
    max_N_hyperplanes = 20
    T = Float64

    # solver config.
    L_min = 4
    L_max = 21

    ITP_config = PowerSeriesIVP.ITPConfig(
        Float64;
        f_tol = 1e-8,
        x_tol = 1e-15,
        k1 = 0.1,
        k2 = 0.98*(1+MathConstants.golden), # see ITP paper, equation 24.
        n0 = 0,
    )

    strategy = PowerSeriesIVP.ContinuitySecondDerivative(
        PowerSeriesIVP.GuentherWolfStep(),
    )

    step_config = PowerSeriesIVP.StepConfig(
        Float64,
        strategy;
        #ϵ = 1e-13,# increase this to improve chance that the piece-wise solution is continuous at boundaries.
        ϵ = 1e-6,
        h_max = Inf,
        reduction_factor = 1,
        discount_factor = 0.9,        
    )

    adaptive_order_config = PowerSeriesIVP.AdaptOrderConfig(
        step_config;
        L_min = L_min,
        L_max = L_max, # increase this for maximum higher-order power series.
        order_increase_factor = 1.35,
        max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    )
    #= L_fixed = 4
    fixed_order_config = PowerSeriesIVP.FixedOrderConfig(
        step_config;
        L = L_fixed,
        max_pieces = 100000, # maximum number of pieces in the piece-wise solution.
    ) =#

    # TODO: test random starting point in the tuture.
    t_start = zero(T)

    config = adaptive_order_config

    # tests.
    for _ in N_tests

        D = rand(2:max_D)
        N_hyperplanes = rand(1:max_N_hyperplanes)

        a = abs(randn()*100)
        b = abs(randn()*100)

        t_fin = abs(randn()*100)

        # generate constraints.
        #as, bs = generateHyperplaneConstraints(T, N_hyperplanes, D)
        as, bs = generatecvxpolyhedron(T, N_hyperplanes, D; interior_pt = zeros(T, D))
        lbs, ubs = generateBoundConstraints(T, D)

        bound = PowerSeriesIVP.BoundConstraints(lbs, 1:D, ubs, 1:D)
        hyperplane = PowerSeriesIVP.HyperplaneConstraints(as, bs)
        constraints_info = PowerSeriesIVP.ConstraintsContainer(
            hyperplane,
            bound;
            complex_zero_tol = 1e-8,
            L_min = 4,
            L_max = 10,
            max_divisions = 0,
            solver_config = PowerSeriesIVP.ITPConfig(T),
        )
        constraints = constraints_info.constraints

        # generate initial conditions.

        x0 = randn(D) # starting position.
        feasible_flags = falses(PowerSeriesIVP.getNconstraints(constraints))
        
        while !(all(feasible_flags))
            # if any of the constraints are violated, generate another starting position.
            PowerSeriesIVP.checkconstraints!(feasible_flags, x0, constraints)
            x0 = randn(D)
        end
        u0 = randn(D) .* 100 # starting velocity.

        # transport 2 vector fields.
        N_parallel_vector_fields = 1

        ## set to same as u.
        v0_set = collect( rand(D) for _ = 1:N_parallel_vector_fields)

        # solve IVP
        metric_params = PowerSeriesIVP.RQ22Metric(a,b)
        prob_params = PowerSeriesIVP.GeodesicIVPStatement(metric_params, x0, u0, v0_set)
        sol, exit_flag = PowerSeriesIVP.solveIVP(
            PowerSeriesIVP.getsoltype(prob_params),
            prob_params,
            PowerSeriesIVP.EnableParallelTransport(),
            t_start,
            t_fin,
            config,
            constraints_info,
        )


        # forward seek along sol trajectory.
        constraints = constraints_info.constraints
        end_time = PowerSeriesIVP.getendtime(sol)
        t2, delta_t = PowerSeriesIVP.forwardseekintersection(
            sol,
            constraints,
            t_start,
            end_time,
            PowerSeriesIVP.DisableParallelTransport();
            N_samples = 10000,
        )
        #@show t2, end_time
        #@show x0, as, bs, lbs, ubs

        @test abs(t2-end_time) < eps(T)*2
    end

end
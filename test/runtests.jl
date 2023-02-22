using Test

import PowerSeriesIVP

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
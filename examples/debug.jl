
#### test roots.

Random.seed!(25)

(all_roots0, c) = (ComplexF64[5.683212038777695 + 3.580435972611175im, 5.683212038777695 - 3.580435972611175im, -5.0882430607223075 - 5.417383511233811im, -5.0882430607223075 + 5.417383511233811im, 6.95212520727513e-310 + 6.9521252072783e-310im, 6.95212573078155e-310 + 6.9521245363743e-310im, 6.95212453637825e-310 + 6.9521245363822e-310im, 6.95212453638616e-310 + 6.95212520730675e-310im, 6.95212520711703e-310 + 6.95212453639406e-310im, 6.9521252073099e-310 + 6.9521252073099e-310im, 6.9521245363901e-310 + 6.95212520731307e-310im, 6.9521252073036e-310 + 5.0e-324im, 6.9521257309452e-310 + 6.95212573094756e-310im, 6.95212573094993e-310 + 6.95212520731624e-310im, 6.9521252073194e-310 + 6.95212520732256e-310im, 6.9521257309523e-310 + 6.95212573095705e-310im, 6.9521257309618e-310 + 6.9521252073257e-310im, 6.9521252073289e-310 + 6.95212520733205e-310im, 6.9521252073352e-310 + 2.0e-323im, 2.0e-323 + 1.0e-323im, 5.0e-324 + 0.0im], [2492.26313074836, -168.71452270073695, -15.313574491243756, -1.1899379561107761, 1.0])


all_roots = Vector{Complex{Float64}}(undef, 4)

PowerSeriesIVP.solvequarticequation!(all_roots, c)

function evalpolynomial(c::Vector{T}, t)::Complex{T} where T

    out = c[begin]
    for n = 1:length(c)-1
        out += c[begin+n]*t^n
    end

    return out
end

f = tt->evalpolynomial(c, tt)

f(1.0)

f.(all_roots)

# c_test = [0.1; 0.3; -2; -0.1; 1.0]
# z1, z2, z3, z4 = PowerSeriesIVP.solvequarticequation(c_test)
# @show z1, z2, z3, z4
# @show PowerSeriesIVP.isapproxreal(z1), PowerSeriesIVP.isapproxreal(z2), PowerSeriesIVP.isapproxreal(z3), PowerSeriesIVP.isapproxreal(z4)
# println()

# f = tt->evalpolynomial(c_test, tt)
# @show f(z1), f(z2), f(z3), f(z4)
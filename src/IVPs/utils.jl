
function extractcoefficients(sol::PiecewiseTaylorPolynomial{T,RQGeodesicPiece{T}}) where T

    c_x = collect( sol.coefficients[k].x for k in eachindex(sol.coefficients) )
    c_u = collect( sol.coefficients[k].u for k in eachindex(sol.coefficients) )

    return c_x, c_u
end
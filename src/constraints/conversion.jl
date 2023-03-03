######## prepare the intersection polynomials from the given Taylor polynomial x and the affine cosntraints as, bs.

function updateintersectionpolynomials!(
    cs::Vector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    constraints::SingleTypeConstraints,
    ) where T
    
    N_constraints = getNconstraints(constraints)
    resize!(cs, N_constraints)

    return updateintersectionpolynomials!(cs, x, constraints, 1, N_constraints)
end

# mutates cs.
# the intersection polynomial operate on the variable t-t0, t0 is the expansion point, a constant.
# no error checking on the fields of constraints.
function updateintersectionpolynomials!(
    cs::Vector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    constraints::AllAffineConstraints{T},
    ) where T

    N_affine_constraints = getNconstraints(constraints.affine)
    N_bound_constraints = getNconstraints(constraints.bound)
    N_constraints = N_bound_constraints + N_affine_constraints
    resize!(cs, N_constraints)

    #cs_affine = view(cs, 1:N_affine_constraints)
    updateintersectionpolynomials!(cs_affine, x, constraints.affine, 1, N_affine_constraints)
    
    #cs_bound = view(cs, 1:N_bound_constraints)
    updateintersectionpolynomials!(cs_bound, x, constraints.bound, N_affine_constraints+1, N_constraints)

    return nothing
end

function updateintersectionpolynomials!(
    cs::Vector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    constraints::AffineConstraints{T},
    st_ind::Integer,
    fin_ind::Integer,
    ) where T

    as = constraints.normals # [constraints][variable]
    bs = constraints.offsets # [constraints]

    L_p1 = length(x[begin])

    for m = st_ind:fin_ind
    
    # resize!(cs, length(bs))
    # for m in eachindex(cs)
        
        #cs[m] = zeros(T, L_p1)
        resize!(cs[m], L_p1)
        fill!(cs[m], zero(T))

        a = as[m]
        b = bs[m]

        # this assumes all variable dims partake the m-th constraint, so we'd need a[d] set to 0 if x[d] does not partake in the m-th constraint.
        for d in eachindex(x)

            for l in eachindex(x[d])
                cs[m][l] += x[d][l]*a[d]
            end
        end
        cs[m][begin] -= b
    end

    return nothing
end

function updateintersectionpolynomials!(
    cs::AbstractVector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    constraints::BoundConstraints{T},
    st_ind::Integer,
    fin_ind::Integer,
    ) where T

    # order + 1.
    L_p1 = length(x[begin])

    # [constraints]
    lbs = constraints.lbs 
    lb_dims = constraints.lb_dims

    ubs = constraints.ubs 
    ub_dims = constraints.ub_dims
    
    N_constraints_lb = length(lbs)
    N_constraints_ub = length(ubs)

    # lb <= x becomes -x +lb <= 0.
    for m in Iterators.first(st_ind:fin_ind, N_constraints_lb)
        resize!(cs[m], L_p1)
        
        d = lb_dims[m]
        for l in eachindex(x[d])
            cs[m][l] = -x[d][l]
        end
        cs[m][begin] += lbs[m]
    end

    # x <= ub becomes x -ub <= 0.
    for m in Iterators.drop(st_ind:fin_ind, N_constraints_ub)
        resize!(cs[m], L_p1)

        #for d in eachindex(x) # this assumes all variable dims partake the m-th constraint, so we'd need a[d] set to 0 if x[d] does not partake in the m-th constraint.
        d = ub_dims[m]
        for l in eachindex(x[d])
            cs[m][l] = x[d][l]
        end
        cs[m][begin] -= ubs[m]
    end

    return nothing
end

# c[i] is the coefficient of the i-th degree monomial. i ∈ {0,1,2,...,L}.
# modifies c such that the last entry, i.e. the coefficient for the L-th degree, is 1.
function standardizecoefficients!(c::Vector{T}) where T

    Z = c[end]
    for i in eachindex(c)
        c[i] = c[i]/Z
    end

    return nothing
end
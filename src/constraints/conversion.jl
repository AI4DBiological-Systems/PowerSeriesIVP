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
    constraints::AffineConstraints,
    ) where T

    N_hyperplane_constraints = getNconstraints(constraints.affine)
    N_bound_constraints = getNconstraints(constraints.bound)
    N_constraints = N_bound_constraints + N_hyperplane_constraints
    resize!(cs, N_constraints)

    #cs_affine = view(cs, 1:N_hyperplane_constraints)
    updateintersectionpolynomials!(cs, x, constraints.affine, 1, N_hyperplane_constraints)
    
    #cs_bound = view(cs, 1:N_bound_constraints)
    updateintersectionpolynomials!(cs, x, constraints.bound, N_hyperplane_constraints+1, N_constraints)

    return nothing
end

function updateintersectionpolynomials!(
    cs::Vector{Vector{T}}, # [constraints][order]
    x::Vector{Vector{T}}, # [variable][order]
    constraints::HyperplaneConstraints{T},
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
    constraints::BoundConstraints,
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
    
    N_lb = length(lbs)
    #@show st_ind:fin_ind

    # lb <= x becomes -x +lb <= 0.
    k = 0 # k is an index for bound constraints.
    for m in Iterators.first(st_ind:fin_ind, N_lb) # m is an index for all constraints.
        resize!(cs[m], L_p1)
        
        k += 1
        d = lb_dims[k]
        for l in eachindex(x[d])
            cs[m][l] = -x[d][l]
        end
        cs[m][begin] += lbs[k]
    end

    # x <= ub becomes x -ub <= 0.
    k = 0 # k is an index for bound constraints.
    for m in Iterators.drop(st_ind:fin_ind, N_lb) # m is an index for all constraints.
        resize!(cs[m], L_p1)

        k += 1
        d = ub_dims[k]
        for l in eachindex(x[d])
            cs[m][l] = x[d][l]
        end
        cs[m][begin] -= ubs[k]
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

### convert box bounds on all variables to hyperplanes.
# this assumes all variables are box bounded by lower bounds lbs and upper bounds ubs.
# In addition, there are ordering constrains on the variables, e.g.,
#   x[i] <= x[j] 
function boxtohyperplanes(
    lbs::Vector{T},
    ubs::Vector{T}
    ) where T

    N_vars = length(lbs)
    @assert N_vars == length(ubs)
    
    N_constraints = 2*N_vars
    as = Vector{Vector{T}}(undef, N_constraints)
    bs = Vector{T}(undef, N_constraints)
    
    m = 0
    for d in eachindex(lbs)
        m += 1
        as[m] = zeros(T, N_vars)
        as[m][d] = -one(T)
        bs[m] = lbs[d]
    end

    for d in eachindex(ubs)
        m += 1
        as[m] = zeros(T, N_vars)
        as[m][d] = one(T)
        bs[m] = -ubs[d]
    end

    @assert m == N_constraints

    return as, bs
end

    
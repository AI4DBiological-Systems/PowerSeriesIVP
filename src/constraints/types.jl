
############# root finding, for constraint intersection detection.

struct PositiveRealTrait end


struct RootsBuffer{T<:AbstractFloat}
    intersection_coefficients::Vector{Vector{T}} # [constraints][order]
    all_roots::Vector{Complex{T}} # [order]
    smallest_positive_roots::Vector{T} # [constraints], real roots.
    zero_tol::T # for deciding wheather a complex number variable is a real number.
end

function RootsBuffer(zero_tol::T, order::Integer, N_constraints::Integer)::RootsBuffer{T} where T
    return RootsBuffer(
        collect( zeros(T, order) for _ = 1:N_constraints ),
        Vector{Complex{T}}(undef, order),
        Vector{T}(undef, N_constraints),
        zero_tol,
    )
end

struct RootsUpperBoundBuffer{T}
    # [constraints][order]
    cs::Vector{Vector{T}}
    cs_right::Vector{Vector{T}}
    
    ubs::Vector{T} # buffer, [constraints].

    #cs_left::Vector{Vector{T}}
    #cs_center::Vector{Vector{T}}
end

function RootsUpperBoundBuffer(::Type{T}, N_constraints::Integer, order::Integer)::RootsUpperBoundBuffer{T} where T
    
    return RootsUpperBoundBuffer(
        collect( zeros(T, order) for _ = 1:N_constraints ),
        collect( zeros(T, order) for _ = 1:N_constraints ),
        ones(T, N_constraints) .* NaN,
    )
end


#### numerical solver types (put into separate library later.)
struct ITPConfig{T}
    f_tol::T
    x_tol::T
    k1::T
    k2::T
    n0::Int
end

function ITPConfig(
    ::Type{T};
    f_tol::T = convert(T, 1e-8),
    x_tol::T = convert(T, 1e-15),
    k1::T = convert(T, 0.1),
    k2::T = convert(T, 0.98*(1+MathConstants.golden)), # see equation 24.
    n0::Int = convert(Int, 0),
    )::ITPConfig{T} where T

    @assert k1 > zero(T)
    @assert one(T) <= k2 < one(T) + MathConstants.golden
    @assert n0 >= 0
    @assert x_tol > 0
    @assert f_tol > 0

    return ITPConfig(f_tol, x_tol, k1, k2, n0)
end

#### constraints.
# see conversion.jl for converting the constraint intersection problem to polynomial root problems.

abstract type ConstraintType end

abstract type SingleTypeConstraints <: ConstraintType end

struct HyperplaneConstraints{T} <: SingleTypeConstraints

    # # affine constraints that aren't bound constraints.
    normals::Vector{Vector{T}}
    offsets::Vector{T}
end

function getNconstraints(C::HyperplaneConstraints)::Int
    return length(C.offsets)
end

function HyperplaneConstraints(
    ordering_constraints::Vector{Tuple{Int,Int}}, # if the m-th entry is (i,j), then it means x[i] <= x[j] + b[m]. for the m-th cosntraint.
    bs::Vector{T}, # the m-th entry is bs[m].
    D::Integer,
    )::HyperplaneConstraints{T} where T

    M = length(ordering_constraints)
    @assert length(bs) == M

    as = Vector{Vector{T}}(undef, M)
    for m in eachindex(as)
        as[m] = zeros(T, D)
        
        i, j = ordering_constraints[m]
        as[m][i] = one(T)
        as[m][j] = -one(T)
    end

    return HyperplaneConstraints(as, bs)
end

# although bound constraints can be expressed as affine constraints, 
# we use a new data type so its intersection problem conversion to
#   polynomial root problem is efficient.
# typeof(1:3) <: AbstractVector, so is typeof([1; 3; 44]) <: AbstractVector.
struct BoundConstraints{T, RT<:AbstractArray} <: SingleTypeConstraints
    # lower bound: each entry is a lower bound constraint, i.e. lb <= x[d].
    lbs::Vector{T} # contain the value lb.
    lb_dims::RT # contain the dimension d.

    # upper bound: each entry is a lower bound constraint, i.e. x[d] <= ub.
    ubs::Vector{T}
    ub_dims::RT
end

function getNconstraints(C::BoundConstraints)::Int
    return length(C.lbs) + length(C.ubs)
end

##
abstract type MultiTypeConstraints <: ConstraintType end

struct AffineConstraints{T,RT} <: MultiTypeConstraints
    hyperplane::HyperplaneConstraints{T} # general hyperplane constraints.
    bound::BoundConstraints{T,RT} # bound constraints.
end

function getNconstraints(C::AffineConstraints)::Int
    return getNconstraints(C.hyperplane) + getNconstraints(C.bound)
end

## container, for interfacing with the engine.

abstract type ConstraintsTrait end

struct NoConstraints <: ConstraintsTrait end # this means there are no constraints for all variables.


struct ConstraintsContainer{T,CT<:ConstraintType} <: ConstraintsTrait
    #
    constraints::CT

    # buffers
    explicit_roots_buffer::RootsBuffer{T}
    upperbound_buffer::RootsUpperBoundBuffer{T}
    bino_mat::Matrix{Int}
    max_divisions::Int

    # configs
    solver_config::ITPConfig{T}
end

# this assumes all variables are box bounded by lower bounds lbs and upper bounds ubs.
# In addition, there are ordering constrains on the variables, e.g.,
#   x[i] <= x[j] 
function ConstraintsContainer(
    lbs::Vector{T},
    ubs::Vector{T},
    ordering_constraints::Vector{Tuple{Int,Int}}, # if the m-th entry is (i,j), then it means x[i] <= x[j] + b[m]. for the m-th cosntraint.
    bs::Vector{T}; # the m-th entry is bs[m].
    complex_zero_tol::Real = 1e-8,
    L_min::Integer = 4,
    L_max::Integer = 10,
    max_divisions = 0,
    solver_config = ITPConfig(T),
    )::ConstraintsContainer{T, AffineConstraints{T,UnitRange{Int}}} where T

    D = length(lbs)
    @assert length(ubs) == D

    bound = BoundConstraints(lbs, 1:D, ubs, 1:D)
    hyperplane = HyperplaneConstraints(ordering_constraints, bs, D)
    
    return ConstraintsContainer(
        hyperplane,
        bound;
        complex_zero_tol = complex_zero_tol,
        L_min = L_min,
        L_max = L_max,
        max_divisions = max_divisions,
        solver_config = solver_config,
    )
end

function ConstraintsContainer(
    hyperplane::HyperplaneConstraints{T},
    bound::BoundConstraints{T,RT};
    complex_zero_tol::Real = 1e-8,
    L_min::Integer = 4,
    L_max::Integer = 10,
    max_divisions = 0,
    solver_config = ITPConfig(T),
    )::ConstraintsContainer{T, AffineConstraints{T,RT}} where {T,RT}

    constraints = AffineConstraints(hyperplane, bound)
    N_constraints = getNconstraints(constraints)

    return ConstraintsContainer(
        constraints,
        RootsBuffer( # quartic solver.
            complex_zero_tol,
            L_min,
            N_constraints,
        ),
        RootsUpperBoundBuffer(T, N_constraints, L_max),
        setupbinomialcoefficients(L_max),
        max_divisions,
        solver_config,
    )
end

# this assumes all variables are box bounded by lower bounds lbs and upper bounds ubs.
# In addition, there are ordering constrains on the variables, e.g.,
#   x[i] <= x[j] 
function ConstraintsContainer(
    lbs::Vector{T},
    ubs::Vector{T};
    complex_zero_tol::Real = 1e-8,
    L_min::Integer = 4,
    L_max::Integer = 10,
    max_divisions = 0,
    solver_config = ITPConfig(T),
    )::ConstraintsContainer{T, BoundConstraints{T,UnitRange{Int}}} where T

    D = length(lbs)
    @assert length(ubs) == D

    constraints = BoundConstraints(lbs, 1:D, ubs, 1:D)
    N_constraints = getNconstraints(constraints)

    return ConstraintsContainer(
        constraints,
        RootsBuffer( # quartic solver.
            complex_zero_tol,
            L_min,
            N_constraints,
        ),
        RootsUpperBoundBuffer(T, N_constraints, L_max),
        setupbinomialcoefficients(L_max),
        max_divisions,
        solver_config,
    )
end

function ConstraintsContainer(
    constraints::BoundConstraints{T,RT};
    complex_zero_tol::Real = 1e-8,
    L_min::Integer = 4,
    L_max::Integer = 10,
    max_divisions = 0,
    solver_config = ITPConfig(T),
    )::ConstraintsContainer{T, BoundConstraints{T,RT}} where {T,RT}

    N_constraints = getNconstraints(constraints)

    return ConstraintsContainer(
        constraints,
        RootsBuffer( # quartic solver.
            complex_zero_tol,
            L_min,
            N_constraints,
        ),
        RootsUpperBoundBuffer(T, N_constraints, L_max),
        setupbinomialcoefficients(L_max),
        max_divisions,
        solver_config,
    )
end

function ConstraintsContainer(
    as::Vector{Vector{T}},
    bs::Vector{T};
    complex_zero_tol::Real = 1e-8,
    L_min::Integer = 4,
    L_max::Integer = 10,
    max_divisions = 0,
    solver_config = ITPConfig(T),
    )::ConstraintsContainer{T, HyperplaneConstraints{T}} where T
    
    @assert length(as) == length(bs)
    
    constraints = HyperplaneConstraints(as, bs)

    N_constraints = getNconstraints(constraints)

    return ConstraintsContainer(
        constraints,
        RootsBuffer( # quartic solver.
            complex_zero_tol,
            L_min,
            N_constraints,
        ),
        RootsUpperBoundBuffer(T, N_constraints, L_max),
        setupbinomialcoefficients(L_max),
        max_divisions,
        solver_config,
    )
end
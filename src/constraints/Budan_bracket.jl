

### Budan's upper bound on the number of real roots of a polynomial equation on an interval.




# # no error checking.
# function copybuffer!(x::Vector{Vector{T}}, src::Vector{Vector{T}}) where T
    
#     resize!(x, length(src))

#     for m in eachindex(x)
#         resize!(x[m], length(src[m]))
#         for l in eachindex(x)
#             x[m][l] = src[m][l]
#         end
#     end

#     return nothing
# end


################### Budan's theorem for smallest positive real root isolation.

function setupbinomialcoefficients(order::Integer)::Matrix{Int}
    @assert order > 0
    
    A = Matrix{Int}(undef, order, order)
    for m in axes(A, 2)
        for n in axes(A, 1)
            A[n,m] = binomial(n,m)
        end
    end
    
    return A
end

function updateshiftedpolynomial!(
    out::Vector{T},
    c::Vector{T},
    b::T,
    a::T,
    binomial_mat::Matrix{Int}, # output of setupbinomialcoefficients()
    ) where T

    resize!(out, length(c))
    L = length(c) - 1

    out[begin] = evaltaylor(c, b, a)

    x = b-a
    for m = 1:L
        # TODO: use Horner's method to speed up the evaluation here.
        out[m+1] = sum( c[n+1] * binomial_mat[n,m] * x^(n-m) for n = m:L )
    end
    return nothing
end

function updateshiftedpolynomial!(
    cs_shifted::Vector{Vector{T}},
    cs::Vector{Vector{T}},
    b::T,
    a::T,
    binomial_mat::Matrix{Int}, # output of setupbinomialcoefficients()
    ) where T

    @assert length(cs) == length(cs_shifted)

    for m in eachindex(cs)
        updateshiftedpolynomial!(cs_shifted[m], cs[m], b, a, binomial_mat)
    end

    return nothing
end

# mutates A.
function upperboundintersections!(
    A::RootsUpperBoundBuffer{T},
    sol::PiecewiseTaylorPolynomial{T},
    constraints,
    bino_mat::Matrix{Int},
    ) where T
    
    N_pieces = length(sol.coefficients)
    @assert length(sol.steps) == N_pieces

    ubs = Vector{T}(undef, N_pieces)
    inds = Vector{Int}(undef, N_pieces)

    for n in eachindex(sol.coefficients)

        upperboundNroots!(
            A,
            sol.coefficients[n].x,
            sol.steps[n],
            constraints,
            bino_mat,
        )
        ubs[n], inds[n] = findmax(A.ubs)
    end

    return ubs, inds
end


# based on findfirstroot!()
function upperboundNroots!(
    A::RootsUpperBoundBuffer{T},
    x::Vector{Vector{T}},
    h::T,
    constraints::ConstraintType,
    bino_mat::Matrix{Int},
    ) where T

    updateintersectionpolynomials!(
        A.cs,
        x, 
        constraints,
    )

    cs = A.cs
    cs_right = A.cs_right

    # initialize.
    #copybuffer!(A.cs_left, A.cs)
    #cs_left = A.cs_left
    cs_left = cs
    
    x_right = h

    # get shifted polynomial.
    updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

    ubs = A.ubs
    resize!(ubs, length(cs_left))

    for m in eachindex(cs_left)
        ubs[m] = upperboundroots(cs_left[m], cs_right[m]) 
    end
    
    return nothing
end


function upperboundroots(c_left::Vector{T}, c_right::Vector{T})::Int where T
    return countsignflips(c_left) - countsignflips(c_right)
end

function countsignflips(c::Vector{T})::Int where T
    prev_sign = signbit(c[begin]) # true for negative sign.

    N_flips = 0
    for i in Iterators.drop(eachindex(c),1)
        
        current_sign = signbit(c[i])
        
        if abs(c[i]) > eps() # ignore coefficients with very small magnitude.
            if prev_sign != current_sign
                N_flips += 1
                
                prev_sign = current_sign
            end
        end
    end

    return N_flips
end

## new.

function refinestepnumerical!(
    ::NoConstraints,
    h::T,
    args...
    )::Tuple{T,Int} where T <: AbstractFloat

    return h, 0
end

# returns h, and an integer that signify the hyperplane index if intersected. Return 0 for the integer if no intersection (no roots).
function refinestepnumerical!(
    C::ConstraintsContainer{T},
    h::T,
    h_prev::T, # previously known good step size with no intersections.    
    x::Vector{Vector{T}},
    )::Tuple{T,Int} where T <: AbstractFloat
    
    @assert h > 0
    #@assert h_prev > 0

    A = C.upperbound_buffer
    solver_config = C.solver_config

    # get upperbounds.
    upperboundNroots!(
        A,
        x,
        h,
        C.constraints,
        C.bino_mat,
    )
    ubs = A.ubs
    max_ub, max_ind = findmax(ubs)

    x_right = h

    if max_ub == 0
        # guaranteed no roots.
        return x_right, 0

    elseif max_ub == 1

        # could have multiple constraints intersecting on [0, x_right].
        min_t = convert(T, Inf)
        valid_step = true
        constraint_ind = -1
        for m in eachindex(ubs)
            if ubs[m] == 1
                t, success_flag = runITP(
                    A.cs,
                    zero(T),
                    x_right,
                    solver_config,
                )
                
                valid_step = valid_step & success_flag
                if min_t > t
                    min_t = t
                    constraint_ind = m
                end
            end
        end
        h_new = min_t

        if valid_step
            return h_new, constraint_ind
        end
    
    end

    return convert(T, NaN), max_ind
end


# # treplace this with the ANewDsc algorithm.
# function refinestepnumerical!(
#     C::ConstraintsContainer{T},
#     h::T,
#     h_prev::T, # previously known good step size with no intersections.    
#     x::Vector{Vector{T}},
#     )::Tuple{T,Int} where T <: AbstractFloat
    
#     @assert h > 0
#     @assert h_prev > 0

#     # get upperbounds.
#     upperboundNroots!(
#         C.upperbound_buffer,
#         x,
#         h,
#         C.constraints,
#         C.bino_mat,
#     )
    
#     # we haven't returned yet, so we subdivide.
#     # upperboundNroots!() updated A.cs to have the right content, and A.ubs to have the right length.
#     A = C.upperbound_buffer
#     cs = A.cs
#     cs_right = A.cs_right
#     max_divisions = A.max_divisions
    
#     ubs = A.ubs
#     bino_mat = A.bino_mat

#     x_right = h # x_right = h/2 *2 since the if-else is before updateshiftedpolynomial!().
    
#     max_ind = -1 # allocate for scope reasons. initialize to non-sense value.

#     for _ = 1:max_divisions
        
#         # process upper bounds.

#         max_ub, max_ind = findmax(ubs)
#         # we could have multiple maximums, but the max_ind is not important, so longa s it is not 0 when we have a single root.
#         # we'll force a constraintd etection anyways in the usage of PowerSeriesIVP under the RCG optimization setting.
        
#         if max_ub == 0
#             # guaranteed no roots.
#             return x_right, 0

#         elseif max_ub == 1

#             # could have multiple constraints intersecting on [0, x_right].
#             min_t = convert(T, Inf)
#             valid_step = true
#             constraint_ind = -1
#             for m in eachindex(ubs)
#                 if ubs[m] == 1
#                     t, success_flag = runITP(A, zero(T), x_right)
                    
#                     valid_step = valid_step & success_flag
#                     if min_t > t
#                         min_t = t
#                         constraint_ind = m
#                     end
#                 end
#             end
#             h_new = min_t

#             if valid_step
#                 return h_new, constraint_ind
#             end
        
#         else
#             # in conclusive upperbound of roots.

#             if iter == max_divisions || x_right < h_prev
#                 # give up refining after current refined step is smaller than 1/2^max_divisions times h.

#                 # tell the calling function to decreaseo order to revert to the previously known accepted step size and solutio order.
#                 return convert(T, NaN), max_ind
#             else
#                 # condintue the subdivision. Use bisection rule.
#                 x_right = x_right/2
#             end
#         end
        
#         # get shifted polynomial.
#         updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

#         # get upper bounds for each constraint.
#         for m in eachindex(cs_left)
#             ubs[m] = upperboundroots(cs_left[m], cs_right[m]) 
#         end
        
#     end

#     # tell the calling function to decreaseo order to revert to the previously known accepted step size and solutio order.
#     return convert(T, NaN), max_ind
# end
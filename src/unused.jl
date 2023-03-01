
# # in development. not used for now.
# function findfirstroot!(
#     A::BudanIntersectionBuffers{T},
#     h::T,
#     bino_mat::Matrix{Int};
#     min_bracket_len = 1e-4,
#     max_iters = 100000,
#     ) where T

#     # initialize.
#     #allocatebuffer!(A, cs)
#     copybuffer!(A.cs_left, A.cs)

#     cs = A.cs
#     cs_left = A.cs_left
#     cs_right = A.cs_right
#     cs_center = A.cs_center

#     x_left = zero(T)
#     x_right = h

#     # get shifted polynomial.
#     updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

#     max_val, max_ind = findmax( upperboundroots(cs_left[m], cs_right[m]) for m in eachindex(cs_left))
#     if max_val < 0
#         return h, :redo_quartic # error! contracticts Budan's theorem. Redo this solution piece using quartic polynomials.
#     end

#     if max_val == 0
#         # Guaranteed no roots in this interval.
#         return h, :continue_simulation # safe to declare no intersection on current piece's interval.
        
#     elseif max_val == 1
#         # there is a root in this interval.
#         x_intersection = runITP(x_left, x_right)

#         if isinite(x_intersection)
#             return x_intersection, :end_simulation
#         end

#     end
    
#     # more than one root in this interval, or we failed to find the root in the interval.
#     # we shall subdivide the interval further, and attempt to find the root, or get a new
#     # starting time that is less than h for the next solution polynomial for the IVP simulation.
    
#     new_h, instruction = bracketBudan(
#         cs_center,
#         cs_left,
#         cs_right,
#         x_left,
#         x_right,
#         h,
#         cs,
#         bino_mat;
#         min_bracket_len = min_bracket_len,
#         max_iters = max_iters,
#     )
    
#     return new_h, instruction
# end

# next simulation is either use quartic degree, or we end simulation on the current polynomial.
# function bracketBudan(
#     cs_center::Vector{Vector{T}},
#     cs_left::Vector{Vector{T}},
#     cs_right::Vector{Vector{T}},
#     x_left::T,
#     x_right::T,
#     x_ub::T,
#     cs::Vector{Vector{T}},
#     bino_mat::Matrix{Int};
#     min_bracket_len = min_bracket_len,
#     max_iters = max_iters,
#     ) where T

#     # when we return with an indication to use quartic polynomials in the next solution piece:
#     # set the starting time as x_left for the next polynomial piece, because
#     # we hope to use explicit quartic root solve formula to get the intersection in the next solution piece.

#     for _ = 1:max_iters

#         x_center = (x_left+x_right)/2
#         ub_roots, m_ind, region = splitquery!(
#             cs_center,
#             cs_left,
#             cs_right,
#             x_left,
#             x_right,
#             x_center,
#             cs,
#             bino_mat;
#             min_bracket_len = min_bracket_len,
#         )

#         # see if we can attempt to solve for the root.
#         if ub_roots == 1
            
#             x_intersection = runITP()
#             if isinite(x_intersection)
#                 return x_intersection, :end_simulation
#             end

#             # in case we couldn't find a root using ITP, for whatever reason.
#             if region == :right

#                 return x_center, :use_quartic_next
#             else
#                 return x_left, :use_quartic_next
#             end
        

#         elseif ub_roots == 0
#             # this only happens if instructions == :right: see splitquery!()
            
#             # no roots in the [x_left, x_right] interval.
#             if (x_ub - x_right) < min_bracket_len
#                 return x_right, :continue_simulation
#             end

#             # next interval is: [x-right, x_ub].
#             x_left = x_right
#             x_right = x_ub
        
#         # we shouldn't solve for the root, so we subdivide the interval further.
#         elseif region == :tol_reached
#             # the region between [x_left, x_right] does not have Budan's root upperbound of 1
#             # but, the length of this region is smaller than min_bracket_len.
#             return x_right, :continue_simulation

#         elseif region == :right
#             x_left = x_center
#             #cs_left, cs_center = cs_center, cs_left
#         else
#             x_right = x_center
#             #cs_right, cs_center = cs_center, cs_right
#         end

#     end
    
#     return x_left, :use_quartic_next
# end


# # the assumption is that Budan's test for the interval [x_left, x_right] is non-zero.
# function splitquery!(
#     cs_center,
#     cs_left,
#     cs_right,
#     x_left, x_right, x_center, cs, bino_mat;
#     min_bracket_len = 1e-4)
    
#     if abs(x_left-x_right) < min_bracket_len
#         return tol_reached
#     end

#     #x_center = (x_left + x_right)/2

#     # # query left.
    
#     updateshiftedpolynomial!(cs_center, cs, x_center, zero(T), bino_mat)
#     updateshiftedpolynomial!(cs_left, cs, x_left, zero(T), bino_mat)

#     max_val1, max_ind1 = findmax( upperboundroots(cs_left[m], cs_center[m]) for m in eachindex(cs_left))
    
#     if max_val1 == 0
#         # # query right.

#         updateshiftedpolynomial!(cs_right, cs, x_right, zero(T), bino_mat)

#         max_val12, max_ind2 = findmax( upperboundroots(cs_center[m], cs_right[m]) for m in eachindex(cs_center))

#         return max_val12, max_ind2, :right
#     end

#     return max_val1, max_ind1, :left
# end



# function allocatebuffer!(
#     A::BudanIntersectionBuffers{T},
#     cs::Vector{Vector{T}},
#     ) where T

#     N_constraints = length(cs)
#     @assert N_constraints > 0

#     resize!(A.cs_left, N_constraints)
#     resize!(A.cs_right, N_constraints)
#     resize!(A.cs_center, N_constraints)

#     L = length(cs[begin]) - 1
    
#     for m in eachindex(cs)
    
#         L_p1 = length(cs[m])
#         @assert L == L_p1-1
    
#         resize!(A.cs_left[m], L_p1)
#         resize!(A.cs_right[m], L_p1)
#         resize!(A.cs_center[m], L_p1)
#     end

#     return nothing
# end




############# adaption strat, old.


# # take care of constraints in the piece-wise solution part, not this inner routine for an individual Taylor solution.
# # final order is N_analysis_terms + L_test.
# # Some functions at particular exapnsion points have Taylor series with vanishing coefficients. Example: https://www.wolframalpha.com/input?i=taylor+expansion+of+%281%2Bx%5E4%29%2F%282%2Bx%5E4%29&key=33rzp
# # - N_analysis_terms > 1 attempts to mitigate this issue.
# # - Since N_analysis_terms cannot be infinite, use h_max to guard the case when zero estimated error is detected. h_max = one(T) means no guard, and that zero error means the solution polynomial is exactly the solution of the DE.
# function computetaylorsolution!(
#     prob::GeodesicIVPBuffer,
#     pt_trait::PT,
#     h_initial::T;
#     ϵ::T = convert(T, 1e-6),
#     L_test_max::Integer = 10,
#     r_order::T = convert(T, 0.3),
#     h_max::T = convert(T, 1),
#     N_analysis_terms = 2, # make larger than 2 if the solution has many consecutive vanishing Taylor coefficients for certain orders.
#     ) where {PT,T}

#     ### set up
    
#     p = prob
#     error_threshold = ϵ/2

#     err_record = ones(T, L_test_max) # first entry is for 1st-order, so on.
#     fill!(err_record, Inf)

#     # intermediate buffers.
#     N_vars = length(prob.x0)
#     error_across_variables = Vector{T}(undef, N_vars)

#     ### get the highest compute order N_analysis_terms + 1. The +1 is since we want a minimum of a 1-st order (line) solution.
#     getfirstorder!(prob, pt_trait)
    
#     ### get to a high enough order so that we can start computing the error.
#     for _ = 1:N_analysis_terms
#         increaseorder!(prob, pt_trait)
#     end

#     # are we using the default initial h?
#     h = h_initial
#     if !isfinite(h)
#         # use the following heurestic in this if-block for initial h:

#         # use the first (arb. chosen) variable to analyze how far to step.
#         expansion_factor = 100.0 # hard code for now.
#         h = choosestepsize(ϵ, p.x.c[begin]; h_max = h_max)*expansion_factor
#     end

#     # error for the l-th order.
#     #err_record[begin] = computeerror(p.x.c, h, N_analysis_terms)

#     error_val, error_var_ind = computemaxerror!(error_across_variables, p.x.c, h, N_analysis_terms)
#     err_record[begin] = error_val

#     if error_val < error_threshold

#         # error is already tolerable. no need for adaptation. find step and exit.
#         return choosestepsize(ϵ, p.x.c[error_var_ind]; h_max = h_max)
#     end
    

#     ### start adaption strategy.
#     for l = (1+N_analysis_terms):L_test_max # l is L_test
#         # increase order.
#         increaseorder!(prob, pt_trait)

#         #err_record[l] = computeerror(p.x.c, h, N_analysis_terms) # error for 1st-order solution.
#         error_val, error_var_ind = computemaxerror!(error_across_variables, p.x.c, h, N_analysis_terms)
        
#         # eq:sufficient_error_reduction_order, eq:error_threshold_condition
#         shrink_change = err_record[l]/err_record[l-1]

#         if error_val < error_threshold || shrink_change < r_order
        
#             # Stop increasing the order. Use highest computed order coefficient to get step size.
#             return choosestepsize(ϵ, p.x.c[error_var_ind]; h_max = h_max)
#         end

#         err_record[l] = error_val # book keep.
#     end

#     # Reached max order. Use highest computed order coefficient to get step size.
#     return choosestepsize(ϵ, p.x.c[error_var_ind]; h_max = h_max)
# end

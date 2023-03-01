
############# for refining step in the cosntrained IVP case.






# function getroot()
#     order = getorder()
#     if order == 4
#         #
#         t = getmin_positive_real_root()
#         return t, status

#     elseif order == 3
#         #
#     elseif order == 2
#         #
#     elseif order == 1
#         #
#     elseif order != 0
#         runITP()
#     end

#     "error_order is 0!!! Redo this piece with quartic"
# end



# ############## might be legacy, since we're detecting as we build the piece.
# function decreaseorder!(sol_piece::GeodesicPiece{T}, amount::Integer; min_order = 4) where T
    
#     decreaseorderpolynomials!(sol_piece.x, amount; min_order = min_order)
#     decreaseorderpolynomials!(sol_piece.u, amount; min_order = min_order)
#     decreaseorderpolynomials!(sol_piece.x, amount; min_order = min_order)

#     return nothing
# end

# # decrease the order of the polynomials in cs.
# function decreaseorderpolynomials!(
#     cs::Vector{Vector{T}}, # [index][order]
#     amount::Integer;
#     min_order = 4) where T
    
#     if length(cs[begin]) < min_order+1
#         println("Warning, decreaseintersectionpolynomial!() received an polynomial less than the supplied minimal order. Did not decrease order further.")
#         return nothing
#     end

#     L_p1 = length(cs[m])
#     L_new_p1 = max(L_p1-amount, min_order+1)

#     for m in eachindex(cs)
#         resize!(cs[m], L_new_p1)
#     end

#     return nothing
# end

# function decreaseorder!(
#     A::BudanIntersectionBuffers{T},
#     amount::Integer;
#     min_order = 4) where T

#     decreaseorderpolynomials!(A.cs, amount; min_order = min_order)
#     #decreaseorderpolynomials!(A.cs_left, amount; min_order = min_order)
#     decreaseorderpolynomials!(A.cs_right, amount; min_order = min_order)
#     #decreaseorderpolynomials!(A.cs_center, amount; min_order = min_order)

#     return nothing
# end

# function refinestep!(
#     sol_piece::GeodesicPiece{T},
#     A::BudanIntersectionBuffers{T},
#     h::T,
#     bino_mat::Matrix{Int},
#     ) where T

#     x_left = zero(T)
#     x_right = h

#     max_requests = 3 # hardcode for now. 1/2^3 times h is the most we'll try to recover before reducing order.

#     N_reduced_step_requests = 0
#     allow_reduce_step = true
    
#     for _ = 1:max_iters
        
#         postprocessstep!(A, sol_piece, x_right, constraints, bino_mat)

#         if instruction == :reduce_step
#             x_right = (x_right + x_left)/2
#             N_reduced_step_requests = max(N_reduced_step_requests + 1, max_requests)

#             if allow_reduce_step == max_requests -1
#                 allow_reduce_step = false
#             end
        
#         elseif instruction == :reduce_order
#             order = getorder()
#             if order > 4
#                 decreaseorder!(A, amount; min_order = 4)

#                 # reset step size and requests.
#                 error_val, error_var_ind = computemaxerror(sol_piece.x, h_test, config.N_analysis_terms)
#                 h = choosestepsize(
#                     ϵ,
#                     p.x.c[error_var_ind];
#                     h_max = h_max,
#                     step_reduction_factor = step_reduction_factor,
#                 )

#                 N_reduced_step_requests = 0
#                 allow_reduce_step = true
#             else

#                 # explicit solve, then exit.
#             end
        
#         elseif instruction == check_tol

#             if pass_constraints
#                 # break out of loop.
#                 :stop_simulation
#             else
#                 :reduce_order
#             end
        
#         elseif instruction == :stop_simulation
#             # break out of loop.

#         elseif instruction == :continue_simulation
#             # break out of loop.
#         end

#     end

#     # force quartic solution.

#     return x_right, outgoing_instruction
# end

# # ## upper bound the left region.
# x_center = (x_left+x_right)/2
# updateshiftedpolynomial!(cs_center, cs, x_center, zero(T), bino_mat)
# ub_m = upperboundroots(cs_left[m], cs_center[m])

# if ub_m == 0            
#     min_t = min(x_center, min_t)
    
#     return :exit_shorten_h

# else
#     # inconclusive upperbound. decrease order and try upper bounding again.
    
#     return :exit_reduce_order_then_try_again
#     # decrease order
#     new_order = clamp(order-1, 4, 13) # upperbound so we don't blow the stack with the recursion.

#     # recursion.
#     postprocessstep!()
# end
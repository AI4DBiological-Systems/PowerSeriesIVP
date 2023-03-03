module PowerSeriesIVP

using LinearAlgebra


include("common_types.jl")

include("./sequences/taylor.jl")
include("./sequences/operators.jl")

include("./composite_funcs/composition_utils.jl")
include("./composite_funcs/affine.jl")
include("./composite_funcs/product.jl")
include("./composite_funcs/quotient.jl")
include("./composite_funcs/multivariate.jl")
include("./composite_funcs/integral_seq.jl")

include("./composite_funcs/geodesic_eqns/geodesic_types.jl")
include("./composite_funcs/geodesic_eqns/RQ22.jl")

include("./IVPs/types.jl")
include("./IVPs/methods.jl")
include("./IVPs/adaptive_strategy.jl")
include("./IVPs/engine.jl") # move contents and rename this file.

include("./constraints/conversion.jl")
include("./constraints/Budan_bracket.jl")
include("./constraints/quartic_solver.jl")
include("./constraints/ITP.jl")

include("./testing/continuity.jl")

include("utils.jl")


 # front end for refining step in the constrained IVP case.
#include("./constraints/intersection.jl")

end # module PowerSeriesIVP

# # Assumptions
# We assume all constraints are compatible, and we assume feasible starting position and velocity.

# # Nomenclature
# IVP: initial value problem. We only deal with differentiable, stable ordinary differential equations in this package.
# variable index` in comments mean the index that loop over the variables the IVP solves.
# order index` in comments mean the index that loop over a truncated power series sequence. The first index points to the 0-th order coefficient.

# # References
# (Rodriguez-Bermúdez, 2022): Rodriguez-Bermudez, Panters. "Division of Power Series: Recursive and Non-Recursive Formulas." Anais da Academia Brasileira de Ciências 94 (2022). DOI 10.1590/0001-3765202220210897
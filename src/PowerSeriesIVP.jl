module PowerSeriesIVP

using LinearAlgebra


include("common_types.jl")
include("./constraints/types.jl")
include("./composite_funcs/geodesic_eqns/geodesic_types.jl")

# polynomial evaluation.
include("./sequences/taylor.jl")
include("./sequences/operators.jl")

# composition function library. See PSM 2019 paper.
include("./composite_funcs/composition_utils.jl")
include("./composite_funcs/affine.jl")
include("./composite_funcs/product.jl")
include("./composite_funcs/quotient.jl")
include("./composite_funcs/multivariate.jl")
include("./composite_funcs/integral_seq.jl")

# composition functions for the specific Riemannian Levi-Civita metrics for the geodesic equations family of IVPs.
include("./composite_funcs/geodesic_eqns/RQ22.jl")
# other metrics go here.
include("./composite_funcs/geodesic_eqns/post_methods.jl")

# engine for solving IVPs via PSM.
include("./IVPs/methods.jl") # merge this to engine.jl and methods.jl in geodesic folder.
include("./IVPs/engine.jl") # move contents and rename this file.

include("./IVPs/geodesic_eqns/methods.jl")
include("./IVPs/geodesic_eqns/adaptive_strategy.jl")

# constraint intersection detection.
include("./constraints/conversion.jl")
include("./constraints/Budan_bracket.jl")
include("./constraints/quartic_solver.jl")
include("./constraints/ITP.jl")

# misc.
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
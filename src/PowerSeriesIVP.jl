module PowerSeriesIVP

using LinearAlgebra



include("./sequences/taylor.jl")
include("./sequences/operators.jl")

include("./composite_funcs/composition_utils.jl")
include("./composite_funcs/affine.jl")
include("./composite_funcs/product.jl")
include("./composite_funcs/quotient.jl")
include("./composite_funcs/multivariate.jl")
include("./composite_funcs/integral_seq.jl")

include("./constraints/root_find.jl")

include("./IVPs/types.jl")
include("./IVPs/RQ22.jl")
include("./IVPs/geodesic_engine.jl")
include("./IVPs/engine.jl") # move contents and rename this file.
include("./IVPs/utils.jl")

end # module PowerSeriesIVP

# # Nomenclature
# IVP: initial value problem. We only deal with differentiable, stable ordinary differential equations in this package.
# variable index` in comments mean the index that loop over the variables the IVP solves.
# order index` in comments mean the index that loop over a truncated power series sequence. The first index points to the 0-th order coefficient.



# # References
# (Rodriguez-Bermúdez, 2022): Rodriguez-Bermudez, Panters. "Division of Power Series: Recursive and Non-Recursive Formulas." Anais da Academia Brasileira de Ciências 94 (2022). DOI 10.1590/0001-3765202220210897
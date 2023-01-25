module PowerSeriesIVP

using LinearAlgebra

include("./eval_series/taylor.jl")
include("./composition/operators.jl")

end # module PowerSeriesIVP

# # Nomenclature
# IVP: initial value problem. We only deal with differentiable, stable ordinary differential equations in this package.
# variable index` in comments mean the index that loop over the variables the IVP solves.
# order index` in comments mean the index that loop over a truncated power series sequence. The first index points to the 0-th order coefficient.



# # References
# (Rodriguez-Bermúdez, 2022): Rodriguez-Bermudez, Panters. "Division of Power Series: Recursive and Non-Recursive Formulas." Anais da Academia Brasileira de Ciências 94 (2022). DOI 10.1590/0001-3765202220210897
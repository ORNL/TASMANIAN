module InterfaceJulia

export ClenshawCurtis
export GlobalGrid, Rule1D
export create_lower_set, is_itd_elem
export get_points, get_quadrature_weights, create_XTheta,
    get_points_and_quadrature_weights, lex_merge_2d_arrays

include("TasData.jl")
include("TasOneDimesionalRule.jl")
include("TasUtil.jl")
include("TasGrid.jl")

end

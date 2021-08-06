# A collection of subroutines, and tests on some subroutines, that appear in
# various key functions.
module proto

include("src/InterfaceJulia.jl")
using .InterfaceJulia

# ==============================================================================
# PROTO 1
# ==============================================================================

# # Specialized permutation function.
# function spec_permute(n)
#     p = [1, n+1]
#     for k=1:(n-1)
#         mask = digits(k, base=2)
#         println(mask)
#         pow2_seq = 2 .^ collect(1:length(mask))
#         println(pow2_seq)
#         pk = n .* sum(mask .// pow2_seq) + 1
#         println(pk)
#         if denominator(pk) == 1
#             append!(p, Int(pk))
#         else
#             error("Permutation index needs to be Int!")
#         end
#     end
#     return(p)
# end

# ==============================================================================
# PROTO 2
# ==============================================================================

# # Clenshaw-Curtis grid of dimension 2 and total degree 4.
# cc_rule = ClenshawCurtis()
# ls = create_lower_set(2, x->TasData.is_itd_elem(4, x))
# g = GlobalGrid([cc_rule, cc_rule], ls)
# pg = get_points(g)
# wg = get_quadrature_weights(g)

# # Use the grid to compute
# #
# #   I := ∫ₐ∫ₐ sin(x²+y²) dx dy
# #
# # where a := [-1, 1].
# int_fn = (t) -> sin(t[1]^2 + t[2]^2)
# I = sum(wg .* int_fn.(pg))
# cmp_mat = Matrix(undef, 2, 1)
# cmp_mat[1, 1] = 2.245161593287624
# cmp_mat[2, 1] = I

# ==============================================================================
# PROTO 2
# ==============================================================================

# # Clenshaw-Curtis grid of dimension 3 and total degree 4.
# dimension = 2
# level = 8
# ls = create_lower_set(dimension, x->is_itd_elem(level, x))
# g = GlobalGrid(ClenshawCurtis(), ls)
# pg = get_point_cache(g)
# wg = get_quad_weight_cache(g)

# pg = get_points(g)
# qwg = get_quadrature_weights(g)
# Profile.print(format=:flat, sortedby=:count, mincount=100)

# ==============================================================================
# PROTO 2
# ==============================================================================

# function cartesian(elem_set::Vector{Vector{T}}) where T<:Any
#     # Efficiently computes the Cartesian product of a set of vectors.
#     num_dims = length(elem_set)
#     num_prod_elems = prod(length.(elem_set))
#     cprod_set = Matrix{T}(undef, num_dims, num_prod_elems)
#     rep_elems = num_prod_elems
#     for k=1:num_dims
#         num_elems = length(elem_set[k])
#         rep_elems = rep_elems ÷ num_elems
#         for i=1:num_prod_elems
#             # Get the index of the element in elem_set[k] to grab.
#             j = (((i-1) ÷ rep_elems) % num_elems) + 1
#             cprod_set[k,i] = elem_set[k][j]
#         end
#     end
#     return(cprod_set)
# end

# ==============================================================================
# PROTO 3
# ==============================================================================

# a = [[1, 2] [1, 3] [3, 3]]
# b = [[1, 1] [1, 3] [2, 3] [4, 1]]
# va = [[1.1, 2.3] [2.8, 0.1] [3.1, 1.3]]
# vb = [[1.0, 0.1] [0.2, 0.3] [2.0, 0.3] [0.4, 1.1]]

# vc, c = lex_merge_2d_arrays(va, vb, a, b)

# ==============================================================================
# PROTO 4
# ==============================================================================

# Utility function that generates a global grid of a given dimension and degree
# with isotropic total degree lower set and each dimension uses the nested
# Clenshaw-Curtis interpolation rule.
function create_cc_grid(dimension::Int, degree::Int)
    grid_lower_set = create_lower_set(dimension, x->is_itd_elem(degree, x))
    return GlobalGrid(ClenshawCurtis(), grid_lower_set)
end

# Utility function that computes the integral given by a function applied
# to a set of weights and points (listed in column major order).
function compute_integral(fn::Function,
                          pts::Matrix{Float64},
                          wts::Vector{Float64})
    return sum(wts .* map(fn, eachcol(pts)))
end

# # Create a grid of total degree 4 and dimension 1.
# grid = create_cc_grid(1, 4)
# cc_points = get_points(grid)
# cc_quad_weights = get_quad_weights(grid)
# # Compute ∫ₐ sin(x²) dx dy where a=[-1,1].
# integrand = z -> sin(z[1]^2)
# integral = compute_integral(integrand, cc_points, cc_quad_weights)
# println([integral, 0.6205366034467622])

# Create a grid of dimension 3 and some large degree
grid = create_cc_grid(3, 12)
cc_points = get_points(grid)
cc_quad_weights = get_quad_weights(grid)
display(size(cc_points))

# ==============================================================================
# PROTO 5
# ==============================================================================

# dimension = 2
# level = 4
# ls = create_lower_set(dimension, x->is_itd_elem(level, x))
# g = GlobalGrid(ClenshawCurtis(), ls)
# is = g.index_set
# ps = get_points(g)
# ws = get_quad_weights(g)

end # proto module

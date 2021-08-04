# A collection of subroutines, and tests on some subroutines, that appear in
# various key functions.

using Profile
using InterfaceJulia.TasOneDimensionalRule
using InterfaceJulia.TasData
using InterfaceJulia.TasGrid

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
# cc_rule = ClenshawCurtis()
# dimension = 2
# level = 4
# ls = create_lower_set(dimension, x->TasData.is_itd_elem(level, x))
# g = GlobalGrid([cc_rule, cc_rule], ls)
# qwg = get_quadrature_weights(g)
# pg = get_points(g)

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

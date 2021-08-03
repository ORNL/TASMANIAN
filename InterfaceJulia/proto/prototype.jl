# A collection of subroutines, and tests on some subroutines, that appear in
# various key functions.
module TasProto

include("../src/TasData.jl")
using .TasData
include("../src/TasOneDimesionalRule.jl")
using .TasOneDimensionalRule
include("../src/TasGrid.jl")
using .TasGrid
include("../src/TasUtil.jl")
using .TasUtil


# ==============================================================================
# PROTO 1
# ==============================================================================

# Specialized permutation function.
function spec_permute(n)
    p = [1, n+1]
    for k=1:(n-1)
        mask = digits(k, base=2)
        println(mask)
        pow2_seq = 2 .^ collect(1:length(mask))
        println(pow2_seq)
        pk = n .* sum(mask .// pow2_seq) + 1
        println(pk)
        if denominator(pk) == 1
            append!(p, Int(pk))
        else
            error("Permutation index needs to be Int!")
        end
    end
    return(p)
end

# ==============================================================================
# PROTO 2
# ==============================================================================

# Clenshaw-Curtis grid of dimension 2 and total degree 4.
cc_rule = Rule1D(true, l->2^l+1, clenshaw_curtis)
ls = create_lower_set(2, x->TasData.is_itd_elem(7, x))
g = GlobalGrid([cc_rule, cc_rule], ls)
wcg = generate_weight_cache(g)
scg = generate_quad_surplus_cache(g)
pg = get_points(g)
wg = get_quadrature_weights(g)

# Use the grid to compute
#
#   I := ∫ₐ∫ₐ sin(x²+y²) dx dy
#
# where a := [-1, 1].
int_fn = (t) -> sin(t[1]^2 + t[2]^2)
I = sum(wg .* int_fn.(pg))
cmp_mat = Matrix(undef, 2, 1)
cmp_mat[1, 1] = 2.245161593287624
cmp_mat[2, 1] = I

end # module

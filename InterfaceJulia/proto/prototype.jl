# A collection of subroutines, and tests on some subroutines, that appear in
# various key functions.

module TasProto

include("../src/TasData.jl")
using .TasData
include("../src/TasOneDimesionalRule.jl")
using .TasOneDimensionalRule
include("../src/TasGrid.jl")
using .TasGrid

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
g = GlobalGrid([true, true],
               [clenshaw_curtis, clenshaw_curtis],
               [[1, 1] [1, 2]])
wg = generate_weight_cache(g)
sg = generate_surplus_cache(g)
print(sg)

end # module

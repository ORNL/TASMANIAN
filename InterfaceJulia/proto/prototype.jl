# A collection of subroutines that appear in various key functions.

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

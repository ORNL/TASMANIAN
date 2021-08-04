# A collection of helpful utility functions.
module TasUtil

export cartesian_product

# Takes in a list of vectors and computes the full cartesian product of them.
# The product should be in lexicographical order.
function cartesian_product(elem_set::Vector{Vector{T}}) where T<:Any
    # Efficiently computes the Cartesian product of a set of vectors.
    num_dims = length(elem_set)
    num_prod_elems = prod(length.(elem_set))
    cprod_set = Matrix{T}(undef, num_dims, num_prod_elems)
    rep_elems = num_prod_elems
    for k=1:num_dims
        num_elems = length(elem_set[k])
        rep_elems = rep_elems ÷ num_elems
        for i=1:num_prod_elems
            # Get the index of the element in elem_set[k] to grab.
            j = (((i-1) ÷ rep_elems) % num_elems) + 1
            cprod_set[k,i] = elem_set[k][j]
        end
    end
    return(cprod_set)
end


end # module

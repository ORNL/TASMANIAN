# A collection of helpful utility functions and types.

# Union types
SubMatrix{T} = Union{Matrix{T}, SubArray{T, 2}} where T

function cartesian_product(elem_set::Vector{Vector{T}}) where T<:Any
    # Takes in a list of vectors and computes the full cartesian product of them.
    # The product should be in lexicographical order.
    num_dims = length(elem_set)
    num_prod_elems = prod(length.(elem_set))
    cprod_set = Matrix{T}(undef, num_dims, num_prod_elems)
    rep_elems = num_prod_elems
    for k=1:num_dims
        num_elems = length(elem_set[k])
        rep_elems = rep_elems ÷ num_elems
        for i=1:num_prod_elems
            # Get the index of the element in elem_set[k] to grab.
            j = (((i - 1) ÷ rep_elems) % num_elems) + 1
            cprod_set[k,i] = elem_set[k][j]
        end
    end
    return(cprod_set)
end

function lex_merge_matrices(ind1::Matrix{Int}, ind2::Matrix{Int})
    # Takes two multi-index arrays, whose columns are in lexicographical order,
    # and merges their columns together. The output array is also in
    # lexicographical order.
    if size(ind1, 1) != size(ind2, 1)
        error("All index matrices have to have the same number of rows.")
    end
    lex_ind_mat = Matrix{Int}(undef, size(ind1, 1), 0)
    (i1, i2) = (1, 1)
    while i1 <= size(ind1, 2) || i2 <= size(ind2, 2)
        if i1 > size(ind1, 2)
            next_ind_elem = ind2[:,i2]
            i2 += 1
        elseif i2 > size(ind2, 2)
            next_ind_elem = ind1[:,i1]
            i1 += 1
        else
            if  ind1[:,i1] > ind2[:,i2]
                next_ind_elem = ind2[:,i2]
                i2 += 1
            elseif ind1[:,i1] < ind2[:,i2]
                next_ind_elem = ind1[:,i1]
                i1 += 1
            else
                next_ind_elem = ind1[:,i1]
                (i1, i2) = (i1 + 1, i2 + 1)
            end
        end
        lex_ind_mat = [lex_ind_mat next_ind_elem]
    end
    return lex_ind_mat
end

function lex_find_idx(a::SubMatrix{Int}, x::Vector{Int})
    # Find a the column index multi-index in a lexicographically sorted
    # multi-index array, listed in column major order.
    function binary_search(a::SubMatrix{Int}, x::Vector{Int}, lo::Int, hi::Int)
        # Binary search helper function.
        if lo <= hi
            mid = lo + (hi - lo) ÷ 2
            if a[:, mid] == x
                return mid
            elseif a[:, mid] > x
                binary_search(a, x, lo, mid-1)
            else
                binary_search(a, x, mid+1, hi)
            end
        else
            error("Index not found!")
        end
    end
    return binary_search(a, x, 1, size(a, 2))
end




# A collection of helpful utility functions.

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

function lex_merge_2d_arrays(val1::Array{Float64,2}, val2::Array{Float64,2},
                             ind1::Array{Int,2}, ind2::Array{Int,2},
                             op::Function)
    # Takes two value matrices, two index arrays, whose columns are in
    # lexicographical order, and a binary operator. Merges columns in the value
    # and index arrays so that: (i) output index array is also lexicographically
    # sorted; (ii) the output value array has the same ordering as the output
    # index array; values of the same index are resolved using the binary
    # operator.
    if size(ind1, 1) != size(ind2, 1)
        error("All index matrices have to have the same number of rows.")
    end
    if size(val1, 2) != size(ind1, 2) || size(val2, 2) != size(ind2, 2)
        error("Value and index array pairs do not have matching number of " *
              "columns!")
    end
    lex_val_mat = Matrix{Float64}(undef, size(val1, 1), 0)
    lex_ind_mat = Matrix{Int}(undef, size(ind1, 1), 0)
    c1 = size(val1, 2)
    c2 = size(val2, 2)
    i1 = 1;
    i2 = 1;
    while i1 <= c1 || i2 <= c2
        if i1 > c1
            next_val_elem = val2[:,i2]
            next_ind_elem = ind2[:,i2]
            i2 += 1
        elseif i2 > c2
            next_val_elem = val1[:,i1]
            next_ind_elem = ind1[:,i1]
            i1 += 1
        else
            if ind1[:,i1] == ind2[:,i2]
                next_val_elem = op(val1[:,i1], val2[:,i2])
                next_ind_elem = ind1[:,i1]
                i1 += 1
                i2 += 1
            elseif  ind1[:,i1] > ind2[:,i2]
                next_val_elem = val2[:,i2]
                next_ind_elem = ind2[:,i2]
                i2 += 1
            elseif ind1[:,i1] < ind2[:,i2]
                next_val_elem = val1[:,i1]
                next_ind_elem = ind1[:,i1]
                i1 += 1
            else
                error("Cannot compare elements!")
            end
        end
        lex_val_mat = [lex_val_mat next_val_elem]
        lex_ind_mat = [lex_ind_mat next_ind_elem]
    end
    return lex_val_mat, lex_ind_mat
end


function lex_merge_matrices(ind1::Array{Int,2}, ind2::Array{Int,2})
    # Takes two index matrices, whose columns are in lexicographical order.
    # Merges columns in the arrays so that: (i) output index array is also
    # lexicographic order and (ii) duplicates are removed
    if size(ind1, 1) != size(ind2, 1)
        error("All matrices have to have the same number of rows.")
    end
    merged_mat = Matrix{Int}(undef, size(ind1, 1), 0)
    i1 = 1;
    i2 = 1;
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
            else ind1[:,i1] == ind2[:,i2]
                next_ind_elem = ind1[:,i1]
                i1 += 1
                i2 += 1
            end
        end
        merged_mat = [merged_mat next_ind_elem]
    end
    return merged_mat
end

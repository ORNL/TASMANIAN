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

function lex_merge_2d_arrays(mat1::Array{T,2}, mat2::Array{T,2}) where T<:Integer
    # Takes two integer matrices, whose columns are in lexicographical order, and
    # merges their columns into another matrix, whose columns are also in
    # lexicographical order.
    if size(mat1, 1) != size(mat2, 1)
        error("Both matrices have to have the same number of rows.")
    end
    lex_mat = Matrix{T}(undef, size(mat1, 1), 0)
    c1 = size(mat1, 2)
    c2 = size(mat2, 2)
    i1 = 1;
    i2 = 1;
    while i1 <= c1 || i2 <= c2
        if i1 > c1
            next_elem = mat2[:,i2]
            i2 += 1
        elseif i2 > c2
            next_elem = mat1[:,i1]
            i1 += 1
        else
            if mat1[:,i1] == mat2[:,i2]
                next_elem = mat1[:,i1]
                i1 += 1
                i2 += 1
            elseif  mat1[:,i1] > mat2[:,i2]
                next_elem = mat2[:,i2]
                i2 += 1
            elseif mat1[:,i1] < mat2[:,i2]
                next_elem = mat1[:,i1]
                i1 += 1
            else
                error("Cannot compare elements!")
            end
        end
        lex_mat = [lex_mat next_elem]
    end
    return lex_mat
end

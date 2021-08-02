# A collection of functions for generating key data structures.
module TasData
export create_lower_set, is_itd_elem

function create_lower_set(num_dimensions::Int, is_elem::Function)
    #=
    Function that returns a column-sorted matrix whose columns satisfy
    a conditional *is_elem*.
    =#
    if num_dimensions == 0
        return Int[]
    end
    num_entries = 0
    indexes = Int[]
    p = ones(Int, num_dimensions)
    p_idx = num_dimensions
    while(is_elem(p) || p_idx > 1)
        if is_elem(p)
            append!(indexes, p)
            num_entries += 1
            p_idx = num_dimensions
            p[p_idx] += 1
        else
            p[p_idx:end] .= 1
            p_idx -= 1
            p[p_idx] += 1
        end
    end
    return reshape(indexes, num_dimensions, num_entries)
end

function is_itd_elem(m::Int, index::Vector{Int})
    #=
    Test if an index is part of the isotropic total degree lower set of
    degree m.
    =#
    return sum(index) <= m
end

end # module

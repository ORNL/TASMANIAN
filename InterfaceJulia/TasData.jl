module TasData
export MultiIndexSet, StorageSet, create_lower_set

# ==============================================================================
# Key Composite Types
# ==============================================================================

struct MultiIndexSet{T<:Integer}
    # Representation of a collection of multidimensional indexes of dimensions
    # *num_dimensions*. The indexes are assumed to be sorted and the data
    # is stored in column major order.

    # Fields.
    num_dimensions::T # rows
    cache_num_indexes::T # columns
    indexes::Array{T}

    # Constructors.
    MultiIndexSet{T}() where T<:Integer = new{T}(0, 0, [])
    MultiIndexSet() = MultiIndexSet{Int}(0, 0, [])
    function MultiIndexSet(input_indexes::Matrix{T}) where T<:Integer
        if length(size(input_indexes)) != 2
            error("Matrices must be 2D!")
        end
        return new{T}(
            size(input_indexes, 2), size(input_indexes, 1), input_indexes)
    end
    function MultiIndexSet(ndim::T, cnum::T, input_indexes::Vector{T}) where T<:Integer
        return MultiIndexSet(reshape(input_indexes, ndim, cnum))
    end
end

struct StorageSet{T<:Integer}

    # Fields.
    num_outputs::T # rows
    num_values::T # columns
    values::Array{S,2} where S<:Real

    # Constructors.
    StorageSet{T}() where T<:Integer = new{T}(0, 0, [])
    StorageSet() = StorageSet{Int}(0, 0, [])
    function StorageSet(input_indices::Array{T}) where T<:Integer
        if length(size(input_indices)) != 2
            error("Input must be a 2D array!")
        end
        if stride(input_indices, 1) != 1
            error("Elements of the array must be stored in column major order!")
        end
        return new{T}(
            size(input_indices, 2), size(input_indices, 1), input_indices)
    end
end

# ==============================================================================
# Key Functions
# ==============================================================================

function create_lower_set(num_dimensions::Int, is_elem::Function)
    # Function that returns a (sorted) MultiIndexSet whose elements satisfy
    # a conditional *is_elem*.

    num_entries = 0
    indexes = Int[]
    p = ones(Int, num_dimensions)
    p_is_elem = is_elem(p)
    p_idx = num_dimensions
    while(p_is_elem || p_idx > 1)
        if p_is_elem
            append!(indexes, p)
            num_entries += 1
            p_idx = num_dimensions
            p[p_idx] += 1
        else
            p[p_idx:end] .= 1
            p_idx -= 1
            p[p_idx] += 1
        end
        p_is_elem = is_elem(p) # Check for the next multi-index
    end
    return MultiIndexSet(num_dimensions, num_entries, indexes)
end

function is_itd_elem(m::Int, index::Vector{Int})
    # Test if an index is part of the isotropic total degree lower set of
    # degree m
    return sum(index) <= m
end

# Test to create a lower set of indices of total degree 3
foo = create_lower_set(3, x -> is_itd_elem(5, x))


end

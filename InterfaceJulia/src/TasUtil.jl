# A collection of helpful utility functions.
module TasUtil

export cartesian_product

# Takes in a list of vectors and computes the full cartesian product of them.
# The product should be in lexicographical order.
function cartesian_product(ll::Vector{Vector{T}}) where T <: Any
    return vec(collect(Base.product(ll...)))
end

end # module

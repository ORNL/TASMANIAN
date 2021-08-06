# A collection of 1D interpolation functions that return points and
# weights in column major format.

struct Rule1D
    #=
    Container for a 1D interpolation rule.
    =#

    # Indicator for if the rule generates nested points.
    nested::Bool
    # Domain of interpolation [a,b] as a tuple (a,b).
    domain::Tuple{Float64,Float64}
    # Function that returns the number of nodes for a given level, i.e., m(⋅).
    num_nodes::Function
    # Function that returns the points and weights in a pair (points, weights)
    # for a given level. Indices for points and weights should align!
    points_and_weights::Function

    # Constructors
    Rule1D() = new(false, [-1.0, 1.0], l->0, l->([],[]))
    Rule1D(n::Bool, d::Tuple{Float64,Float64}, nn::Function, pw::Function) =
        new(n, d, nn, pw)
end

# Clenshaw-Curtis interpolation rule.
function clenshaw_curtis_points_and_weights(l::Int)
    #=
    Generates (*points*, *weights*) following the Clenshaw-Curtis rule over
    the interval [-1,1] for a given degree l, i.e. m(i). The construction
    is taken from (Sommariva, 2013) and m(l) = 2^l + 1 points are generated.
    =#
    if l == 0
        points = zeros(1, 1)
        weights = zeros(1, 1)
    else
        n = 2 ^ l # equal to m(l) - 1
        n_div_2 = Int(n / 2)
        rg_n_div_2 = float(1:n_div_2)
        v = (pi / n) * float(0:n)
        b = [2.0 * ones(n_div_2 - 1); 1.0]
        c = [1.0; 2.0 * ones(n - 1); 1.0]
        points = cos.(v)
        weights = Vector{Float64}(undef, n + 1)
        for i = 1:(n+1)
            sub_arr =
                b ./ (4.0 * rg_n_div_2 .^ 2.0 .- 1.0) .*
                cos.(2.0 * rg_n_div_2 * v[i])
            weights[i] = c[i] / n * (1.0 - sum(sub_arr))
        end
    end
    #=
    Permute the points so that they are nested as you increase l. That is,
    if CC(l) is the output for degree l, then CC(l+1) = [CC(l) t(l+1)]
    where t(l+1) are the additional outputs in CC(l+1) that are not in CC(l).
    =#
    p = [1, n + 1]
    for k=1:(n-1)
        mask = digits(k, base=2)
        pow2_seq = 2 .^ collect(1:length(mask))
        pk = n .* sum(mask .// pow2_seq) + 1
        if denominator(pk) == 1
            append!(p, Int(pk))
        else
            error("Permutation index needs to be Int!")
        end
    end
    return points[p], weights[p]
end
ClenshawCurtis() = Rule1D(true,
                          (-1.0, 1.0),
                          l -> l == 0 ? 1 : 2^l+1,
                          clenshaw_curtis_points_and_weights)

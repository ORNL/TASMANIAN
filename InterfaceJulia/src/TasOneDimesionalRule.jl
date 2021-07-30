module TasOneDimensionalRule

export clenshaw_curtis

function clenshaw_curtis(n::Int)
    # Generates (*points*, *weights*) following the Clenshaw-Curtis rule over
    # the interval [-1,1] for a given level *n*, i.e. m(i). The construction
    # is taken from (Sommariva, 2013).
    if isodd(n)
        throw(DomainError("n is not even!"))
    end
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
    return points, weights

end

end

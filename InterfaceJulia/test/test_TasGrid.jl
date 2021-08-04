using XUnit
using InterfaceJulia

# Utility function that generates a global grid of a given dimension and degree
# with isotropic total degree lower set and each dimension uses the nested
# Clenshaw-Curtis interpolation rule.
function create_cc_grid(dimension::Int, degree::Int)
    cc_rule = ClenshawCurtis()
    grid_rule_vec = Array{Rule1D}(undef, 0)
    for _=1:dimension
        push!(grid_rule_vec, cc_rule)
    end
    grid_lower_set = create_lower_set(
        dimension, x->is_itd_elem(degree, x))
    return GlobalGrid(grid_rule_vec, grid_lower_set)
end

# Utility function that computes the integral given by a function applied
# to a set of weights and points (listed in column major order).
function compute_integral(fn::Function,
                          pts::Matrix{Float64},
                          wts::Vector{Float64})
    return sum(wts .* map(fn, eachcol(pts)))
end

# Test the accuracy of a grid with an isotropic total degree lower set and
# a Clenshaw-Curtis 1D rule.
@testset "TasGrid" begin
    @testcase "Clenshaw-Curtis 1D Quadrature" begin
        # Create a grid of total degree 4 and dimension 1.
        grid = create_cc_grid(1, 4)
        cc_points = get_points(grid)
        cc_quad_weights = get_quadrature_weights(grid)
        # Compute ∫ₐ sin(x²) dx dy where a=[-1,1].
        integrand = z -> sin(z[1]^2)
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 0.6205366034467622
        # Compute ∫ₐ exp(-x²) dx dy where a=[-1,1].
        integrand = z -> exp(-z[1]^2)
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 1.493648265624854
        # Compute ∫ₐ sqrt(arccos(x^2/)) dx dy where a=[-1,1].
        integrand = z -> sqrt(acos(z[1]^2/2))
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 2.3634629916571772
    end
    @testcase "Clenshaw-Curtis 2D Quadrature" begin
        # Create a grid of total degree 6 and dimension 2.
        grid = create_cc_grid(2, 6)
        cc_points = get_points(grid)
        cc_quad_weights = get_quadrature_weights(grid)
        # Compute ∫ₐ∫ₐ sin(x²+y²) dx dy where a=[-1,1].
        integrand = z -> sin(z[1]^2 + z[2]^2)
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 2.245161593287624
        # Compute ∫ₐ∫ₐ exp(-(x²+y²)) dx dy where a=[-1,1].
        integrand = z -> exp(-z[1]^2 - z[2]^2)
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 2.230985141404135
        # Compute ∫ₐ∫ₐ exp(-x²) * cos(y) dx dy where a=[-1,1] (taken from
        # example_sparse_grids_01.py).
        integrand = z -> exp(-z[1]^2) * cos(z[2])
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 2.513723354063905
    end
    @testcase "Clenshaw-Curtis 3D Quadrature" begin
        # Create a grid of total degree 8 and dimension 3.
        grid = create_cc_grid(3, 8)
        cc_points = get_points(grid)
        cc_quad_weights = get_quadrature_weights(grid)
        # Compute ∫ₐ∫ₐ∫ₐ exp(-(x²+y²+z²)) dx dy dz where a=[-1,1].
        integrand = z -> exp(-z[1]^2-z[2]^2-z[3]^2)
        integral = compute_integral(integrand, cc_points, cc_quad_weights)
        @test integral ≈ 3.332307087093105
    end
end


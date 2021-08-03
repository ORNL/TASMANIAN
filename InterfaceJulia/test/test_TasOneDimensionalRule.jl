using XUnit
using InterfaceJulia.TasOneDimensionalRule

@testset "TasOneDimensionalRule" begin
    @testcase "Clenshaw-Curtis 1D Quadrature" begin
	      points, weights = clenshaw_curtis(3)
        # Compute ∫ₐ arccos(x) dx where a = [-1,1].
        @test sum(acos.(points) .* weights) ≈ pi
        # Compute ∫ₐ cos(πx/2) dx where a = [-1,1].
        @test sum(cos.(pi * points / 2) .* weights) ≈ 4 / pi
    end
end

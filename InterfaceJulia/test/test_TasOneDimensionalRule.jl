using Test
using InterfaceJulia.TasOneDimensionalRule

@testset "Clenshaw-Curtis Quadrature" begin
	  points, weights = clenshaw_curtis(3)
    @test sum(acos.(points) .* weights) ≈ pi
    @test sum(cos.(pi * points / 2) .* weights) ≈ 4 / pi
end

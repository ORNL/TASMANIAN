using Test
using InterfaceJulia.TasOneDimensionalRule

@testset "Clenshaw-Curtis Quadrature" begin
	  points, weights = clenshaw_curtis(10)
    @test sum(acos.(points) .* weights) ≈ pi
end

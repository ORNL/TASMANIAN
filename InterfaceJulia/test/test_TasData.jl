using XUnit
using InterfaceJulia

# Create a lower sets of various degrees.
@testset "TasData" begin
    @testcase "Isotropic Total Degree Set" begin
        @test isempty(create_lower_set(0, x -> is_itd_elem(0, x)))
        @test isempty(create_lower_set(0, x -> is_itd_elem(1, x)))
        @test isempty(create_lower_set(1, x -> is_itd_elem(0, x)))
        @test create_lower_set(1, x -> is_itd_elem(1, x)) == ones(Int, 1, 1)
        @test create_lower_set(1, x -> is_itd_elem(2, x)) ==
            reshape([1, 2], 1, 2)
        @test isempty(create_lower_set(2, x -> is_itd_elem(1, x)))
        @test create_lower_set(2, x -> is_itd_elem(2, x)) ==
            reshape([1; 1], 2, 1)
        @test create_lower_set(2, x -> is_itd_elem(4, x)) ==
            [[1, 1] [1, 2] [1, 3] [2, 1] [2, 2] [3, 1]]
    end
end

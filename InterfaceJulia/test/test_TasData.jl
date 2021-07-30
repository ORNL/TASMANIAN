using Test
using InterfaceJulia.TasData

# TODO Test to create a lower set of indices of total degree 5
@testset "Isotropic Total Degree Set" begin
    foo = create_lower_set(3, x -> TasData.is_itd_elem(5, x))
end

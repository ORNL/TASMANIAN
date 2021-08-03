using XUnit
using InterfaceJulia.TasOneDimensionalRule

@time @testset runner=ParallelTestRunner() "Tasmanian [InterfaceJulia]" begin

    include("test_TasData.jl")
    include("test_TasOneDimensionalRule.jl")
    include("test_TasGrid.jl")

end

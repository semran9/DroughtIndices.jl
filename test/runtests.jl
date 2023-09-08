using DroughtIndices
using Test

@testset "DroughtIndices.jl" begin
    # Write your tests here.
    @test floor(thornthwaite([1,2,3], 1)[1]) == 59
    @test floor(moms([1,2,3], 3)[1][1]) == 2.0
    @test params([1,2,3], [0.5, 0.3, 0.2])[3] == -.2
    @test floor(cdfglo([4,1],[1.2,1,0.2])[1]) == 2
end

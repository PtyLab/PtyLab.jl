@testset "Params" begin
    p = Params()
    @test p.fftshiftFlag == false
    @test p.propagatorType == Fraunhofer


    # add more for missing
end

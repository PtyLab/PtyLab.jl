@testset "Params" begin
    p = Params()
    @test p.fftshiftSwitch == false
    @test p.propagatorType == Fraunhofer


    # add more for missing
end

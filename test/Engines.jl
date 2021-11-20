@testset "Engines" begin

    @testset "DataType Checks" begin
        e = ePIE()
        @test typeof(e) <: PtyLab.Engines
        @test e.betaProbe == 0.25f0
        @test e.betaObject == 0.25f0
        @test e.numIterations == 50 
    end


    @testset "Constraints" begin
        @test PtyLab.centerOfMassStabilizationOffset([0 0 0 0; 0 0 0 0; 0 0 1 0;  0 0 0 0]) == (3,3) 
        @test PtyLab.centerOfMassStabilizationOffset([0 1 0 0; 0 0 0 0; 0 0 0 0;  0 0 0 0]) == (1, 2)
    end

end

@testset "Test parameter dtype setting" begin
    @test typeof(initializeExperimentalData()) === PtyLab.ExperimentalData{Float32} 
    @test typeof(initializeExperimentalData(dtype=Float64)) === PtyLab.ExperimentalData{Float64} 
end

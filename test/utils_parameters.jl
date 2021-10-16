@testset "Test parameter setting" begin
    @test typeof(init_Parameters()) === Parameters{Float32} 
end

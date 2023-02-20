@testset "Utils display" begin
    grr = grid_regular_rand((4,4), (2,2), (2,2), 0)

    @test show_grid(grr) == Gray{Float64}[Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(0.0); Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(0.0); Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(1.0) Gray{Float64}(0.0); Gray{Float64}(0.0) Gray{Float64}(0.0) Gray{Float64}(0.0) Gray{Float64}(0.0)]

end

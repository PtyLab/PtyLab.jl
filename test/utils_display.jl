@testset "Utils display" begin
    @test PtyLab.complex_show([1 + 0im, 0 + 0im, -1im, 1im]) == ColorTypes.HSV{Float64}[HSV{Float64}(0.0,1.0,1.0), HSV{Float64}(0.0,1.0,0.0), HSV{Float64}(-64.0,1.0,1.0), HSV{Float64}(64.0,1.0,1.0)]
    @test PtyLab.gray_show([1, 0]) == ColorTypes.Gray{Float64}[Gray{Float64}(1.0), Gray{Float64}(0.0)]
    @test PtyLab.gray_show([2, -1], set_one = true) == ColorTypes.Gray{Float64}[Gray{Float64}(1.0), Gray{Float64}(-0.5)]
    @test PtyLab.gray_show([1, -1], set_zero = true) == ColorTypes.Gray{Float64}[Gray{Float64}(1.0), Gray{Float64}(0.0)]
    @test PtyLab.gray_show([1, -1], set_one = true) == ColorTypes.Gray{Float64}[Gray{Float64}(1.0), Gray{Float64}(-1.0)]
end

@testset "Utils" begin
    @testset "Random positions order" begin
        arr = [1 2 3 4 5; 1 2 3 4 5]
        out, r = PtyLab.randpermPositions(arr)
        @test size(out) == size(arr)
        @test typeof(out) == typeof(arr)
        
        @test in(1, out[2, :])
        @test in(2, out[2, :])
        @test in(3, out[2, :])
        @test in(4, out[2, :])
        @test in(5, out[2, :])
        @test in(1, out[1, :])
        @test in(2, out[1, :])
        @test in(3, out[1, :])
        @test in(4, out[1, :])
        @test in(5, out[1, :])
    end
end

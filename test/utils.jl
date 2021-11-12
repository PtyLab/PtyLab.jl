@testset "Utils" begin

    @testset "circ" begin
        @test PtyLab.circ((1, 5), 2, 0.5) == Bool[0 0 1 0 0]
        @test PtyLab.circ((1, 5), 2, 1.01) == Bool[0 1 1 1 0]
        @test PtyLab.circ((1, 5), 2, 1) == Bool[0 0 1 0 0]
        @test PtyLab.circ((5,), 2, 1) == Bool[0, 0, 1, 0, 0]
    end


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
    

    @testset "Maybe fftshift" begin
        maybe_fftshift, maybe_ifftshift = PtyLab.getMaybeFftshifts(true)
        @test maybe_fftshift([1, 2, 3, 4]) == [3, 4, 1, 2]
        @test maybe_ifftshift([1, 2, 3, 4]) == [3, 4, 1, 2]
        @test maybe_fftshift([1, 2, 3, 4, 5]) == [4, 5, 1, 2, 3]
        @test maybe_ifftshift([-1, -2, -3, 10, 20]) == [-3, 10, 20, -1, -2]

        maybe_fftshift, maybe_ifftshift = PtyLab.getMaybeFftshifts(false)
        @test maybe_fftshift([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]
        @test maybe_ifftshift([1, 2, 3, 4, -1]) == [1, 2, 3, 4, -1]
    end
end

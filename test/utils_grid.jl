@testset "Grid Regular Random" begin

    @testset "Random wiggle" begin
        ks = 1:10:101
        for i = 1:20
            @test 1 <= PtyLab.random_constrained_wiggle(1, ks, 2)  <= 3
            @test 99 <= PtyLab.random_constrained_wiggle(101, ks, 2) <= 101
            @test 48 <= PtyLab.random_constrained_wiggle(50, ks, 2) <= 52  
        end
    end

    @testset "grid regular rand" begin
        grid_size = (120, 100)
        tile_size = (30, 20)
        N, M = (7, 9)
        grr = PtyLab.grid_regular_rand(grid_size, tile_size, (N, M))
 
        for t in grr.tiles
            @test t.i₁ >= 1
            @test t.j₁ >= 1
            @test t.i₂ <= grid_size[1]
            @test t.j₂ <= grid_size[2]
            @test (-t.i₁ + t.i₂ + 1) == tile_size[1]
            @test (-t.j₁ + t.j₂ + 1) == tile_size[2]
        end
    end


    @testset "calc_overlap" begin
        @test PtyLab.calc_overlap(5, 10) == 0.5
        @test PtyLab.calc_overlap(1, 10) == 0.9
    end
end

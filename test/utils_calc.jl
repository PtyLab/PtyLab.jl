@testset "Utils calc" begin

    @test PtyLab.calc_Nd(zeros((5,3))) == 5
    # @test PtyLab.calc_xd(10, 0.1) == []
    # @test PtyLab.calc_xo(10, 0.1) == []
    # @test PtyLab.calc_xp(10, 0.1) == []

    @test PtyLab.calc_Ld(42, 0.13) == 5.46
    @test PtyLab.calc_Lp(42, 0.13) == 5.46
    @test PtyLab.calc_Lo(42, 0.13) == 5.46
    @test PtyLab.calc_Nd(zeros((5,3))) == 5


    # @test PtyLab.calc_Np(10) == 10
    # @test PtyLab.calc_No(10) == 40


    @test PtyLab.calc_dxo(1e-5) == 1e-5

    @test PtyLab.calc_numFrames(randn((11,))) == 11
    
    x = zeros((2,2,2))
    x[:, :, 1] = [1 1; 0  -1]
    x[:, :, 2] = [0 0; 0 2]
    @test PtyLab.calc_energyAtPos(x) ≈ [3, 2] 

    x = zeros((2,2,2))
    x[:, :, 1] = [1 1; 0  -1]
    x[:, :, 2] = [0 0; 0 2]
    @test PtyLab.calc_maxProbePower(x) ≈ 3 


    @test PtyLab.calc_NAd(10e-3, 50e-3) ≈ 0.1
    @test PtyLab.calc_DoF(550e-9, 0.1) ≈ 55e-6


end

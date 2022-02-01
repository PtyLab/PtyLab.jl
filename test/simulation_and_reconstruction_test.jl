function simulate()
    img_abs = Float32.(testimage("fabio_gray_512"))[200:380, 200:380]
    img_phase = Float32.(testimage("resolution_test_512"))[200-30:380-30, 200+30:380+30]
    object = img_abs .* cispi.(2 .* img_phase)

    grid_size = size(object)
    tile_size = (70, 70)
    grr = PtyLab.grid_regular_rand(grid_size, tile_size, (18, 18), 30);

    probe = IndexFunArrays.gaussian(Float32, tile_size, scale=0.010) .* cis.(Float32(2π) .* 
     4 .* gaussian(Float32, tile_size, scale=0.003));

    ptychogram = zeros(Float32, (tile_size..., length(grr.tiles)));
    p = Params()
    o2d, d2o = Fraunhofer(probe, fftshiftFlag=true);
    
    for (i, t) in enumerate(grr.tiles)
        ptychogram[:, :, i] = abs2.(o2d(view(object, t.i₁:t.i₂,  t.j₁:t.j₂) .* probe))#, 200000)
    end
    
    
    lambda = 633f-9
    z = 50f-3
    dxd = 10f-6
    Nd = size(ptychogram, 1)
    dxo = lambda * z / (Nd * dxd)
    
    fid_new = h5open("simulated_ptychography.hdf5", "w");
    fid_new["Nd"] = Nd
    fid_new["No"] = size(img_abs, 1)
    fid_new["dxd"] = 10f-6
    fid_new["encoder"] = PtyLab.encoder(grr, dxo, offset=(50, 50))
    fid_new["wavelength"] = lambda
    fid_new["entrancePupilDiameter"] = dxo * 30
    fid_new["zo"] = z
    fid_new["ptychogram"] = ptychogram
    close(fid_new)

    return object, probe
end


function reconstruct()
    
    experimentalData = ExperimentalDataCPM("simulated_ptychography.hdf5");
    
    reconstruction = ReconstructionCPM(experimentalData);
    reconstruction = PtyLab.initializeObjectProbe!(reconstruction);
    
    engine = PtyLab.ePIE()
    params2 = Params(fftshiftFlag=false, transposePtychogram=false, comStabilizationSwitch=true)

    #engines.
    engine.betaProbe = 0.75f0
    engine.betaObject = 0.75f0
    
    reconstruction = PtyLab.initializeObjectProbe!(reconstruction);
    engine.numIterations = 50
    @time p, o = PtyLab.reconstruct(engine, params2, reconstruction);

    return o, p

 end


@testset "Test forward simulation with reconstruction" begin
    object_gt, probe_gt = simulate()
    object, probe = reconstruct()


    # probe preparation for comparison
    p_1D_1 = real.(probe[:, 35, 1,1,1,1] ./ sum(probe))
    p_gt_1D_1 = real.(probe_gt[:, 35] ./ sum(probe_gt))

    p_1D_2 = real.(probe[30, :, 1,1,1,1] ./ sum(probe))
    p_gt_1D_2 = real.(probe_gt[30, :] ./ sum(probe_gt))



    @test p_1D_1 .+ 10 ≈ p_gt_1D_1 .+ 10
    @test p_1D_2 .+ 10 ≈ p_gt_1D_2 .+ 10

    object_comp = object[50+50:150+50, 100-20:150-10, 1,1,1,1] 
    
    object_comp ./= sum(object_comp)
    
    
    object_comp = real.(object_comp)
    object_comp ./= sum(object_comp)
    
    o1, o2 = (45, 45 - 10)
    
    object_gt_comp = object_gt[o1:o1+100, o2:o2+60]
    object_gt_comp ./= sum(object_gt_comp)
    
    object_gt_comp = real.(object_gt_comp)
    object_gt_comp ./= sum(object_gt_comp)
    
    
    @test all(≈(1 .+ real.(object_gt_comp), 1 .+ real.(object_comp), rtol=2.1e-4))
end

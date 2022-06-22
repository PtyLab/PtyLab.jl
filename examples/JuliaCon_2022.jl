### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5efb56cc-81c7-11ec-2b28-e71f25e0d443
begin
	using Pkg;
	Pkg.activate(".")
end

# ╔═╡ 8b2e821a-cde0-4eba-ba1b-b176456cf667
begin
	using Revise, EllipsisNotation, ImageShow
	using PtyLab, TestImages, IndexFunArrays, FFTW, HDF5, Noise, FourierTools, CUDA
	using PlutoUI, NDTools, Plots
end

# ╔═╡ 520b96b2-e35d-4820-8d4b-fbedd962fd92
md"## PtyLab.jl - Ptychography Reconstruction

Julia implementation of PtyLab.

### Felix Wechsler
* PhD student at the IPHT in Jena (Rainer Heintzmann's Lab)
* [felixwechsler.science](https://felixwechsler.science)

### Dr. Lars Loetgering
* Research scientist at Carl Zeiss AG
"

# ╔═╡ 89773bf9-3c7d-400a-9cef-46d6d97a4635
md"
## What is Ptychography?

Schick mir ein Bild, dann kann ich es einfügen

"

# ╔═╡ db2c45d1-a86b-428d-a696-df01637d4ad9


# ╔═╡ 0f1f8cfa-477c-4133-b2c4-0e39138893d1
md"## PtyLab.jl
* Implementation of Iterative Ptychographic Engine (PIE)
* CUDA.jl support
* Flexible type hierarchy to add more solvers
"

# ╔═╡ ba3e98c4-8e25-405c-b342-e2aae01ec045
md"## Allocation free style of programming

* the PIE engine needs often $≈ 50 - 300$ of iterations to converge
* the arrays are 6 dimensions but their shapes can vary dynamically on the problem
* we need at lot (!) memory buffers to handle
* To not confuse the buffers we use a functional style of programming

For example, the operation *detector to object* and *object to detector* is required in each iteration.

We create new functions which capture the outer variables.
Because of [#15276](https://github.com/JuliaLang/julia/issues/15276) we use `let` blocks.
"

# ╔═╡ a0545e2c-f666-4bf3-af67-336ccc0efe1d
function Fraunhofer(arr::T) where T
	ss = sqrt(size(arr, 1) * size(arr, 2))
	p = plan_fft!(arr)
	
    o2d! = let  p=p
				buffer_shift = similar(arr)
        function o2d!(x)
            p * x
            x ./= ss
			fftshift!(buffer_shift, x)
        end
    end

    d2o! = let  p_inv=inv(p)
		 		buffer_shift = similar(arr)
        function d2o!(x)
            p_inv * ifftshift!(buffer_shift, x)
            buffer_shift .*= ss
        end
    end

	return o2d!, d2o!
end

# ╔═╡ a9499a36-a3a3-41d1-ae52-0f18db1de62a
arr = randn(ComplexF32, (1024, 1024));

# ╔═╡ e887ec5e-75ac-4d62-a84f-8e1f4e96bc98
o2d_, d2o_ = Fraunhofer(arr);

# ╔═╡ 798ea080-3a25-452f-880b-ae9ecc57e16a
# ╠═╡ show_logs = false
@code_warntype o2d_(arr);

# ╔═╡ 87f623ca-c761-4a54-8eea-ac06c10ea7d1
md"### Usage
The new functions use the buffers and precalculated objects.

The function call is allocation free!
"

# ╔═╡ 48081a97-4106-4ec5-84f7-410c3c335672
@time o2d_(d2o_(arr));

# ╔═╡ 53a974c3-82b0-4195-a777-c6e2ec8af7cd
arr ≈ o2d_(d2o_(copy(arr)))

# ╔═╡ d2ce7638-9765-4752-b4b5-953ef7d43003
md"# Load a Testimage
Let's assemble an image with an absolute intensity and a phase term.
Absolute value is _Fabio_ whereas the phase term is a siemens star and some other features.
"

# ╔═╡ 83dd43ee-aa6d-414d-b4a2-dec98afc8aa1
begin
	img_abs = Float32.(testimage("fabio_gray_512"))[200:400, 200:400]
	img_phase = Float32.(testimage("resolution_test_512"))[140:340, 250:450]

	object = img_abs .* cispi.(2 .* img_phase)
	complex_show(object);
end

# ╔═╡ b0c38140-9e58-4f55-83fb-2e0632a83c63


# ╔═╡ a86151ed-c915-465b-af67-dee7a552e58c
md"## Make a probe
Let's choose a gaussian probe despite it's not very optimal for Ptychography

"

# ╔═╡ 3394318c-7b7c-4720-8822-f7cc1bf3645f
#tile_size = (250, 250)
tile_size = (80, 80)

# ╔═╡ 73b8c142-df52-4ed2-b35d-3d9f1b1ae707
begin
	# scale=0.0008
	# scale =0.008
	probe = IndexFunArrays.gaussian(Float32, tile_size, scale=0.008) .* cis.(Float32(2π) .* 
	      rr2(Float32, tile_size, scale=0.045));
	complex_show(probe)
end

# ╔═╡ 3860c141-45cd-4035-a8bd-a51a67463c46
md"# Scanning Grid
PtyLab.jl supports only regular grid with random disturbance.
Also not very optimal, but is simple to start with.
"

# ╔═╡ 130d273c-9ad3-41f0-bee9-737136e9fd35
begin
	grid_size = size(object)
	
	grr = PtyLab.grid_regular_rand(grid_size, tile_size, (14, 14), 10);
	show_grid(grr, only_points=true, thickness=1)
end

# ╔═╡ 2c9c4ce7-2efa-4086-ab65-f0ee98b5404b
grr.overlap

# ╔═╡ b7da887b-f826-4eb2-8d8a-80c265ec1eed
md"# Simulate dataset

Procedure:

* extract a tile of the object which fits the size of the probe
* multiply with probe
* apply Fraunhofer propagation on the field
* `abs2` for intensity

"

# ╔═╡ a129a90c-cfd7-4a8a-8c37-6a5a5f4aadc2
md"
For performance, it's better to choose `Float32` (4 Byte) over `Float64` (8 Byte).
Especially on GPUs, that's critical since they perform poor with `Float64`.

"

# ╔═╡ 095f67cd-04be-431b-a51e-ec4e6a5e16cb
ptychogram = zeros(Float32, (tile_size..., length(grr.tiles)));

# ╔═╡ 3ab61039-acb1-4994-8d5c-88b1613700bf
size(ptychogram)

# ╔═╡ 4c47631a-1a7f-46f1-8d11-fb3d5d2f99de
# object2detector, detector2object
o2d, d2o = PtyLab.Fraunhofer(probe, fftshiftFlag=true);

# ╔═╡ 6cb98beb-eeaf-48d2-9dbd-70e81015f09b
md"## Simulate!"

# ╔═╡ f276663f-21b7-45e6-97f6-4b0b841dcdde
for (i, t) in enumerate(grr.tiles)
	tile = view(object, t.i₁:t.i₂,  t.j₁:t.j₂)
    ptychogram[:, :, i] = poisson(abs2.(o2d(tile .* probe)), 1000)
end

# ╔═╡ 5b7f8494-a417-4f0c-a9c3-5fd2a87df673
@bind slice Slider(1:length(grr.tiles), show_value=true)

# ╔═╡ 5f1d8f65-655d-4854-b820-be834ab3b1f6
gray_show(ptychogram[:, :, slice])

# ╔═╡ 7193a29e-a6a5-4c60-b68b-dc4446121bc5
Plots.heatmap(ptychogram[:, :, slice]);

# ╔═╡ 2622d993-5ae9-4072-9775-0ac448029bc6
md"# Storing the dataset using HDF5
That format is compatible to PtyLab.m and PtyLab.py.
"

# ╔═╡ 3dceb274-7a87-426b-92a3-8cc319d442e7
begin
	lambda = 633f-9
	z = 50f-3
	dxd = 10f-6
	Nd = size(ptychogram, 1)
	dxo = lambda * z / (Nd * dxd)
	
	fid = h5open("simulated_ptychography.hdf5", "w");
	fid["Nd"] = Nd
	fid["No"] = size(img_abs, 1)
	fid["dxd"] = 10f-6
	fid["encoder"] = PtyLab.encoder(grr, dxo, offset=(50, 50))
	fid["wavelength"] = lambda
	fid["entrancePupilDiameter"] = dxo * Nd / 2
	fid["zo"] = z
	fid["ptychogram"] = ptychogram
	close(fid)
end

# ╔═╡ 5a750501-c58a-4566-9981-ec0bd71cd50e


# ╔═╡ 5d026111-f8e3-4e36-9d81-2970c6e74ac1
md"# Load Dataset
From an `.hdf5` dataset. Compatible within the PtyLab ecosystem.

"

# ╔═╡ d26dc7af-d4a9-45c9-9929-0017d9c82a6c
experimentalData = ExperimentalDataCPM("simulated_ptychography.hdf5");

# ╔═╡ 1d91d409-6307-4bd7-a4a1-c26cadd009c5
md"
Create a `reconstruction` struct containing all the physical parameters.

`ReconstructionCPM` is a constructor creating a struct containing all data needed for reconstruction. This struct is mutable (e.g. change `z`).
"

# ╔═╡ a17b1514-30dd-4dc4-b0cf-d9d83e4a8466
reconstruction = ReconstructionCPM(experimentalData, cuda=false);

# ╔═╡ 80881c2c-15c5-4dae-9b5c-04b0817be995
md"
Initialize the probe
"

# ╔═╡ 21759a3d-8b59-453d-9432-21b33a63f563
PtyLab.initializeObjectProbe!(reconstruction);

# ╔═╡ 550d4838-da12-41b7-be8e-0c9a4ed0497a
complex_show(Array(reconstruction.probe)[:, :, 1,1,1,1])

# ╔═╡ b20cd8d6-6ca9-47b4-84af-6154e2bc313f
reconstruction;

# ╔═╡ 1a74c5b8-c92b-4113-a4d6-256936988533
md"# Params for reconstruction"

# ╔═╡ 2c6563e4-d848-43c2-b123-fe32649401f9
params = Params(transposePtychogram = false,
				comStabilizationSwitch = true);

# ╔═╡ db12f8ed-afbb-4205-9542-bd99f5f60e4a
params

# ╔═╡ aac02813-e0f7-4b57-9df7-e2e768c8dbac
md"# Select an Engine"

# ╔═╡ 8d395f98-3138-4084-b96e-805ab12335a3
engine = PtyLab.ePIE(numIterations=50, betaObject=0.75f0, betaProbe=0.75f0)

# ╔═╡ db005e37-ed08-4964-a806-e53bdc09f527


# ╔═╡ d18269ae-2006-4ba9-b81e-edb2f504f208
md"# Run the reconstruction"

# ╔═╡ 50c991b3-02b6-48ae-8eeb-72d9321a9990
@time p, o = PtyLab.reconstruct(engine, params, reconstruction);

# ╔═╡ b7d2a8da-7f5c-4edd-b80d-59fef05bf692
complex_show(select_region(Array(o)[:, :, 1,1,1,1], new_size=size(img_abs)))

# ╔═╡ e04e73a2-28b0-4594-8a5e-eb97ded5b1a4
complex_show(Array(p)[:, :, 1,1,1,1])

# ╔═╡ ec87b2bd-d812-4cde-8cf2-70f16dc600cc
md"Ptychogram size: $(round(sizeof(reconstruction.ptychogram) / 2^20, digits=2)) MiB"

# ╔═╡ 9effcfc6-bf47-4b39-bfa3-b94959f843d6


# ╔═╡ 748d5712-4609-4446-ac32-d45359d26000


# ╔═╡ Cell order:
# ╠═5efb56cc-81c7-11ec-2b28-e71f25e0d443
# ╠═8b2e821a-cde0-4eba-ba1b-b176456cf667
# ╟─520b96b2-e35d-4820-8d4b-fbedd962fd92
# ╟─89773bf9-3c7d-400a-9cef-46d6d97a4635
# ╠═db2c45d1-a86b-428d-a696-df01637d4ad9
# ╟─0f1f8cfa-477c-4133-b2c4-0e39138893d1
# ╟─ba3e98c4-8e25-405c-b342-e2aae01ec045
# ╠═a0545e2c-f666-4bf3-af67-336ccc0efe1d
# ╠═a9499a36-a3a3-41d1-ae52-0f18db1de62a
# ╠═e887ec5e-75ac-4d62-a84f-8e1f4e96bc98
# ╠═798ea080-3a25-452f-880b-ae9ecc57e16a
# ╟─87f623ca-c761-4a54-8eea-ac06c10ea7d1
# ╠═48081a97-4106-4ec5-84f7-410c3c335672
# ╠═53a974c3-82b0-4195-a777-c6e2ec8af7cd
# ╟─d2ce7638-9765-4752-b4b5-953ef7d43003
# ╠═83dd43ee-aa6d-414d-b4a2-dec98afc8aa1
# ╠═b0c38140-9e58-4f55-83fb-2e0632a83c63
# ╟─a86151ed-c915-465b-af67-dee7a552e58c
# ╠═3394318c-7b7c-4720-8822-f7cc1bf3645f
# ╠═73b8c142-df52-4ed2-b35d-3d9f1b1ae707
# ╟─3860c141-45cd-4035-a8bd-a51a67463c46
# ╠═130d273c-9ad3-41f0-bee9-737136e9fd35
# ╠═2c9c4ce7-2efa-4086-ab65-f0ee98b5404b
# ╟─b7da887b-f826-4eb2-8d8a-80c265ec1eed
# ╟─a129a90c-cfd7-4a8a-8c37-6a5a5f4aadc2
# ╠═095f67cd-04be-431b-a51e-ec4e6a5e16cb
# ╠═3ab61039-acb1-4994-8d5c-88b1613700bf
# ╠═4c47631a-1a7f-46f1-8d11-fb3d5d2f99de
# ╟─6cb98beb-eeaf-48d2-9dbd-70e81015f09b
# ╠═f276663f-21b7-45e6-97f6-4b0b841dcdde
# ╟─5b7f8494-a417-4f0c-a9c3-5fd2a87df673
# ╠═5f1d8f65-655d-4854-b820-be834ab3b1f6
# ╠═7193a29e-a6a5-4c60-b68b-dc4446121bc5
# ╟─2622d993-5ae9-4072-9775-0ac448029bc6
# ╠═3dceb274-7a87-426b-92a3-8cc319d442e7
# ╠═5a750501-c58a-4566-9981-ec0bd71cd50e
# ╟─5d026111-f8e3-4e36-9d81-2970c6e74ac1
# ╠═d26dc7af-d4a9-45c9-9929-0017d9c82a6c
# ╟─1d91d409-6307-4bd7-a4a1-c26cadd009c5
# ╠═a17b1514-30dd-4dc4-b0cf-d9d83e4a8466
# ╟─80881c2c-15c5-4dae-9b5c-04b0817be995
# ╠═21759a3d-8b59-453d-9432-21b33a63f563
# ╠═550d4838-da12-41b7-be8e-0c9a4ed0497a
# ╠═b20cd8d6-6ca9-47b4-84af-6154e2bc313f
# ╟─1a74c5b8-c92b-4113-a4d6-256936988533
# ╠═2c6563e4-d848-43c2-b123-fe32649401f9
# ╠═db12f8ed-afbb-4205-9542-bd99f5f60e4a
# ╟─aac02813-e0f7-4b57-9df7-e2e768c8dbac
# ╠═8d395f98-3138-4084-b96e-805ab12335a3
# ╠═db005e37-ed08-4964-a806-e53bdc09f527
# ╟─d18269ae-2006-4ba9-b81e-edb2f504f208
# ╠═50c991b3-02b6-48ae-8eeb-72d9321a9990
# ╠═b7d2a8da-7f5c-4edd-b80d-59fef05bf692
# ╠═e04e73a2-28b0-4594-8a5e-eb97ded5b1a4
# ╟─ec87b2bd-d812-4cde-8cf2-70f16dc600cc
# ╠═9effcfc6-bf47-4b39-bfa3-b94959f843d6
# ╠═748d5712-4609-4446-ac32-d45359d26000

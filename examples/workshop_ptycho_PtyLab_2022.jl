### A Pluto.jl notebook ###
# v0.17.7

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
md"# PtyLab.jl

Julia implementation of PtyLab.

## Felix Wechsler
* PhD student at the IPHT in Jena (Rainer Heintzmann's Lab)
* [felixwechsler.science](https://felixwechsler.science)


## What is [Julia](https://julialang.org/)?

Relatively new dynamic, high level programming language

* 1.0 released in 2018
* high level and dynamic (like Matlab, Python, ...)
* just-in-time (JIT) compilation
* core paradigm: multiple dispatch
* from beginning on designed for performance
* first class support for CUDA (AMDGPU works partially)
* we use [Pluto.jl](https://github.com/fonsp/Pluto.jl) as notebook

## PtyLab.jl
Current state

* conventional Ptychography (however, interfaces are flexible to add Fourier Ptychography quickly)
* CPU support mostly single threaded
* CUDA support for highly parallelized execution (1 GPU)
* ePIE
* multiple modes
* focus on performance
* currently around 400 LLOC
* flexible interface to add more functionality
* tested with GitHub actions and Codecov


"

# ╔═╡ 13e8b8ed-3297-41dd-a0b7-b993f669fe47
md" # Basics Structure
Similar to the Python structure

## `ExperimentalData` 
* struct containing the data loaded from a dataset

## `Reconstruction` 
* struct containing the data used within the reconstruction
* those fields are mutable!

## `Engine` 
* containing the data for the algorithm such as `numIterations`.

## `Params`
* containing meta information for the reconstruction 
* such as: `transposePtychogram` or `propagatorType`

## Visualization
* currently very rudimentary
* In long-term based on [napari](https://github.com/napari/napari).
* napari has no support for complex numbers yet

"

# ╔═╡ ff799f59-38d8-4b6d-a0ed-5c1b6f4059a1
md"# Load Packages

"

# ╔═╡ d2ce7638-9765-4752-b4b5-953ef7d43003
md"# Load a Testimage
Let's assemble an image with an absolute intensity and a phase term.
Absolute value is _Fabio_ whereas the phase term is a siemens star and some other features.
"

# ╔═╡ 83dd43ee-aa6d-414d-b4a2-dec98afc8aa1
begin
	img_abs = Float32.(testimage("fabio_gray_512"))[200:400, 200:400]
	img_phase = Float32.(testimage("resolution_test_512"))[140:340, 250:450]
	
	#img_abs = Float32.(testimage("fabio_gray_512"))
	#img_phase = Float32.(testimage("resolution_test_512"))
	
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

# ╔═╡ d7bafd2b-313d-41b1-8782-4789258fd29b
Plots.plot(abs.(probe[:, 51]))

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

* extract a tile of the object
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

# ╔═╡ e49cd149-eb2d-4b6e-9530-da00b6d79826
md"### Functional Style of programming
`Fraunhofer` assembles two functions `object2detector` and `detector2object` which (under the hood) store the phase functions, the FFT plan, buffer, etc.

Hence, an application works completely allocation free and very efficiently.

An iterative reconstruction, such as ePIE, needs very efficient functions in the hot for loops. Therefore, allocate for all operations memory buffers. That prevents, that the language has to allocate and free memory from time to time. Also GPUs work more efficiently if used with buffers.
"

# ╔═╡ 4c47631a-1a7f-46f1-8d11-fb3d5d2f99de
# object2detector, detector2object
object2detector, detector2object = Fraunhofer(probe, fftshiftSwitch=true);

# ╔═╡ 6cb98beb-eeaf-48d2-9dbd-70e81015f09b
md"## Simulate!"

# ╔═╡ f276663f-21b7-45e6-97f6-4b0b841dcdde
for (i, t) in enumerate(grr.tiles)
	tile = view(object, t.i₁:t.i₂,  t.j₁:t.j₂)
    ptychogram[:, :, i] = poisson(abs2.(object2detector(tile .* probe)), 1000)
end

# ╔═╡ 5b7f8494-a417-4f0c-a9c3-5fd2a87df673
@bind slice Slider(1:length(grr.tiles), show_value=true)

# ╔═╡ 5f1d8f65-655d-4854-b820-be834ab3b1f6
gray_show(ptychogram[:, :, slice])

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

# ╔═╡ 5045ea4a-33c5-4741-854e-d9c379450259
params

# ╔═╡ aac02813-e0f7-4b57-9df7-e2e768c8dbac
md"# Select an Engine"

# ╔═╡ 8d395f98-3138-4084-b96e-805ab12335a3
engine = PtyLab.ePIE(numIterations=50, betaObject=0.95f0, betaProbe=0.95f0)

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
# ╟─520b96b2-e35d-4820-8d4b-fbedd962fd92
# ╟─13e8b8ed-3297-41dd-a0b7-b993f669fe47
# ╠═5045ea4a-33c5-4741-854e-d9c379450259
# ╟─ff799f59-38d8-4b6d-a0ed-5c1b6f4059a1
# ╠═8b2e821a-cde0-4eba-ba1b-b176456cf667
# ╟─d2ce7638-9765-4752-b4b5-953ef7d43003
# ╠═83dd43ee-aa6d-414d-b4a2-dec98afc8aa1
# ╠═b0c38140-9e58-4f55-83fb-2e0632a83c63
# ╟─a86151ed-c915-465b-af67-dee7a552e58c
# ╠═3394318c-7b7c-4720-8822-f7cc1bf3645f
# ╠═73b8c142-df52-4ed2-b35d-3d9f1b1ae707
# ╠═d7bafd2b-313d-41b1-8782-4789258fd29b
# ╟─3860c141-45cd-4035-a8bd-a51a67463c46
# ╠═130d273c-9ad3-41f0-bee9-737136e9fd35
# ╠═2c9c4ce7-2efa-4086-ab65-f0ee98b5404b
# ╟─b7da887b-f826-4eb2-8d8a-80c265ec1eed
# ╟─a129a90c-cfd7-4a8a-8c37-6a5a5f4aadc2
# ╠═095f67cd-04be-431b-a51e-ec4e6a5e16cb
# ╠═3ab61039-acb1-4994-8d5c-88b1613700bf
# ╟─e49cd149-eb2d-4b6e-9530-da00b6d79826
# ╠═4c47631a-1a7f-46f1-8d11-fb3d5d2f99de
# ╟─6cb98beb-eeaf-48d2-9dbd-70e81015f09b
# ╠═f276663f-21b7-45e6-97f6-4b0b841dcdde
# ╟─5b7f8494-a417-4f0c-a9c3-5fd2a87df673
# ╠═5f1d8f65-655d-4854-b820-be834ab3b1f6
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

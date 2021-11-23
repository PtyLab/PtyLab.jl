module PtyLab

using MicroscopyTools

using NDTools

 # more fancy struct filling
using Parameters

 # we need fft and variants
using FFTW, FourierTools
FFTW.set_num_threads(4)
 # for displaying
using ColorTypes
 # data loading
using HDF5 
 # makes .. possible of array
using EllipsisNotation

 # stuff like rr2,... 
using IndexFunArrays

 # for randperm
using Random

using LoopVectorization

 #
using Tullio

 # progress bar for loops
using ProgressMeter

 # CUDA acceleration
using CUDA

 # @kwdef macro
import Base.@kwdef

 # some constants
const DEFAULT_WAVELENGTH = 632.8e-9 :: Float64


 # a few types for either Fourier or Classical Ptychography
abstract type ModePtychograpy end
abstract type CPM <: ModePtychograpy end
abstract type FPM <: ModePtychograpy end

include("utils_grid.jl")


 # basic parts
include("ExperimentalData.jl")

 # types to store physical quantities for reconstruction
include("Reconstruction.jl")

 # parameters for reconstruction
include("Params.jl")

 # operators for field propagation
include("Operators.jl")

 # different reconstruction engines
include("Engines.jl")
include("Engines/ePIE.jl")


 # some utilities
include("utils.jl")
include("utils_calc.jl")
include("utils_display.jl")

end

module PtyLab

 # more fancy struct filling
using Parameters

 # we need fft and variants
using FFTW, FourierTools
 # for displaying
using ColorTypes
 # data loading
using HDF5 
 # makes .. possible of array
using EllipsisNotation

 # stuff like rr2,... 
using IndexFunArrays

 # @kwdef macro
import Base.@kwdef

 # some constants
const DEFAULT_WAVELENGTH = 632.8e-9 :: Float64


 # a few types for either Fourier or Classical Ptychography
abstract type ModePtychograpy end
abstract type CPM <: ModePtychograpy end
abstract type FPM <: ModePtychograpy end

include("utils_grid.jl")


include("ExperimentalData.jl")
include("Reconstruction.jl")
include("Params.jl")
include("Operators.jl")
include("Engines.jl")

include("utils.jl")
include("utils_calc.jl")
include("utils_display.jl")

end

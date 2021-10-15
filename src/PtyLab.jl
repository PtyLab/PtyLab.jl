module PtyLab

 # we need fft and variants
using FFTW, FourierTools
 # for displaying
using ColorTypes


 # some constants
const DEFAULT_WAVELENGTH = 632.8e-9 :: Float64


include("utils_parameters.jl")
include("utils_grid.jl")
include("utils_display.jl")


end

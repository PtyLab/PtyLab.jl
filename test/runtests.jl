using PtyLab 
using Test
using ColorTypes
using FFTW
using IndexFunArrays
using TestImages
using HDF5
using Random
Random.seed!(42)

 # tests
include("utils.jl")
include("utils_grid.jl")
include("utils_display.jl")
include("utils_calc.jl")
include("Operators.jl")
include("Params.jl")
include("Engines.jl")


 # final simulation and reconstruction comparison!
include("simulation_and_reconstruction_test.jl")

# include("ExperimentalData.jl")

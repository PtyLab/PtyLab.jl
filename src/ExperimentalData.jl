export ExperimentalDataCPM

"""
    ExperimentalData

Abstract type for [`ExperimentalDataCPM`](@ref ExperimentalDataCPMa) 
and [`ExperimentalDataCPM`](@ref ExperimentalDataFPM)
"""
abstract type ExperimentalData{T} end


"""
    ExperimentalDataCPM{T}

`mutable struct` stores all experimental parameters for a CPM dataset
"""
@kwdef mutable struct ExperimentalDataCPM{T} <: ExperimentalData{T}
    ptychogram::Union{Nothing, Array{T, N}} where N
    numFrames::Int
    energyAtPos::Vector{T}
    maxProbePower::T
    wavelength::T
    encoder::Union{Nothing, Array{T, 2}}
    # detector sampling
    dxd::T
    Nd::Int
    xd::Vector{T}
    Ld::T
    # distance to detector 
    zo::T
    # optional parameters
    entrancePupilDiameter::Union{Nothing, T}
    spectralDensity::Union{Nothing, T}
    theta::Union{Nothing, T}
end

"""
    ExperimentalDataFPM{T}

`mutable struct` stores all parameters for a FPM dataset

"""
@kwdef mutable struct ExperimentalDataFPM{T} <: ExperimentalData{T}
end


"""
    ExperimentalDataCPM(fileName::String)

Fill the `ExperimentalDataCPM` struct with data from a `*.hdf5` file.
"""
function ExperimentalDataCPM(fileName::String)
    # open h5 file
    fid = HDF5.h5open(fileName)
    @info "Reading $fid was successful"
    # function to extract number from a vector wrap
    r_number(x) = haskey(fid, x) ? read(fid, x)[begin] : nothing
    # read arrays, take full
    r_array(x) = haskey(fid, x) ? read(fid, x) : nothing

    # we work here with a dict since that makes passing
    # the variables around more convenient
    d = Dict()
    # data
    d[:ptychogram] = r_array("ptychogram")
    d[:encoder] = r_array("encoder")

    # physics
    d[:wavelength] =r_number("wavelength")
   
    # object
    d[:zo] = r_number("zo")
    # detector
    d[:dxd] = r_number("dxd")

    # optional
    d[:entrancePupilDiameter] = r_number("entrancePupilDiameter")
    d[:spectralDensity] = r_number("spectralDensity")
    d[:theta] = r_number("theta")

    # close fid!
    close(fid)

    # derived properties
    d[:Nd] = calc_Nd(d[:ptychogram])
    d[:xd] = calc_xd(d[:Nd], d[:dxd])
    d[:Ld] = calc_Ld(d[:Nd], d[:dxd])
    d[:numFrames] = calc_numFrames(d[:ptychogram])
    d[:energyAtPos] = calc_energyAtPos(d[:ptychogram])
    d[:maxProbePower] = calc_maxProbePower(d[:ptychogram])

    # create ExperimentalDataCPM struct
    expData = ExperimentalDataCPM(;d...)
    return expData
end


function Base.getproperty(expData::ExperimentalData, sym::Symbol)
    return getfield(expData, sym)
end

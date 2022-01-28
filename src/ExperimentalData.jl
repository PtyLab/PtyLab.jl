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
    # if there are multiple sensors, the size in each dimensions
    # for the smal sensors are listed here
    Nd1::Union{Nothing, Int}
    Nd2::Union{Nothing, Int}
    # positions of the detectors. The position indicates the top left corner
    # of the respective detector
    posDetectors::Union{Nothing, Vector{Tuple{Float32, Float32}}}
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
    ExperimentalDataCPM(fileName::String, [T=Float32])

Fill the `ExperimentalDataCPM` struct with data from a `*.hdf5` file.
`T` converts all float numbers (including array) to type of `T`.
"""
function ExperimentalDataCPM(fileName::String, T=Float32; Nd=nothing)
    # open h5 file
    fid = HDF5.h5open(fileName)
    @info "Reading $fid was successful"
    T_conv(x::Integer) = x 
    T_conv(x::Nothing) = x 
    T_conv(x::T2) where T2 = T == T2 ? x : T(x)
    T_conv(x::AbstractArray{T2, N}) where {T2, N} = T2 == T ? x : T.(x)
    # function to extract number from a vector wrap
    r_number(x) = haskey(fid, x) ? T_conv(read(fid, x)[begin]) : nothing
    # read arrays, take full
    r_array(x) = haskey(fid, x) ? T_conv(read(fid, x)) : nothing

    # we work here with a dict since that makes passing
    # the variables around more convenient
    d = Dict()
    # data
    d[:ptychogram] = r_array("ptychogram")

    d[:encoder] = r_array("encoder")

    # physics
    d[:wavelength] = r_number("wavelength")

    # object
    d[:zo] = r_number("zo")
    # detector
    d[:dxd] = r_number("dxd")
    d[:Nd1] = r_number("Nd1")
    d[:Nd2] = r_number("Nd2")
    d[:posDetectors] = r_array("posDetectors")

    # change to tuple format
    d[:posDetectors] = let
        if !isnothing(d[:posDetectors])
            e = d[:posDetectors]
            [(e[1,i], e[2,i]) for i in axes(e, 2)]
        else
            nothing
        end
    end


    # optional
    d[:entrancePupilDiameter] = r_number("entrancePupilDiameter")
    d[:spectralDensity] = r_number("spectralDensity")
    d[:theta] = r_number("theta")

    # close fid!
    close(fid)
    

    # derived properties
    d[:Nd] = isnothing(Nd) ? calc_Nd(d[:ptychogram]) : Nd
    d[:Ld] = calc_Ld(d[:Nd], d[:dxd])

    # change ptychogram in case of multiple sensors
    if !isnothing(d[:Nd1]) && !isnothing(d[:Nd2])
        d[:Ld] = calc_Ld(d[:Nd], d[:dxd])
        d[:ptychogram] = assembleMultiSensorPtychogram(d[:ptychogram], d[:Nd], d[:Nd1], 
                                                       d[:Nd2], d[:posDetectors], d[:Ld], d[:dxd])
    end




    d[:xd] = calc_xd(d[:Nd], d[:dxd])
    d[:numFrames] = calc_numFrames(d[:ptychogram])
    d[:energyAtPos] = calc_energyAtPos(d[:ptychogram])
    d[:maxProbePower] = calc_maxProbePower(d[:ptychogram])

    # create ExperimentalDataCPM struct
    expData = ExperimentalDataCPM(;d...)
    return expData
end

"""
    assembleMultiSensorPtychogram(ptychogram_gaps, Nd, Nd1, Nd2, posDetectors, Ld, dxd)

In the case of multiple sensors, we assemble here a normal ptychogram from the patches we have.
The gaps are filled with `0`.
"""
function assembleMultiSensorPtychogram(ptychogram_gaps, Nd, Nd1, Nd2, posDetectors, Ld, dxd)
    ptychogram = similar(ptychogram_gaps, Nd, Nd, size(ptychogram_gaps)[end])
    fill!(ptychogram, 0)

    xs = range(-Ld/2, Ld/2, length=size(ptychogram, 1))
    f(x) = 1 + (x - (-Ld / 2)) / Ld * size(ptychogram, 1)
    inds = [round.(Ref(Int), f.(ind)) for ind in posDetectors]
        
    for (i, d) in enumerate(inds)
        ptychogram[d[1]:d[1]+Nd1 - 1, 
                   d[2]:d[2]+Nd2 - 1, ..] .= view(ptychogram_gaps, :, :, i, ..)
    end
    
    return ptychogram
end



function Base.getproperty(expData::ExperimentalData, sym::Symbol)
    return getfield(expData, sym)
end

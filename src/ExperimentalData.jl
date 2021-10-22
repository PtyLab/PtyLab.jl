export ExperimentalDataCPM

"""
    ExperimentalDataCPM{T}

`mutable struct` stores all parameters for a CPM dataset
"""
@kwdef mutable struct ExperimentalDataCPM{T}
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
@kwdef mutable struct ExperimentalDataFPM{T}
    ptychogram::Union{Nothing, Array{T, N}} where N
    wavelength::T
    encoder::Union{Nothing, Array{T, 2}}
    # detector sampling
    dxd::T
    # distance to detector 
    zo::T
    zled::T
    # optional parameters
    magnfication::Union{Nothing, T}
    dxp::Union{Nothing, T}
end


"""
    ExperimentalDataCPM(fileName::String)

Fill the `ExperimentalDataCPM` struct with data from a `*.hdf5` file.
"""
function ExperimentalDataCPM(fileName::String)
    # open h5 file
    fid = HDF5.h5open(fileName)
   
    # function to extract number from a vector wrap
    r_number(x) = haskey(fid, x) ? read(fid, x)[begin] : nothing
    # read arrays, take full
    r_array(x) = haskey(fid, x) ? read(fid, x) : nothing

    # data
    ptychogram = r_array("ptychogram")
    encoder = r_array("encoder")

    # physics
    wavelength=r_number("wavelength")
   
    # object
    zo = r_number("zo")
    # detector
    dxd = r_number("dxd")

    # optional
    entrancePupilDiameter = r_number("entrancePupilDiameter")
    spectralDensity = r_number("spectralDensity")
    theta = r_number("theta")

   

    # derived properties
    Nd = calc_Nd(ptychogram)
    xd = calc_xd(Nd, dxd)
    Ld = calc_Ld(Nd, dxd)
    numFrames = calc_numFrames(ptychogram)
    energyAtPos = calc_energyAtPos(ptychogram)
    maxProbePower = calc_maxProbePower(ptychogram)

    # create ExperimentalDataCPM struct
    expData = ExperimentalDataCPM(;
            ptychogram,
            numFrames,
            energyAtPos,
            maxProbePower,
            wavelength,
            encoder,
            dxd,
            Nd,
            xd,
            Ld,
            zo,
            entrancePupilDiameter,
            spectralDensity,
            theta 
        )

    return expData
end



 # some helper functions
"""
    calc_Nd(ptychogram)

How many pixels per dimension.
"""
calc_Nd(ptychogram) = size(ptychogram, 1)
        

"""
    calc_xd(Nd, dxd)

Calculate detector coordinates in 1D.
"""
calc_xd(Nd, dxd) = typeof(dxd).(dxd .* range(-Nd / 2, Nd / 2, length=Nd))
        

"""
    calc_Ld(Nd, dxd)

Calculate the size of the detector.
"""
calc_Ld(Nd, dxd) = Nd * dxd


"""
    calc_numFrames(ptychogram)

Calculate the number of frames.
"""
calc_numFrames(ptychogram) = size(ptychogram)[end]


"""
    calc_energyAtPos(ptychogram)

Calculate the energy at each position.
"""
calc_energyAtPos(ptychogram) = sum(abs, ptychogram, dims=(1,2))[1, 1, ..]


"""
    calc_maxProbePower(ptychogram)

Calculate max probe power at each position.
"""
calc_maxProbePower(ptychogram) = sqrt(maximum(sum(ptychogram, dims=(1,2))))


"""
    calc_Lp(Np, dxp)

Calculate probe length
"""
calc_Lp(Np, dxp) = Np * dxp



"""
    calc_Np(Nd)

Calculate probe pixel
"""
calc_Np(Nd) = Nd

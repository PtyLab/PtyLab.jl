export ExperimentalData
export initializeExperimentalData

"""
    ExperimentalData{T}

`mutable struct` stores all parameters.

# Parameters
 
## Physical Properties
* `wavelength` of the laser

## General Properties

* `dtype`: datatype used for reconstructions. `Float32` is usually faster, especially on GPUs.
* `propagator`: 

## Probe
* `dxp`: pixel size
* `Np`: number of pixels
* `xp`: 1D coordinates 
* `Lp`: field of view probe (field size)
* `zp`: distance to next plane

## Object 
* `dxo`: pixel size
* `No`: number of pixels
* `xo`: 1D coordinates 
* `Lo`: field of view probe (field size)
* `zo`: distance to next plane

## Detector
* `dxd`: pixel size
* `Nd`: number of pixels
* `xd`: 1D coordinates 
* `Ld`: field of view probe (field size)


"""
@kwdef mutable struct ExperimentalData{T}
    wavelength::T
    # probe sampling
    dxp::T
    Np::Int
    xp::Vector{T}
    Lp::T
    zp::Union{T, Nothing}
    entrancePupilDiameter::Union{T, Nothing}
    # object sampling
    dxo::T
    No::Int
    xo::Vector{T}
    Lo::T
    zo::T
    # detector sampling
    dxd::T
    Nd::Int
    xd::Vector{T}
    Ld::T
    # the general data type which is enforced
    dtype::Type{T}
    ptychogram::Union{Nothing, Array{T, N}} where N
    encoder::Union{Nothing, Array{T, 2}}
end

"""
    ExperimentalData(fileName::String)

Fill the `ExperimentalData` struct with data from a `*.hdf5` file.
"""
function ExperimentalData(fileName::String, mode=CPM::CPM)
    # open h5 file
    fid = HDF5.h5open(fileName)
   
    # function to extract number from a vector wrap
    r_number(x) = read(fid, x)[begin]
    # read arrays, take full
    r_array(x) = read(fid, x)
    # call a simple constructor and fill with data
    expData = initializeExperimentalData(
            # those numbers are stored as arrays, take first element
            wavelength=r_number("wavelength"),
            No=r_number("No"),
            Nd=r_number("Nd"),
            dxd=r_number("dxd"),
            zo=r_number("zo"),
            entrancePupilDiameter=r_number("encoder"),
            # those are arrays, keep them
            ptychogram=r_array("ptychogram"),
            encoder=r_array("encoder"),
        )

    return expData
end

"""
    initializeExperimentalData(<kwargs>)

Function to return parameter object `mutable struct` storing
all meta information needed for reconstruction.
Note that you can access all those parameters but many
of them are connected hence to fill the struct we only need a few of
them.
"""
function initializeExperimentalData(;
        wavelength=DEFAULT_WAVELENGTH, 
        No=2^7,
        Nd=2^9,
        dxd=4.5e-6,
        zo=50e-3,
        dtype=Float32,
        zp=nothing,
        entrancePupilDiameter=nothing,
        ptychogram=nothing,
        encoder=nothing,
        )

    # detector
    Ld = Nd * dxd
    xd = FourierTools.fftpos(Ld, Nd) 
    # probe 
    dxp = wavelength * zo / Ld       
    Np = Nd
    Lp = dxp*Np
    xp = FourierTools.fftpos(Lp, Np) 
    # object 
    dxo = dxp
    Lo = dxo * No
    xo = FourierTools.fftpos(Lo, No) 

    # anonymous function to convert to the dtype
    _dtype(x) = isnothing(x) ? x : dtype(x)
    

    # fill struct
    ExperimentalData(
        wavelength=_dtype(wavelength),
        # probe sampling
        dxp=_dtype(dxp),
        Np=Np,
        xp=_dtype.(xp),
        Lp=_dtype(Lp),
        zp=_dtype(zp),
        # object sampling,
        dxo=_dtype(dxo),
        No=No,
        xo=_dtype.(xo),
        Lo=_dtype(Lo),
        zo=_dtype(zo),
        # detector sampling,
        dxd=_dtype(dxd),
        Nd=Nd,
        entrancePupilDiameter=dtype(entrancePupilDiameter),
        xd=_dtype.(xd),
        Ld=_dtype(Ld),
        # the general data type which is enforced,
        dtype=dtype,
        ptychogram=dtype.(ptychogram)),
        encoder=dtype.(encoder)
end

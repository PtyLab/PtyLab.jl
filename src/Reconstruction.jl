export ReconstructionCPM


 # maybe move to a different file
export InitialObject, InitialObjectOnes
export InitialProbe, InitialProbeCirc

abstract type InitialObject end
abstract type InitialObjectOnes <:InitialObject end

abstract type InitialProbe end
abstract type InitialProbeCirc <: InitialProbe end


@kwdef mutable struct ReconstructionCPM{T, N}
    data::ExperimentalDataCPM{T}
    wavelength::T
    zo::T
    xdx::T
    # reconstruction parameters
    nlambda::Int
    nosm::Int
    npsm::Int
    nslice::Int
    # optional parameters
    entrancePupilDiameter::Union{Nothing, T}
    spectralDensity::Union{Nothing, T}
    theta::Union{Nothing, T}
    # reconstructions
    object::Union{Nothing, Array{T, 2}}
    probe::Union{Nothing, Array{T, N} where N}
end


"""
    ReconstructionCPM(data::ExperimentalDataCPM)

Fill the `ReconstructionCPM` struct with `experimentData` and some 
initial guesses.
"""
function ReconstructionCPM(data::ExperimentalDataCPM)
    nlambda = 1
    ambda = 1
    nosm = 1
    npsm = 1
    nslice = 1

    # beam and object purity
    purityProbe = 1
    purityObject = 1

    dxp = data.wavelength * data.zo / data.Ld

    # if entrancePupilDiameter is not provided in the hdf5 file, set it to be one third of the probe FoV.
    entrancePupilDiameter = 
        let 
            if isnothing(data.entrancePupilDiameter)
                calc_Lp(calc_Np(Nd), dxp)
            else 
                data.entrancePupilDiameter
            end
        end
            
    # wrap in array if not provided
    spectralDensity = isnothing(data.spectralDensity) ? [data.wavelength] : data.spectralDensity

    self.initialObject = InitialObjectOnes
    self.initialProbe = InitialProbeCirc




end

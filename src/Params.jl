export Params


abstract type IntensityConstraint end
abstract type IntensityConstraintStandard <: IntensityConstraint end


"""
Some settings are shared in between different optimizers, such as the type of propagatorType that you intend to use,
if you want to use probe orthogonalization, etc. These are stored in the reconstruction_parameters object.

This ensures that code like this will work as expected:
"""
@with_kw mutable struct Params
    # Default settings for switches, settings that involve how things are computed
    fftshiftSwitch::Bool = false 
    fftshiftFlag::Bool = false
    FourierMaskSwitch::Bool = false
    CPSCswitch::Bool = false
    CPSCupsamplingFactor::Bool = false
    intensityConstraint::Type{<:IntensityConstraint} = IntensityConstraintStandard 
    propagatorType::Function = Fraunhofer
    momentumAcceleration::Bool = false
    adaptiveMomentumAcceleration::Bool = false
    positionOrder::Type{<:PositionOrder} = GridRegularRand
end

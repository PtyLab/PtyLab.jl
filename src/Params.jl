export Params


abstract type IntensityConstraint end
abstract type IntensityConstraintStandard <: IntensityConstraint end


"""
@with_kw mutable struct Params

Those are some settings which are shared in between different optimizers, such as the type of propagatorType that you intend to use.

 # TODO

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

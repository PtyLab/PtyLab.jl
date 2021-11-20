export Params


abstract type IntensityConstraint end
struct IntensityConstraintStandard <:IntensityConstraint

end


"""
@with_kw mutable struct Params

Those are some settings which are shared in between different optimizers, such as the type of propagatorType that you intend to use.

 # Fields

* `fftshiftFlag::Bool = false`: `false` means that the ptychogram is already in the center 
* `propagatorType = Fraunhofer`: Default is `Fraunhofer`. See other Operator options. 
* `transposePtychogram::Bool = true`: swap the first two dimensions of the ptychogram (transpose the sensor)
* `randPositionOrder::Bool = true` randomly draw the encoder positions during reconstruction
* `comStabilizationSwitch::Bool = true` center of mass -> forces probe to be in the center
"""
@with_kw mutable struct Params
    # Default settings for switches, settings that involve how things are computed
    fftshiftFlag::Bool = false
    transposePtychogram::Bool = true
    intensityConstraint::Type{<:IntensityConstraint} = IntensityConstraintStandard 
    propagatorType::Function = Fraunhofer
    randPositionOrder::Bool = true 
    comStabilizationSwitch::Bool = true
end

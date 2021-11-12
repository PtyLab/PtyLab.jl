export Params


abstract type IntensityConstraint end
struct IntensityConstraintStandard <:IntensityConstraint

end


"""
@with_kw mutable struct Params

Those are some settings which are shared in between different optimizers, such as the type of propagatorType that you intend to use.

 # Fields

* `fftshiftFlag`: TODO
* `propagatorType`: Default is `Fraunhofer`. See other Operator options. 
"""
@with_kw mutable struct Params
    # Default settings for switches, settings that involve how things are computed
    fftshiftFlag::Bool = false
    propagatorType::Function = Fraunhofer
end

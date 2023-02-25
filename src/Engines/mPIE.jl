"""
    @with_kw struct mPIE{T} <: Engines

`mPIE` is a struct to store important parameters for the `mPIE::Engine`.
In most cases you want to initialize the numbers as `Float32` for higher performance.
Needs to be the same datatype as your object you want to reconstruct.


# Parameters
* `numIterations::Int=50`
* `betaProbe=0.25f0`
* `betaObject=0.25f0`
* `alphaProbe=0.1f0`
* `alphaObject=0.1f0`
* `feedbackM=0.3f0`
* `frictionM=0.7f0`
"""
@with_kw struct mPIE{T} <: Engines
    numIterations::Int = 50
    betaProbe::T = 0.25f0
    betaObject::T = 0.25f0
    # probe regularization
    alphaProbe::T = 0.1f0 
    # object regularization
    alphaObject::T = 0.1f0
    feedbackM::T = 0.3f0
    frictionM::T = 0.7f0
end



"""
    createUpdateFunctions(engine::mPIE{T}, objectPatch, probe, DELTA) where T

Return the functions `probeUpdate` and `objectPatchUpdate` which are used
to update the `probe` and `objectPatch`.

We create a lot of memory buffers.
Note, for the `sum` operations, the dimensions over summation is already
specified by the buffer.
"""
function createUpdateFunctions(engine::mPIE{T}, objectPatch, probe, object, DELTA) where T
    # the let statement is used such that the variables are available in the
    # functions as closures. The actualy operation is not meaningful but
    # instead we want to be sure to have the right shape 
    probeUpdate = let   fracProbe = similar(objectPatch), 
                        newProbe = similar(probe), 
                        fracProbeDELTA = fracProbe .* DELTA,
                        sumBufferNewProbe = sum(fracProbeDELTA, dims=(5,)),
                        abs2objectPatch = similar(objectPatch, real(eltype(objectPatch))),
                        sumabs2objectPatch = similar(abs2objectPatch, size(abs2objectPatch, 1), size(abs2objectPatch, 2),
                                                                      1, 1, 1, 1),
                        engine = engine
        function probeUpdate(objectPatch, probe, DELTA)
            abs2objectPatch .= abs2.(objectPatch)
            Omax = maximum(sum(abs2objectPatch, dims=(3, 4, 5, 6)), dims=(1, 2))
            fracProbe .= conj.(objectPatch) ./ maximum(sum!(sumabs2objectPatch, abs2objectPatch))
            fracProbe .= conj.(objectPatch) ./ 
                            (engine.alphaProbe .* Omax .+ (1 .- engine.alphaProbe) .* abs2objectPatch)
            fracProbeDELTA .= fracProbe .* DELTA
            newProbe .= probe .+ engine.betaProbe .* sum!(sumBufferNewProbe, fracProbeDELTA)

            return newProbe
        end
    end
   
    objectPatchUpdate = let engine = engine,
                            fracObject = similar(probe),
                            newObject = similar(objectPatch),
                            fracObjectDELTA = fracObject .* DELTA,
                            sumBufferNewObject = sum(fracObjectDELTA, dims=(4,)),
                            abs2probe = similar(probe, real(eltype(probe))),
                            sumabs2probe = similar(abs2probe, size(abs2probe, 1), 
                                                   size(abs2probe, 2),
                                                   1, 1, 1, 1)
        function objectPatchUpdate(objectPatch, probe, DELTA) 
            abs2probe .= abs2.(probe)
            Pmax = maximum(sum(abs2probe, dims=(3,4,5,6)), dims=(1,2))
            fracObject .= conj.(probe) ./ (engine.alphaObject .* Pmax .+ (1 .- engine.alphaObject) .* abs2probe)
            fracObjectDELTA .= fracObject .* DELTA
            newObject .= objectPatch .+ engine.betaObject .* sum!(sumBufferNewObject, fracObjectDELTA)
            return newObject
        end
    end


    objectMomentumUpdate! = let objectBuffer = similar(object),
                                engine = engine
        function objectMomentumUpdate!(object)
            gradient = reconstruction.objectBuffer - object
            objectMomentum = (gradient .+ self.frictionM .* objectMomentum)
            object .= (object .- engine.feedbackM .* self.reconstruction.objectMomentum)
            objectBuffer .= object 
        end
    end

    probeMomentumUpdate! = let probeBuffer = similar(probe)
                               engine = engine
        
        function probeMomentum!(probe)
            gradient = probeBuffer .- reconstruction.probe
            probeMomentum = gradient .+ engine.frictionM .* self.reconstruction.probeMomentum
            probe .= probe .- engine.feedbackM .* engine.probeMomentum
            probeBuffer .= probe
        end
    end

    return (; objectPatchUpdate, probeUpdate, objectMomentumUpdate!, probeMomentumUpdate!) 
end

export ePIE


"""
    @with_kw struct ePIE{T} <: Engines

`ePIE` is a struct to store important parameters for the `ePIE::Engine`.
In most cases you want to initialize the numbers as `Float32` for higher performance.
Needs to be the same datatype as your object you want to reconstruct.


# Parameters
* `betaProbe=0.25f0`
* `betaObject=0.25f0`
* `numIterations::Int=50`
"""
@with_kw struct ePIE{T} <: Engines
    betaProbe::T=0.25f0
    betaObject::T=0.25f0
    numIterations::Int=50
end




"""
    createUpdateFunctions(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Return the functions `probeUpdate` and `objectPatchUpdate` which are used
to update the `probe` and `objectPatch`.

We create a lot of memory buffers.
Note, for the `sum` operations, the dimensions over summation is already
specified by the buffer.
"""
function createUpdateFunctions(engine::Union{ePIE{T}}, objectPatch, probe, object, DELTA) where T
    # the let statement is used such that the variables are available in the
    # functions as closures. The actualy operation is not meaningful but
    # instead we want to be sure to have the right shape 
    probeUpdate = let   fracProbe = similar(objectPatch), 
                        engine = engine,
                        newProbe = similar(probe), 
                        fracProbeDELTA = fracProbe .* DELTA,
                        sumBufferNewProbe = sum(fracProbeDELTA, dims=(3,4,6)),
                        abs2objectPatch = similar(objectPatch, real(eltype(objectPatch))),
                        sumabs2objectPatch = similar(abs2objectPatch, size(abs2objectPatch, 1), size(abs2objectPatch, 2),
                                                                      1, 1, 1, 1)
        function probeUpdate(objectPatch, probe, DELTA)
            abs2objectPatch .= abs2.(objectPatch)
            fracProbe .= conj.(objectPatch) ./ maximum(sum!(sumabs2objectPatch, abs2objectPatch))
            fracProbeDELTA .= fracProbe .* DELTA
            newProbe .= probe .+ engine.betaProbe .* sum!(sumBufferNewProbe, fracProbeDELTA)
            return newProbe
        end
    end
   
    objectPatchUpdate = let fracObject = similar(probe),
                            engine = engine,
                            newObject = similar(objectPatch),
                            fracObjectDELTA = fracObject .* DELTA,
                            sumBufferNewObject = sum(fracObjectDELTA, dims=(3,5,6)),
                            abs2probe = similar(probe, real(eltype(probe))),
                            sumabs2probe = similar(abs2probe, size(abs2probe, 1), 
                                                   size(abs2probe, 2),
                                                   1, 1, 1, 1)
        function objectPatchUpdate(objectPatch, probe, DELTA) 
            abs2probe .= abs2.(probe)
            fracObject .= conj.(probe) ./ maximum(sum!(sumabs2probe, abs2probe))
            fracObjectDELTA .= fracObject .* DELTA
            newObject .= objectPatch .+ engine.betaObject .* sum!(sumBufferNewObject, fracObjectDELTA)
            return newObject
        end
    end

    return (; objectPatchUpdate, probeUpdate)
end

"""
    _loopUpdatePIE!(<args>)

Internal function which serves as the inner for loop of the PIE based algorithms
like `ePIE` or `zPIE`.
    
"""
function _loopUpdatePIE!(randPositionOrder, positions, Np, ptychogram, 
                         object, oldObjectPatch, oldProbe, probe, esw, DELTA,
                         intensityProjection!, probeUpdate, objectPatchUpdate)


    # random order of positions or not
    posList = randPositionOrder ? randperm(size(positions, 1)) : (1:size(positions, 1))

    for positionIndex in posList 
        row, col = view(positions, positionIndex, :)
            
        sy = row:(row + Np - 1)
        sx = col:(col + Np - 1)
        # no copy necessary -> because esw immediately calculated
        objectPatch = view(object, sy, sx, ..)
        
        # save old state since we need that later for probeUpdate and objectPatchUpdate
        oldObjectPatch .= objectPatch
        oldProbe .= probe
        
        # exit surface wave,
        esw .= objectPatch .* probe
        
        # already store esw in DELTA, since intensityProjection is going to change esw
        DELTA .= -1 .* esw
        eswUpdate = intensityProjection!(esw, view(ptychogram, :, :, positionIndex))
        
         # difference term
        DELTA .+= eswUpdate
        
        # update newProbe und newObjectPatch
        newObjectPatch = objectPatchUpdate(oldObjectPatch, oldProbe, DELTA)
        newProbe = probeUpdate(oldObjectPatch, oldProbe, DELTA)
        probe .= newProbe
        object[sy, sx, ..] .= newObjectPatch
    end 

    return nothing
end

"""
    _prepareBuffersAndFunctionsPIE(rec, params, engine)

Prepate buffers and functions needed for PIE based reconstructions
"""
function _prepareBuffersAndFunctionsPIE(rec, params, engine::Union{ePIE, zPIE})
    # create intensityProjection function 
    intensityProjection! = IntensityProjection(rec, params)

    # get two function which maybe shift depending on the flag             \
    maybe_fftshift, maybe_ifftshift = getMaybeFftshifts(! params.fftshiftFlag) 

    # ptychogram = rec.ptychogram
    ptychogram = maybe_ifftshift(rec.ptychogram)

    # if that flag is true, swap first and second dimension (transpose)
    if params.transposePtychogram
        ptychogram = permutedims(ptychogram, (2,1,3))
    end

    # those are buffers
    oldProbe = copy(rec.probe)
    oldObjectPatch = rec.object[1:rec.Np, 1:rec.Np, ..]
    # just a buffer with correct shape
    esw = rec.object[1:rec.Np, 1:rec.Np, ..] .* rec.probe
    DELTA = similar(esw)
    # performant update functions for probe and patch without memory usage.
    objectPatchUpdate, probeUpdate = createUpdateFunctions(engine, oldObjectPatch, oldProbe, DELTA)

    return (intensityProjection!, objectPatchUpdate, probeUpdate, 
            ptychogram, oldProbe, oldObjectPatch, esw, DELTA)
end



"""
    createUpdateFunctions(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Return the functions `probeUpdate` and `objectPatchUpdate` which are used
to update the `probe` and `objectPatch`.

We create a lot of memory buffers.
Note, for the `sum` operations, the dimensions over summation is already
specified by the buffer.
"""
function createUpdateFunctions(engine::Union{ePIE{T}, zPIE{T}}, objectPatch, probe, DELTA) where T
    # the let statement is used such that the variables are available in the
    # functions as closures. The actualy operation is not meaningful but
    # instead we want to be sure to have the right shape 
    probeUpdate = let   fracProbe = similar(objectPatch), 
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

    return objectPatchUpdate, probeUpdate
end


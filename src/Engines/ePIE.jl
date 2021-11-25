export ePIE

"""
    @with_kw mutable struct ePIE{T} <: Engines

`ePIE` is a struct to store important parameters for the `ePIE::Engine` 
"""
@with_kw mutable struct ePIE{T} <: Engines
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
function createUpdateFunctions(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracProbe = similar(objectPatch)
    newProbe = similar(probe)
    fracProbeDELTA = fracProbe .* DELTA
    sumBufferNewProbe = fracProbeDELTA[:, :, 1, :, 1, 1] 
    abs2objectPatch = real.(objectPatch)
    sumabs2objectPatch = real(similar(objectPatch, size(objectPatch, 1), size(objectPatch, 2)))


    fracObject = similar(probe)
    newObject = similar(objectPatch)
    fracObjectDELTA = fracObject .* DELTA
    sumBufferNewObject = fracObjectDELTA[:, :, 1, 1, :, 1] 
    abs2probe = real.(probe)
    sumabs2probe = real(similar(probe, size(probe, 1), size(probe, 2)))

    probeUpdate = let   fracProbe = fracProbe
                        newProbe = newProbe
                        fracProbeDELTA = fracProbeDELTA 
                        sumBufferNewProbe = sumBufferNewProbe
                        abs2objectPatch = abs2objectPatch
                        sumabs2objectPatch = sumabs2objectPatch 
        function probeUpdate(objectPatch, probe, DELTA) where T
            abs2objectPatch .= abs2.(objectPatch)
            fracProbe .= conj.(objectPatch) ./ maximum(sum!(sumabs2objectPatch, abs2objectPatch))
            fracProbeDELTA .= fracProbe .* DELTA
            newProbe .= probe .+ engine.betaProbe .* sum!(sumBufferNewProbe, fracProbeDELTA)
            return newProbe
        end
    end
   
    objectPatchUpdate = let   fracObject = fracObject
                        newObject = newObject
                        fracObjectDELTA = fracObjectDELTA 
                        sumBufferNewObject = sumBufferNewObject
                        abs2probe = abs2probe
                        sumabs2probe = sumabs2probe
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



"""
    reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 

Reconstruct a CPM dataset.
"""
function reconstruct(engine::ePIE{T}, params::Params, rec::ReconstructionCPM{T}) where T 
    # calculate the positions since rec.positions is in fact a function!
    positions = rec.positions

    # create intensityProjection function 
    intensityProjection! = IntensityProjection(rec, params)

    # get two function which maybe shift depending on the flag             \
    maybe_fftshift, maybe_ifftshift = getMaybeFftshifts(! params.fftshiftFlag) 

    # alias and shift maybe
    object = rec.object
    probe = rec.probe
    # ptychogram = rec.ptychogram
    ptychogram = maybe_ifftshift(rec.ptychogram)

    # if that flag is true, swap first and second dimension (transpose)
    if params.transposePtychogram
        ptychogram = collect(PermutedDimsArray(ptychogram, (2,1,3)))
    end
    Np = rec.Np

    # those are buffers
    oldProbe = copy(probe)
    oldObjectPatch = object[1:Np, 1:Np, ..]
    # just a buffer with correct shape
    esw = object[1:Np, 1:Np, ..] .* probe
    DELTA = similar(esw)
   
    # performant update functions for probe and patch without memory usage.
    objectPatchUpdate, probeUpdate = createUpdateFunctions(engine, oldObjectPatch, oldProbe, DELTA)

    # loop for iterations
    @showprogress for loop in 1:engine.numIterations
    # random order of positions or not
        posList = params.randPositionOrder ? randperm(size(positions, 2)) : (1:size(p, 2))
        for positionIndex in posList 
            # row, col does not work!
            # col, row = positions[:, positionIndex] 
            row, col = positions[:, positionIndex] 
                
            sy = row:(row + Np - 1)
            sx = col:(col + Np - 1)
            # @show col, row
            # no copy necessary -> because esw immediately calculated
            objectPatch = view(object, sy, sx, ..)

            # save old state since we need that later for probeUpdate and objectPatchUpdate
            oldObjectPatch .= objectPatch
            oldProbe .= probe

            # exit surface wave,
            # tullio seems to be slower on such simple operations
            #@tullio esw[i1,i2,i3,i4,i5,i6] = objectPatch[i1,i2,i3,i4,1,i6] .* probe[i1,i2,i3,1,i5,i6]
            esw .= objectPatch .* probe

            # already store esw in DELTA, since intensityProjection is going to change esw
            DELTA .= -1 .* esw
            eswUpdate = intensityProjection!(esw, view(ptychogram, :, :, positionIndex))
            
             # difference term
            DELTA .+= eswUpdate

            # update newProbe und newObjectPatch
            #newProbe, fracProbe = probeUpdate(engine, oldObjectPatch, oldProbe, DELTA, fracProbe, newProbe) 
            #newObjectPatch, fracObject = objectPatchUpdate(engine, oldObjectPatch, oldProbe, DELTA, 
			#				   fracObject, newObjectPatch) 
            #
            newProbe = probeUpdate(oldObjectPatch, oldProbe, DELTA)
            newObjectPatch = objectPatchUpdate(oldObjectPatch, oldProbe, DELTA)
            probe .= newProbe
            object[sy, sx, ..] .= newObjectPatch
        end 

        
        enforceConstraints!(rec, params)

    end
    return probe, object
end

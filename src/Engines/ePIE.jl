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
    # the let statement is used such that the variables are available in the
    # functions as closures. The actualy operation is not meaningful but
    # instead we want to be sure to have the right shape 
    probeUpdate = let   fracProbe = similar(objectPatch) 
                        newProbe = similar(probe) 
                        fracProbeDELTA = fracProbe .* DELTA
                        sumBufferNewProbe = sum(fracProbeDELTA, dims=(3,4,6)) 
                        abs2objectPatch = real.(objectPatch)
                        sumabs2objectPatch = sum(abs2objectPatch, dims=(3,4,5,6))
        function probeUpdate(objectPatch, probe, DELTA) where T
            abs2objectPatch .= abs2.(objectPatch)
            fracProbe .= conj.(objectPatch) ./ maximum(sum!(sumabs2objectPatch, abs2objectPatch))
            fracProbeDELTA .= fracProbe .* DELTA
            #@show size(newProbe), size(probe), size(sumBufferNewProbe), size(fracProbeDELTA)
            newProbe .= probe .+ engine.betaProbe .* sum!(sumBufferNewProbe, fracProbeDELTA)
            return newProbe
        end
    end
   
    objectPatchUpdate = let   fracObject = similar(probe)
        newObject = similar(objectPatch) 
                        fracObjectDELTA = fracObject .* DELTA
                        sumBufferNewObject = sum(fracObjectDELTA, dims=(3,5,6))
                        abs2probe = real.(probe)
                        sumabs2probe = sum(abs2probe, dims=(3,4,5,6))
        function objectPatchUpdate(objectPatch, probe, DELTA) 
            abs2probe .= abs2.(probe)
            fracObject .= conj.(probe) ./ maximum(sum!(sumabs2probe, abs2probe))
            fracObjectDELTA .= fracObject .* DELTA
            #@show size(newObject), size(objectPatch), size(sumBufferNewObject), size(fracObjectDELTA)
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
        ptychogram = permutedims(ptychogram, (2,1,3))
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

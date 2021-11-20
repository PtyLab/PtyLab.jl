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
    probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Returns the `newProbe`.

TODO memory improvements possible

"""
function probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracProbe = conj.(objectPatch) ./ maximum(sum(abs2.(objectPatch), dims=3:ndims(objectPatch)))
    newProbe = probe .+ engine.betaProbe .* sum(fracProbe .* DELTA, dims=(3, 5, 6))

    return newProbe
end

"""
    probeUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T

Returns the `newobjectPatch`.
TODO memory improvements possible
"""
function objectPatchUpdate(engine::ePIE{T}, objectPatch, probe, DELTA) where T
    fracObject = conj.(probe) ./ maximum(sum(abs2.(probe), dims=3:ndims(probe)))
    newObject = objectPatch .+ engine.betaObject .* sum(fracObject .* DELTA, dims=(3, 4, 6))

    return newObject
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
        ptychogram = PermutedDimsArray(ptychogram, (2,1,3))
    end
    Np = rec.Np

    # copy them!
    # those are buffers
    oldProbe = copy(probe)
    oldObjectPatch = object[1:Np, 1:Np, ..]
    esw = object[1:Np, 1:Np, ..] .* probe
    
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
            # @tullio esw[i1,i2,i3,i4,i5,i6] = objectPatch[i1,i2,i3,i4,1,i6] .* probe[i1,i2,i3,1,i5,i6]
            esw = objectPatch .* probe

            # already store esw in DELTA, since intensityProjection is going to change esw
            DELTA = -1 .* esw
            eswUpdate = intensityProjection!(esw, view(ptychogram, :, :, positionIndex))
            
             # difference term
            DELTA .+= eswUpdate

            # update newProbe und newObjectPatch
            newProbe = probeUpdate(engine, oldObjectPatch, oldProbe, DELTA) 
            newObjectPatch = objectPatchUpdate(engine, oldObjectPatch, oldProbe, DELTA) 
            probe .= newProbe
            object[sy, sx, ..] .= newObjectPatch
        end 

        
        enforceConstraints!(rec, params)

    end
    return probe, object
end

 # some helper functions
"""
    calc_Nd(ptychogram)

How many pixels per dimension.
"""
calc_Nd(ptychogram) = size(ptychogram, 1)
        

"""
    calc_xd(Nd, dxd)

Calculate detector coordinates in 1D.
"""
calc_xd(Nd, dxd) = typeof(dxd).(dxd .* range(-Nd / 2, Nd / 2, length=Nd))

"""
    calc_xo(No, dxo)

Calculate detector coordinates in 1D.
"""
calc_xo(No, dxo) = typeof(dxo).(dxo .* range(-No / 2, No / 2, length=No))


"""
    calc_xp(Np, dxp)

Calculate probe coordinates in 1D.
"""
calc_xp(Np, dxp) = typeof(dxp).(dxp .* range(-Np / 2, Np / 2, length=Np))



"""
    calc_Ld(Nd, dxd)

Calculate the size of the detector.
"""
calc_Ld(Nd, dxd) = Nd * dxd


"""
    calc_numFrames(ptychogram)

Calculate the number of frames.
"""
calc_numFrames(ptychogram) = size(ptychogram)[end]


"""
    calc_energyAtPos(ptychogram)

Calculate the energy at each position.
"""
calc_energyAtPos(ptychogram) = sum(abs, ptychogram, dims=(1,2))[1, 1, ..]


"""
    calc_maxProbePower(ptychogram)

Calculate max probe power at each position.
"""
calc_maxProbePower(ptychogram) = sqrt(maximum(sum(ptychogram, dims=(1,2))))


"""
    calc_Lp(Np, dxp)

Calculate probe length
"""
calc_Lp(Np, dxp) = Np * dxp



"""
    calc_Np(Nd)

Calculate probe pixel
"""
calc_Np(Nd) = Nd



"""
    calc_No(Nd) 

Calculate object pixel number
"""
calc_No(Np) = Np * 4

"""
    calc_NAd(Ld, zo) 

"""
calc_NAd(Ld, zo) = Ld / (2 * zo)



"""
    calc_DoF(lambda, NAd) 

"""
calc_DoF(lambda, NAd) = wavelength / NAd^2

"""
    calc_dxo(dxp) 

"""
calc_dxo(dxp) = dxp 


"""
    calc_Lo(No, dxo) 

"""
calc_Lo(No, dxo) = No * dxo 



"""
    calc_positions(::CPM, encoder, dxo, No, Np)

Calculate the positions in pixels from the `encoder` in meter.
"""
function calc_positions(::Type{<:CPM}, encoder, dxo, No, Np) 
    positions = round.(Int, encoder ./ dxo) .+ No ÷ 2 .- Np ÷ 2
    return positions
end


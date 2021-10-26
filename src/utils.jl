function circ(shape, xi, radius)
    rr2 = IndexFunArrays.rr2(shape, scale=ScaNorm, offset=CtrMid) .* 4 .* xi[end]^2
    return rr2 .< radius^2
end

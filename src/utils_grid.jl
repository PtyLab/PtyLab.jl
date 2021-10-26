export grid_regular_rand

abstract type PositionOrder end 

"""
    Tile(i₁, i₂, j₁, j₂)

Defines a tile starting at `(i₁, j₁)` and ending at `(i₂, j₂)`.
"""
struct Tile 
    i₁::Int
    i₂::Int
    j₁::Int
    j₂::Int
end


"""
    GridRegularRand(grid_size, tile_size, points, overlap)

Defines a regular grid which has random offsets.
"""
struct GridRegularRand
    grid_size::Tuple{Int, Int}
    tile_size::Tuple{Int, Int}
    tiles::Vector{Tile}
    overlap::Tuple{Float64, Float64}
end

"""
    grid_regular_rand(grid_size, tile_size, (N, M), rand_offset=10)

Returns a regular grid with `N` points per first dimension and `M` per second.
With random offsets (from the set `-rand_offset:rand_offset`) to each of the coordinates.

"""
function grid_regular_rand(grid_size, tile_size, (N, M), rand_offset=10)

    # regular positions which need to be randomly disturbed
    is = round.(Ref(Int), range(1, grid_size[1] - tile_size[1], length=N))
    js = round.(Ref(Int), range(1, grid_size[2] - tile_size[2], length=M))
   
    # linear overlap for each direction
    overlap = 1 .- ((grid_size[1] - 1)/(N-1), (grid_size[2] - 1)/(M-1)) ./ tile_size
    # struct storing all the output information
    grr = GridRegularRand(grid_size, tile_size, Vector{Tile}[], overlap)

    # wiggle the regular grids slightly
    for i in is
        for j in js
            # new random positions
            i_pos = clamp(random_constrained_wiggle(i, is, rand_offset), 1, grid_size[1] - tile_size[1])
            j_pos = clamp(random_constrained_wiggle(j, js, rand_offset), 1, grid_size[2] - tile_size[2])
            # save tile 
            push!(grr.tiles, Tile(i_pos, i_pos + tile_size[1] - 1, j_pos, j_pos + tile_size[2] - 1)) 
        end
    end
    return grr 
end


"""
    random_constrained_wiggle(k, ks, rand_offset)

We need this function to prevent that we wiggle out of the area to scan.
`ks` all indices which will be wiggled and hence the information about
the size of the scanning area.

`k` the current index to wiggle.
`rand_offset` wiggle offset

"""
function random_constrained_wiggle(k, ks, rand_offset)
    if k == ks[1]
        return k + rand(0:rand_offset)
    elseif k == ks[end]
        return k + rand(-rand_offset:0)
    else 
        return k + rand(-rand_offset:rand_offset)
    end
end

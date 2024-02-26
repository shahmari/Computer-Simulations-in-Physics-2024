module GameOfLife
export simulate_game
using StatsBase: sample, Weights
using Plots: heatmap, plot, cgrad, Animation, frame

function initialize_grid(rows::Int, cols::Int, dens::Float64)::Matrix{Bool}
    # grid = rand(Bool, rows, cols)
    grid = sample(Bool[true, false], Weights([dens, 1 - dens]), (rows, cols))
    return grid
end

function count_neighbors(grid::Matrix{Bool}, i::Int, j::Int)::Int
    count = 0
    Neighbors = [[1 + i % size(grid, 1), j], [i, 1 + j % size(grid, 2)],
        [(-2 + i + size(grid, 1)) % size(grid, 1) + 1, j],
        [i, (-2 + j + size(grid, 2)) % size(grid, 2) + 1],
        [1 + i % size(grid, 1), 1 + j % size(grid, 2)],
        [(-2 + i + size(grid, 1)) % size(grid, 1) + 1, (-2 + j + size(grid, 2)) % size(grid, 2) + 1],
        [1 + i % size(grid, 1), (-2 + j + size(grid, 2)) % size(grid, 2) + 1],
        [(-2 + i + size(grid, 1)) % size(grid, 1) + 1, 1 + j % size(grid, 2)]]
    @simd for neighbor ∈ Neighbors
        @inbounds count += grid[neighbor...] # modulo to handle edge cases
    end

    # Alternative implementation using array comprehension:
    # return sum([grid[neighbor...] for neighbor ∈ Neighbors])

    return count
end

function update_grid(grid::Matrix{Bool})::Matrix{Bool}
    new_grid = deepcopy(grid)
    @simd for i ∈ 1:size(grid, 1)
        @simd for j ∈ 1:size(grid, 2)
            neighbors = count_neighbors(grid, i, j)
            if @inbounds grid[i, j]
                if neighbors < 2 || neighbors > 3
                    @inbounds new_grid[i, j] = false
                end
            else
                if neighbors == 3
                    @inbounds new_grid[i, j] = true
                end
            end
        end
    end

    # Alternative implementation using CartesianIndices iterating and Tuple destructuring:

    # @simd for index ∈ CartesianIndices(grid)
    #     i, j = Tuple(index)
    #     neighbors = count_neighbors(grid, i, j)
    #     if grid[i, j]
    #         if neighbors < 2 || neighbors > 3
    #             new_grid[i, j] = false
    #         end
    #     else
    #         if neighbors == 3
    #             new_grid[i, j] = true
    #         end
    #     end
    # end

    return new_grid
end

@inline function simulate_game(rows::Int, cols::Int, steps::Int; dens::Float64=0.05)::Animation
    grid = initialize_grid(rows, cols, dens)
    anim = Animation()    
    for step ∈ 1:steps
        heatmap(grid', c=cgrad([:white, :green]), title="Step: $step", axis=nothing, grid=false, legend=false,frame=:none, ratio=1, size=(500, 500))
        frame(anim)
        grid = update_grid(grid)
    end
    return anim
end

end
cd(@__DIR__)

using Plots: gif

include("GameOfLife-Modules.jl")
using .GameOfLife

gif(simulate_game(100, 100, 500, dens=0.1), "game_of_life.gif", fps=10)
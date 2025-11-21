include("src/solve.jl")

using JSON
using .PdRans

params = JSON.parsefile("simple-setup.json")

PdRans.solve(params)
include("src/solve.jl")

using JSON
using .PdRans

params = JSON.parsefile("simple_setup.json")

PdRans.solve(params)
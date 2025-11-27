using JSON
using ArgParse
include("src/solve.jl")
using .PdRans

println(@__FILE__)

argset = ArgParseSettings()
@add_arg_table argset begin
    "--shutup"
        help = "do not display convergence history every step"
        action = :store_true
end
args = parse_args(argset)

params = JSON.parsefile("u-sampling-setup.json")

det_params = deepcopy(params)

function get_mean(var)
    if var["type"] == "normal"
        return var["mean"]
    elseif var["type"] == "uniform"
        range = var["range"]
        return 0.5*(range[1] + range[2])
    end
end
folder = "$(params["prefix"])-det"
mkpath(folder)
det_params["inlet u"] = get_mean(params["inlet u"])
det_params["output"] = "$folder/result.csv"

while true
    try
        PdRans.solve(det_params; verbose=!args["shutup"])
        break
    catch e
        @warn e stacktrace(catch_backtrace())
        sleep(1)
    end
end
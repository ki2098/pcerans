using JSON
using ArgParse
include("src/solve.jl")
using .PdRans

argset = ArgParseSettings()
@add_arg_table argset begin
    "--shutup"
        help = "do not display convergence history every step"
        action = :store_true
end
args = parse_args(argset)

params = JSON.parsefile("uti-sampling-setup.json")

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
det_params["inlet I"] = get_mean(params["inlet I"])
det_params["output"] = "$folder/result.csv"

while true
    try
        PdRans.solve(det_params; verbose=!args["shutup"])
        break
    catch e
        bt = catch_backtrace()
        s = sprint(Base.showerror, e, bt)
        @warn s
        sleep(1)
    end
end
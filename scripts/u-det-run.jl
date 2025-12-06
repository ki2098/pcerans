using JSON
using ArgParse
include("../src/solve.jl")
using .PdRans

println(@__FILE__)

argset = ArgParseSettings()
@add_arg_table argset begin
    "--skip-history"
        help = "do not display convergence history every step"
        action = :store_true
    "file"
        help = "setup file path"
        required = true
        arg_type = String
end
args = parse_args(argset)

params = JSON.parsefile(args["file"])

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
        PdRans.solve(det_params; show_history=!args["skip-history"])
        break
    catch e
        bt = catch_backtrace()
        s = sprint(Base.showerror, e, bt)
        @warn s
        sleep(1)
    end
end
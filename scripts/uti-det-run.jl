using JSON
using ArgParse

include("../src/solve-v2.jl")
using .PdRans

println(@__FILE__)

argset = ArgParseSettings()
@add_arg_table argset begin
    "--skip-history"
        help = "do not display convergence history every step"
        action = :store_true
    "case"
        help = "case path"
        required = true
        arg_type = String
end
args = parse_args(argset)

case_path = args["case"]
setup_path = "$case_path/setup.json"
params = JSON.parsefile(setup_path)
folder = "$case_path/data/det"
mkpath(folder)

function get_mean(var)
    if var["type"] == "normal"
        return var["mean"]
    elseif var["type"] == "uniform"
        range = var["range"]
        return 0.5*(range[1] + range[2])
    end
end
det_params = deepcopy(params)
det_params["inlet u"] = get_mean(params["inlet u"])
det_params["inlet I"] = get_mean(params["inlet I"])
det_params["output"] = "$folder/result.csv"

while true
    try
        logstr = PdRans.solve(det_params; show_history=!args["skip-history"])
        print("$logstr\n")
        break
    catch e
        bt = catch_backtrace()
        s = sprint(Base.showerror, e, bt)
        @warn s
        sleep(1)
    end
end
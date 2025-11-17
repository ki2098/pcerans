using JSON
include("src/solve.jl")
using .PdRans

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
folder = "data/uti-det"
mkpath(folder)
det_params["inlet u"] = get_mean(params["inlet u"])
det_params["inlet I"] = get_mean(params["inlet I"])
det_params["output"] = "$folder/result.csv"

PdRans.solve(det_params)
using JSON
using PolyChaos
using DataFrames
using CSV
using ArgParse
using Dates

argset = ArgParseSettings()
@add_arg_table argset begin
    "--samples"
        help = "generate samples"
        action = :store_true
    "--statistics"
        help = "calculate statistics"
        action = :store_true
end
args = parse_args(argset)
if !args["samples"] && !args["statistics"]
    args["samples"] = true
    args["statistics"] = true
end

params = JSON.parsefile("uti_sampling_setup.json")

det_params = deepcopy(params)

degree = 5
n_samples = degree*2
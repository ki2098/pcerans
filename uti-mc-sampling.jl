using JSON
using DataFrames
using CSV
using ArgParse
using Dates
using Distributions
using Random

argset = ArgParseSettings()
@add_arg_table argset begin
    "--samples"
        help = "generate samples"
        action = :store_true
    "--statistics"
        help = "calculate statistics"
        action = :store_true
    "n"
        help = "number of samples per dimension"
        required = true
        arg_type = Int
end
args = parse_args(argset)
if !args["samples"] && !args["statistics"]
    args["samples"] = true
    args["statistics"] = true
end

n_samples = args["n"]
total_samples = n_samples^2

params = JSON.parsefile("uti-sampling-setup.json")
folder = "data/uti-mc-$(total_samples)-samples"
mkpath(folder)

det_params = deepcopy(params)

function get_dist(var)
    if var["type"] == "normal"
        d = Normal(var["mean"], var["sd"])
    elseif var["type"] == "uniform"
        d = Uniform(var["range"]...)
    end
    println("dist = $d")
    x = rand(d, n_samples)
    return x
end

uin = get_dist(params["inlet u"])
tiin = get_dist(params["inlet I"])

if args["samples"]
    include("src/solve.jl")
    using .PdRans
    @time "$total_samples MC samples" for I=1:n_samples,J=1:n_samples
        sample_id = (I-1)*n_samples+J
        println("Monte Carlo sample $sample_id/$total_samples")
        det_params["output"] = "$folder/sample-$sample_id.csv"
        det_params["inlet u"] = uin[I]
        det_params["inlet I"] = tiin[J]
        while true
            try
                PdRans.solve(det_params)
                # println("Lucky!")
                break
            catch e
                @warn "get error=$e"
                sleep(1)
            end
        end
        println()
    end
end

if args["statistics"]
    sz = params["divisions"]
    cnt = sz[1]*sz[2]
    u_mean = zeros(cnt)
    u_var  = zeros(cnt)
    v_mean = zeros(cnt)
    v_var  = zeros(cnt)
    p_mean = zeros(cnt)
    p_var  = zeros(cnt)
    x = y = z = nothing
    for sample_id = 1:total_samples
        filename = "$folder/sample-$sample_id.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]
        du = u - u_mean
        dv = v - v_mean
        dp = p - p_mean
        u_mean .+= du ./ sample_id
        v_mean .+= dv ./ sample_id
        p_mean .+= dp ./ sample_id
        du2 = u - u_mean
        dv2 = v - v_mean
        dp2 = p - p_mean
        u_var .+= du .* du2
        v_var .+= dv .* dv2
        p_var .+= dp .* dp2
        if sample_id == 1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end
        print("\rread $sample_id/$total_samples")
    end
    u_var ./= total_samples
    v_var ./= total_samples
    p_var ./= total_samples
    println()

    stat_df = DataFrame(
        "x"=>x, "y"=>y, "z"=>z,
        "E[u]"=>u_mean, "Var[u]"=>u_var,
        "E[v]"=>v_mean, "Var[v]"=>v_var,
        "E[p]"=>p_mean, "Var[p]"=>p_var
    )

    CSV.write("$folder/statistics.csv", stat_df)
end
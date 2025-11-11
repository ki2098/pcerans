using JSON
using DataFrames
using CSV
using ArgParse
using Dates
using Distributions

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

params = JSON.parsefile("u_sampling_setup.json")

det_params = deepcopy(params)

n_samples_few = 100
n_samples_many = n_samples_few*10

function get_dist(var)
    if var["type"] == "normal"
        d = Normal(var["mean"], var["sd"])
    elseif var["type"] == "uniform"
        d = Uniform(var["range"]...)
    end
    println("dist = $d")
    x = rand(d, n_samples_many)
    return x
end

uin = get_dist(params["inlet u"])
base_filename = "data/u_mc_sample"

if args["samples"]
    include("src/solve.jl")
    using .PdRans

    i_sample = 1
    @time "$n_samples_many MC samples" while i_sample <= n_samples_many
        println("Monte Carlo sample $i_sample/$n_samples_many")
        det_params["output"] = "$base_filename.$i_sample.csv"
        det_params["inlet u"] = uin[i_sample]
        try
            PdRans.solve(det_params)
        catch e
            println("err = %e, retry")
            continue
        end
        global i_sample += 1
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
    for i_sample = 1:n_samples_many
        filename = "$base_filename.$i_sample.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]
        du = u - u_mean
        dv = v - v_mean
        dp = p - p_mean
        u_mean .+= du ./ i_sample
        v_mean .+= dv ./ i_sample
        p_mean .+= dp ./ i_sample
        du2 = u - u_mean
        dv2 = v - v_mean
        dp2 = p - p_mean
        u_var .+= du .* du2
        v_var .+= dv .* dv2
        p_var .+= dp .* dp2
        if i_sample == 1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end
        print("\rread $filename")
        if i_sample == n_samples_few
            u_var_tmp = u_var/n_samples_few
            v_var_tmp = v_var/n_samples_few
            p_var_tmp = p_var/n_samples_few
            few_df = DataFrame(
                "x"=>x, "y"=>y, "z"=>z,
                "E[u]"=>u_mean, "Var[u]"=>u_var_tmp,
                "E[v]"=>v_mean, "Var[v]"=>v_var_tmp,
                "E[p]"=>p_mean, "Var[p]"=>p_var_tmp
            )
            CSV.write("$base_filename.statistics$n_samples_few.csv", few_df)
        end
    end
    u_var ./= n_samples_many
    v_var ./= n_samples_many
    p_var ./= n_samples_many
    println()

    many_df = DataFrame(
        "x"=>x, "y"=>y, "z"=>z,
        "E[u]"=>u_mean, "Var[u]"=>u_var,
        "E[v]"=>v_mean, "Var[v]"=>v_var,
        "E[p]"=>p_mean, "Var[p]"=>p_var
    )
    CSV.write("$base_filename.statistics$n_samples_many.csv", many_df)
end
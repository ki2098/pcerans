using JSON
using DataFrames
using CSV
using ArgParse
using Dates
using Distributions
using Base.Threads

println(@__FILE__)

argset = ArgParseSettings()
@add_arg_table argset begin
    "--samples"
        help = "generate samples"
        action = :store_true
    "--statistics"
        help = "calculate statistics"
        action = :store_true
    "--skip-history"
        help = "do not display convergence history every step"
        action = :store_true
    "case"
        help = "case path"
        required = true
        arg_type = String
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

samples_per_dim = args["n"]
n_samples = samples_per_dim^2

case_path = args["case"]
setup_path = "$case_path/setup.json"
params = JSON.parsefile(setup_path)
folder = "$case_path/data/mc-$(n_samples)-samples"

function get_dist(var, n)
    if var["type"] == "normal"
        d = Normal(var["mean"], var["sd"])
    elseif var["type"] == "uniform"
        d = Uniform(var["range"]...)
    end
    println("dist = $d")
    x = rand(d, n)
    return x
end

u_nodes = get_dist(params["inlet u"], samples_per_dim)
println("uin = $u_nodes")
I_nodes = get_dist(params["inlet I"], samples_per_dim)
println("Iin = $I_nodes")

uin = vec([u_nodes[i] for i=1:samples_per_dim,j=1:samples_per_dim])
Iin = vec([I_nodes[j] for i=1:samples_per_dim,j=1:samples_per_dim])

start_time = now()

if args["samples"]
    using CUDA
    include("../src/solve-v2.jl")
    using .PdRans

    mkpath(folder)

    @threads for i_sample = 1:n_samples
        rank = threadid()
        ngpu = length(CUDA.devices())
        gpuid = (rank)%ngpu
        CUDA.device!(gpuid)
        headstr = "job=$i_sample/$n_samples\ntid=$rank\ngpu=$(CUDA.device())\n"

        det_params = deepcopy(params)
        det_params["output"] = "$folder/sample-$i_sample.csv"
        det_params["inlet u"] = uin[i_sample]
        det_params["inlet I"] = Iin[i_sample]

        while true
            try
                solvestr = PdRans.solve(det_params; show_history=!args["skip-history"])
                print("$headstr$solvestr\n")
                break
            catch e
                bt = catch_backtrace()
                s = sprint(Base.showerror, e, bt)
                @warn "job=$i_sample\n$s"
                sleep(1)
            end
        end
    end
end

end_time = now()
elapse = Dates.value(end_time-start_time)/1000
println("$n_samples MC samples took $(elapse)s")

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
    for i_sample = 1:n_samples
        filename = "$folder/sample-$i_sample.csv"
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
        print("read $filename\n")
    end
    u_var ./= (n_samples - 1)
    v_var ./= (n_samples - 1)
    p_var ./= (n_samples - 1)
    println()

    stat_df = DataFrame(
        "x"=>x, "y"=>y, "z"=>z,
        "E[u]"=>u_mean, "Var[u]"=>u_var,
        "E[v]"=>v_mean, "Var[v]"=>v_var,
        "E[p]"=>p_mean, "Var[p]"=>p_var
    )

    CSV.write("$folder/statistics.csv", stat_df)
end
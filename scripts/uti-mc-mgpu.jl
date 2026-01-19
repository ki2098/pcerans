using JSON
using DataFrames
using CSV
using ArgParse
using Dates
using Distributions
using MPI

include("task-distributor.jl")

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

n_samples = args["n"]
total_samples = n_samples^2

case_path = args["case"]
setup_path = "$case_path/setup.json"
params = JSON.parsefile(setup_path)
folder = "$case_path/data/mc-$(total_samples)-samples"

det_params = deepcopy(params)

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)
root = 0

start_time = now()
if args["samples"]
    using CUDA
    devices = CUDA.devices()
    n_gpu = length(devices)
    i_gpu = rank % n_gpu
    CUDA.device!(i_gpu)
    println("rank $rank -> $(CUDA.device(CUDA.current_context()))")

    include("../src/solve.jl")
    using .PdRans

    if rank == root
        mkpath(folder)
        u_base = get_dist(params["inlet u"], n_samples)
        I_base = get_dist(params["inlet I"], n_samples)
        uin = vec([u_base[i] for i=1:n_samples,j=1:n_samples])
        Iin = vec([I_base[j] for i=1:n_samples,j=1:n_samples])
    else
        uin = nothing
        Iin = nothing
    end
    uin = MPI.bcast(uin, root, comm)
    Iin = MPI.bcast(Iin, root, comm)
    println("u samples = $uin")
    println("I samples = $Iin")

    jobs = split_view([i for i=1:total_samples], size)[rank+1]
    println("local jobs = $jobs")
    for sample_id in jobs
        println("Monte Carlo sample $sample_id/$total_samples")
        det_params["output"] = "$folder/sample-$sample_id.csv"
        det_params["inlet u"] = uin[sample_id]
        det_params["inlet I"] = Iin[sample_id]
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
        println()
    end
end
MPI.Barrier(comm)
end_time = now()
elapse = Dates.value(end_time-start_time)/1000
println("$n_samples MC samples took $(elapse)s")

if args["statistics"] && rank == root
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
    u_var ./= (total_samples - 1)
    v_var ./= (total_samples - 1)
    p_var ./= (total_samples - 1)
    println()

    stat_df = DataFrame(
        "x"=>x, "y"=>y, "z"=>z,
        "E[u]"=>u_mean, "Var[u]"=>u_var,
        "E[v]"=>v_mean, "Var[v]"=>v_var,
        "E[p]"=>p_mean, "Var[p]"=>p_var
    )

    CSV.write("$folder/statistics.csv", stat_df)
end

MPI.Finalize()
using JSON
using PolyChaos
using DataFrames
using CSV
using ArgParse
using Dates
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
    "d"
        help = "degree per dimension"
        required = true
        arg_type = Int
    "n"
        help = "number of nodes per dimension"
        required = true
        arg_type = Int
end
args = parse_args(argset)
if !args["samples"] && !args["statistics"]
    args["samples"] = true
    args["statistics"] = true
end

degree = args["d"]
n_samples = args["n"]

case_path = args["case"]
setup_path = "$case_path/setup.json"
params = JSON.parsefile(setup_path)
folder = "$case_path/data/pce-$(n_samples)-samples"

function get_op_and_transform(var, n)
    if var["type"] == "normal"
        op = GaussOrthoPoly(degree, Nrec=n+1, addQuadrature=true)
        return op, var["sd"], var["mean"]
    elseif var["type"] == "uniform"
        op = Uniform01OrthoPoly(degree, Nrec=n+1, addQuadrature=true)
        range = var["range"]
        return op, (range[2] - range[1]), range[1]
    end
end

op, scale, offset = get_op_and_transform(params["inlet u"], n_samples)
ξ = op.quad.nodes
println("ξ = $ξ")
uin = ξ .* scale .+ offset
println("uin = $(scale)ξ + $offset")
println("uin = $uin")

start_time = now()

if args["samples"]
    using CUDA
    include("../src/solve-v2.jl")
    using .PdRans

    mkpath(folder)
    # mylock = ReentrantLock()

    @threads for i_sample = 1:n_samples
        rank = threadid()
        ngpu = length(CUDA.devices())
        gpuid = (rank)%ngpu
        CUDA.device!(gpuid)
        headstr = "job=$i_sample/$n_samples\ntid=$rank\ngpu=$(CUDA.device())\n"

        det_params = deepcopy(params)
        det_params["output"] = "$folder/sample-$i_sample.csv"
        det_params["inlet u"] = uin[i_sample]

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
println("$n_samples PCE samples took $(elapse)s")

if args["statistics"]
    w = op.quad.weights
    τ = computeSP2(op)

    sz = params["divisions"]
    cnt = sz[1]*sz[2]
    U = zeros(cnt, degree+1)
    V = zeros(cnt, degree+1)
    P = zeros(cnt, degree+1)
    x = y = z = nothing

    for i_sample = 1:n_samples
        filename = "$folder/sample-$i_sample.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]
        ψ = evaluate(ξ[i_sample], op)
        for i=1:cnt
            for k=1:degree+1
                c = ψ[k]*w[i_sample]/τ[k]
                U[i, k] += u[i]*c
                V[i, k] += v[i]*c
                P[i, k] += p[i]*c
            end
        end
        if i_sample == 1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end
        print("read $filename\n")
    end
    println()

    function compute_statistics(v)
        s = zeros(cnt, 2)
        for i=1:cnt
            s[i, 1] = v[i, 1]
            for k=2:degree+1
                s[i, 2] += τ[k]*v[i, k]^2
            end
        end
        return s
    end

    u_stat = compute_statistics(U)
    v_stat = compute_statistics(V)
    p_stat = compute_statistics(P)
    stat_df = DataFrame(
        "x"=>x,
        "y"=>y,
        "z"=>z,
        "E[u]"=>u_stat[:,1],
        "Var[u]"=>u_stat[:,2],
        "E[v]"=>v_stat[:,1],
        "Var[v]"=>v_stat[:,2],
        "E[p]"=>p_stat[:,1],
        "Var[p]"=>p_stat[:,2]
    )

    CSV.write("$folder/statistics.csv", stat_df)
end
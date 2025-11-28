using JSON
using PolyChaos
using DataFrames
using CSV
using ArgParse
using Dates

println(@__FILE__)

argset = ArgParseSettings()
@add_arg_table argset begin
    "--samples"
        help = "generate samples"
        action = :store_true
    "--statistics"
        help = "calculate statistics"
        action = :store_true
    "--shutup"
        help = "do not display convergence history every step"
        action = :store_true
end
args = parse_args(argset)
if !args["samples"] && !args["statistics"]
    args["samples"] = true
    args["statistics"] = true
end

degree = 5
n_samples = degree*2

params = JSON.parsefile("u-sampling-setup.json")
folder = "$(params["prefix"])-nipce-$(n_samples)-samples"
mkpath(folder)

det_params = deepcopy(params)

function get_op_and_transform(var)
    if var["type"] == "normal"
        op = GaussOrthoPoly(degree, Nrec=n_samples+1, addQuadrature=true)
        return op, var["sd"], var["mean"]
    elseif var["type"] == "uniform"
        op = Uniform01OrthoPoly(degree, Nrec=n_samples+1, addQuadrature=true)
        range = var["range"]
        return op, (range[2] - range[1]), range[1]
    end
end

op, scale, offset = get_op_and_transform(params["inlet u"])
ξ = op.quad.nodes
println("ξ = $ξ")
println("uin = $(scale)ξ + $offset")

if args["samples"]
    include("src/solve.jl")
    using .PdRans

    uin = ξ .* scale .+ offset

    start_time = now()
    for i_sample=1:n_samples
        println("niPCE sample $i_sample/$n_samples")
        det_params["output"] = "$folder/sample-$i_sample.csv"
        det_params["inlet u"] = uin[i_sample]
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
        println()
    end
    end_time = now()
    elapse = Dates.value(end_time-start_time)/1000
    println("$n_samples niPCE samples took $(elapse)s")
end

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
        print("\rread $filename")
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

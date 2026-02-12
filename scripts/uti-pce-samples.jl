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
samples_per_dim = args["n"]
IJmap = vec([ [I,J] for I=1:samples_per_dim,J=1:samples_per_dim ])
n_samples = samples_per_dim^2

case_path = args["case"]
setup_path = "$case_path/setup.json"
params = JSON.parsefile(setup_path)
folder = "$case_path/data/pce-$(n_samples)-samples"

function get_op_and_transform(var, d, n)
    if var["type"] == "normal"
        op = GaussOrthoPoly(d, Nrec=n+1, addQuadrature=true)
        return op, var["sd"], var["mean"]
    elseif var["type"] == "uniform"
        op = Uniform01OrthoPoly(d, Nrec=n+1, addQuadrature=true)
        range = var["range"]
        return op, (range[2] - range[1]), range[1]
    end
end

opu, scaleu, offsetu = get_op_and_transform(params["inlet u"], degree, samples_per_dim)
ξu = opu.quad.nodes
opI, scaleI, offsetI = get_op_and_transform(params["inlet I"], degree, samples_per_dim)
ξI = opI.quad.nodes
println("ξu = $ξu")
println("u = $(scaleu) ξu + $offsetu")
println("ξI = $ξI")
println("I = $(scaleI) ξI + $offsetI")

uin = ξu .* scaleu .+ offsetu
Iin = ξI .* scaleI .+ offsetI
println("uin = $uin")
println("Iin = $Iin")

start_time = now()

if args["samples"]
    using CUDA
    include("../src/solve-v2.jl")
    using .PdRans

    mkpath(folder)

    @threads for i_sample = 1:n_samples
        I, J = IJmap[i_sample]
        rank = threadid()
        ngpu = length(CUDA.devices())
        gpuid = (rank)%ngpu
        CUDA.device!(gpuid)
        headstr = "job=$i_sample($I,$J)/$n_samples\ntid=$rank\ngpu=$(CUDA.device())\n"

        det_params = deepcopy(params)
        det_params["output"] = "$folder/sample-$i_sample.csv"
        det_params["inlet u"] = uin[I]
        det_params["inlet I"] = Iin[J]

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

    wu = opu.quad.weights
    τu = computeSP2(opu)
    wI = opI.quad.weights
    τI = computeSP2(opI)

    # w = [wu[I]*wI[J] for I=1:samples_per_dim,J=1:samples_per_dim]
    
    mop = MultiOrthoPoly([opu, opI], degree)
    N = mop.dim
    L = mop.ind

    mt2 = Tensor(2, mop)
    τ = [mt2.get([k, k]) for k=0:N-1]
    show(τ)

    U = zeros(cnt, N)
    V = zeros(cnt, N)
    P = zeros(cnt, N)
    x = y = z = nothing

    for i_sample = 1:n_samples
        I, J = IJmap[i_sample]
        filename = "$folder/sample-$i_sample.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]
        ψu = evaluate(ξu[I], opu)
        ψI = evaluate(ξI[J], opI)
        ψ = [ ψu[ L[k,1]+1 ]*ψI[ L[k,2]+1 ] for k=1:N ]
        w = wu[I]*wI[J]

        for k=1:N
            c = ψ[k]*w/τ[k]

            U[:, k] .+= u .* c
            V[:, k] .+= v .* c
            P[:, k] .+= p .* c
            
        end

        if I==1 && J==1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end

        print("\rread $filename ($I,$J)")
    end
    println()

    function compute_statistics(v)
        s = zeros(cnt, 2)
        s[:, 1] .= v[:, 1]
        for k = 2:N
            s[:, 2] .+= τ[k] .* v[:, k] .^ 2
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

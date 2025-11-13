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
base_filename = "data/uti_nipce_sample"

det_params = deepcopy(params)

degree = 5
n_samples = degree*2

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

opu, scaleu, offsetu = get_op_and_transform(params["inlet u"])
ξu = opu.quad.nodes
opI, scaleI, offsetI = get_op_and_transform(params["inlet I"])
ξI = opI.quad.nodes
println("ξu = $ξu")
println("u = $(scaleu) ξu + $offsetu")
println("ξI = $ξI")
println("I = $(scaleI) ξI + $offsetI")
total_samples = n_samples^2

if args["samples"]
    include("src/solve.jl")
    using .PdRans

    uin = ξu .* scaleu .+ offsetu
    tiin = ξI .* scaleI .+ offsetI
    @time "$(n_samples*n_samples) niPCE samples" for I=1:n_samples, J=1:n_samples
        sample_id = (I-1)*n_samples+J
        println("niPCE sample $sample_id/$total_samples")
        det_params["output"] = "$base_filename.$sample_id.csv"
        det_params["inlet u"] = uin[I]
        det_params["inlet I"] = tiin[J]
        while true
            try
                PdRans.solve(det_params)
                break
            catch e
                println("err = $e, retry")
            end
        end
        println()
    end
end

if args["statistics"]
    sz = params["divisions"]
    cnt = sz[1]*sz[2]

    wu = opu.quad.weights
    τu = computeSP2(opu)
    wI = opI.quad.weights
    τI = computeSP2(opI)

    opm = MultiOrthoPoly([opu, opI], degree)
    dm = opm.dim
    L = opm.ind .+ 1
    w = [wu[I] * wI[J] for I=1:n_samples, J=1:n_samples]
    τ = [τu[L[K,1]] * τI[L[K,2]] for K=1:dm]
    println(τ)

    U = zeros(cnt, dm)
    V = zeros(cnt, dm)
    P = zeros(cnt, dm)
    x = y = z = nothing
    for I=1:n_samples, J=1:n_samples
        sample_id = (I-1)*n_samples+J
        filename = "$base_filename.$sample_id.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]
        ψu = evaluate(ξu[I], opu)
        ψI = evaluate(ξI[J], opI)
        ψ = [ψu[L[K,1]] * ψI[L[K,2]] for K=1:dm]
        for K = 1:dm
            c = ψ[K]*w[I, J]/τ[K]
            for i = 1:cnt
                U[i, K] += u[i] * c
                V[i, K] += v[i] * c
                P[i, K] += p[i] * c
            end
        end
        if sample_id == 1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end
        print("\rread $filename")
    end
    println()

    function compute_statistics(v)
        s = zeros(cnt, 2)
        for i = 1:cnt
            s[i, 1] = v[i, 1]
            for K = 2:dm
                s[i, 2] += τ[K] * v[i, K]^2
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

    CSV.write("$base_filename.statistics.csv", stat_df)
end
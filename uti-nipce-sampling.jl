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
n_samples = 26
total_samples = n_samples^2

params = JSON.parsefile("uti-sampling-setup.json")
folder = "$(params["prefix"])-nipce-$(total_samples)-samples"
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

    start_time = now()
    for I=1:n_samples, J=1:n_samples
        sample_id = (I-1)*n_samples+J
        println("niPCE sample $sample_id/$total_samples")
        det_params["output"] = "$folder/sample-$sample_id.csv"
        det_params["inlet u"] = uin[I]
        det_params["inlet I"] = tiin[J]
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
    println("$total_samples niPCE samples took $(elapse)s")
end

if args["statistics"]
    sz = params["divisions"]
    cnt = sz[1]*sz[2]

    wu = opu.quad.weights
    τu = computeSP2(opu)
    wI = opI.quad.weights
    τI = computeSP2(opI)

    w = [wu[i]*wI[j] for i=1:n_samples,j=1:n_samples]
    τ = [τu[i]*τI[j] for i=1:degree+1,j=1:degree+1]

    cnt = sz[1]*sz[2]
    U = zeros(cnt, degree+1, degree+1)
    V = zeros(cnt, degree+1, degree+1)
    P = zeros(cnt, degree+1, degree+1)
    x = y = z = nothing
    for I=1:n_samples, J=1:n_samples
        sample_id = (I-1)*n_samples+J
        filename = "$folder/sample-$sample_id.csv"
        sample_df = CSV.read(filename, DataFrame)
        u = sample_df[!, "u"]
        v = sample_df[!, "v"]
        p = sample_df[!, "p"]

        for i=1:degree+1, j=1:degree+1
            ψu = evaluate(i-1, ξu[I], opu)
            ψI = evaluate(j-1, ξI[J], opI)
            ψ = ψu*ψI
            c = ψ*w[I,J]/τ[i,j]

            U[:, i, j] .+= u .* c
            V[:, i, j] .+= v .* c
            P[:, i, j] .+= p .* c
            
        end

        if I==1 && J==1
            global x = sample_df[!, "x"]
            global y = sample_df[!, "y"]
            global z = sample_df[!, "z"]
        end

        print("\rread $filename")
    end
    println()

    function compute_statistics(v)
        s = zeros(cnt, 2)
        s[:, 1] .= v[:, 1, 1]
        for i=1:degree+1, j=1:degree+1
            if !(i == 1 && j == 1)
                s[:, 2] .+= (v[:, i, j] .^ 2) .* τ[i, j] 
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
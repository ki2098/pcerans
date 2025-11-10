using JSON, CSV, DataFrames

include("eq.jl")
include("rans.jl")
include("bc.jl")
include("pd.jl")

nthread=(16,16)

struct Ls
    A
    b
    r
    ω
    maxdiag
    maxit
    maxerr
end

struct Solver
    u
    v
    ut
    vt
    uu
    vv
    p
    k
    ω
    kt
    ωt
    nut
    divU
    dfunc
    uin
    kin
    ωin
    ls::Ls
    x
    y
    dx
    dy
    dt
    Re
    sz
    gc
    maxstep
end

function init(path)
    gc = 2
    init_json = JSON.parsefile(path)
    domain_json = init_json["domain"]
    xmin = domain_json["x range"][1]
    xmax = domain_json["x range"][2]
    ymin = domain_json["y range"][1]
    ymax = domain_json["y range"][2]
    nx = domain_json["divisions"][1]
    ny = domain_json["divisions"][2]
    sz = (nx + 2*gc, ny + 2*gc)
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny
    println("DOMAIN INFO")
    println("\tgc = $gc")
    println("\tx range = [$xmin, $xmax]")
    println("\ty range = [$ymin, $ymax]")
    println("\tdivison = ($(sz[1]-2*gc) $(sz[2]-2*gc))")
    println("\tcell size = ($dx, $dy)")

    time_json = init_json["time"]
    maxtime = time_json["end"]
    dt = time_json["dt"]
    maxstep::Int = maxtime/dt
    println("TIME INFO")
    println("\tend time = $maxtime")
    println("\tdt = $dt")
    println("\ttotal steps = $maxstep")

    Re = init_json["cfd"]["Re"]
    println("CFD INFO")
    println("\tRe = $Re")

    inlet_json = init_json["inlet"]
    uin = inlet_json["u"]
    kin = inlet_json["k"]
    ωin = inlet_json["omega"]
    println("INLET INFO")
    println("\tu = $uin")
    println("\tk = $kin")
    println("\tω = $ωin")

    A, b, r, maxdiag = gpu_init_pressure_eq(dx, dy, sz, gc, nthread)
    ls = Ls(A, b, r, 1.2, maxdiag, 1000, 1e-4)
    println("EQ INFO")
    println("\tmax diag(A) = $maxdiag")

    x_h = [dx*(i - gc - 0.5) + xmin for i = 1:sz[1]]
    y_h = [dy*(j - gc - 0.5) + ymin for j = 1:sz[2]]
    x = CuArray(x_h)
    y = CuArray(y_h)
    u    = CUDA.zeros(Float64, sz...)
    v    = CUDA.zeros(Float64, sz...)
    ut   = CUDA.zeros(Float64, sz...)
    vt   = CUDA.zeros(Float64, sz...)
    uu   = CUDA.zeros(Float64, sz...)
    vv   = CUDA.zeros(Float64, sz...)
    p    = CUDA.zeros(Float64, sz...)
    k    = CUDA.zeros(Float64, sz...)
    ω    = CUDA.zeros(Float64, sz...)
    kt   = CUDA.zeros(Float64, sz...)
    ωt   = CUDA.zeros(Float64, sz...)
    nut  = CUDA.zeros(Float64, sz...)
    divU = CUDA.zeros(Float64, sz...)
    fill!(u  , uin)
    fill!(ut , uin)
    fill!(uu , uin)
    fill!(k  , kin)
    fill!(ω  , ωin)
    fill!(nut, kin/ωin)
    dfunc = prepare_dfunc(init_json["wind turbines"], x_h, y_h, dx, dy, sz)
    

    so = Solver(
        u, v, ut, vt, uu, vv, p,
        k, ω, kt, ωt, nut,
        divU, dfunc,
        uin, kin, ωin,
        ls,
        x, y, dx, dy, dt, Re,
        sz, gc, maxstep
    )

    output_path = init_json["output"]
    println("output to $output_path")
    flush(stdout)

    return so, output_path
end

function time_integral!(so::Solver)
    so.ut .= so.u
    so.vt .= so.v
    so.kt .= so.k
    so.ωt .= so.ω
    ls = so.ls

    gpu_pseudo_U!(
        so.u, so.v, so.ut, so.vt, so.uu, so.vv, so.nut, so.dfunc,
        ls.b, so.dx, so.dy, so.dt, 1/so.Re, ls.maxdiag,
        so.sz, so.gc, nthread
    )
    lsit, lserr = gpu_sor!(
        ls.A, so.p, ls.b, ls.r, ls.ω,
        so.sz, so.gc, ls.maxerr, ls.maxit, nthread
    )
    gpu_pbc!(
        so.p, so.sz, so.gc
    )
    gpu_project_p!(
        so.u, so.v, so.uu, so.vv, so.p,
        so.dx, so.dy, so.dt,
        so.sz, so.gc, nthread
    )
    gpu_Ubc!(
        so.u, so.v, so.ut, so.vt,
        so.uin, so.dx, so.dt,
        so.sz, so.gc
    )
    gpu_UUbc!(
        so.uu, so.vv, so.uin,
        so.sz, so.gc
    )
    gpu_kωSST2003!(
        so.k, so.ω, so.kt, so.ωt, so.u, so.v, so.uu, so.vv, so.nut,
        1/so.Re, 1.25, so.dx, so.dy, so.dt,
        so.sz, so.gc, nthread
    )
    gpu_kbc!(
        so.k, so.kt, so.kin, so.u,
        so.dx, so.dt, so.sz, so.gc
    )
    gpu_ωbc!(
        so.ω, so.ωt, so.ωin, so.u,
        so.dx, so.dt, so.sz, so.gc
    )
    divmag = gpu_divU!(
        so.uu, so.vv, so.divU,
        so.dx, so.dy,
        so.sz, so.gc, nthread
    )
    if divmag > 1 || isnan(divmag)
        error("cfd solver failed to converge, |div U|/N = $divmag")
    end
    return lsit, lserr, divmag
end

function write_csv(path::String, so::Solver)
    u = Array(so.u)
    v = Array(so.v)
    p = Array(so.p)
    k = Array(so.k)
    ω = Array(so.ω)
    nut = Array(so.nut)
    divu = Array(so.divU)
    x = Array(so.x)
    y = Array(so.y)
    sz = so.sz
    gc = so.gc
    x_coord = zeros(sz...)
    y_coord = zeros(sz...)
    z_coord = zeros(sz...)
    for j = 1:sz[2], i = 1:sz[1]
        x_coord[i, j] = x[i]
        y_coord[i, j] = y[j]
    end

    df = DataFrame(
        x = vec(x_coord),
        y = vec(y_coord),
        z = vec(z_coord)
    )

    df[!, :u] = vec(u)
    df[!, :v] = vec(v)
    df[!, :p] = vec(p)
    df[!, :k] = vec(k)
    df[!, :omega] = vec(ω)
    df[!, :nut] = vec(nut)
    df[!, :divu] = vec(divu)
    CSV.write(path, df)
    println("written to $path")
end
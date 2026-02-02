# cmd = `sh -c 'bin/csegfault'`
cmd = `sh -c 'julia segfault.jl'`
try
    io = IOBuffer()
    p = run(cmd, devnull, io, io; wait=false)
    wait(p)
    code = p.exitcode
    output = String(take!(io))
    println("code=$code\nmessage=$output")
catch e
    bt = catch_backtrace()
    s = sprint(Base.showerror, e, bt)
    println(s)
end
println("OK")
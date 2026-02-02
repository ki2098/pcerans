# cmd = `sh -c 'bin/csegfault'`
cmd = `sh -c 'julia segfault.jl'`
try
    io = IOBuffer()
    run(cmd, devnull, io, io)
    output = String(take!(io))
    print("cmd's message: $output")
catch e
    bt = catch_backtrace()
    s = sprint(Base.showerror, e, bt)
    println(s)
end
println("OK")
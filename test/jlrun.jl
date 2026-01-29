cmd = `ls xxx`
try
    io = IOBuffer()
    p = run(cmd, devnull, io, io; wait=false)
    wait(p)
    output = String(take!(io))
    print("cmd's message: $output")
catch e
    bt = catch_backtrace()
    s = sprint(Base.showerror, e, bt)
    println(s)
end
println("OK")
include("solve.jl")

using Dates

cd(@__DIR__)
println("now working in $(pwd())")

solver, output_path = init("setup.json")
println(now())
for step = 1:solver.maxstep
    try
        lsit, lserr, divmag = time_integral!(solver)
        print("\rIt=$step, divU=$divmag, LsIt=$lsit, LsErr = $lserr ")
        flush(stdout)
    catch e
        println(e)
        break
    end
end
println()
println(now())

write_csv(output_path, solver)
while true
    a = 0
    try
        a = randn()
        @assert a<0
        println("$a < 0")
        break
    catch e
        println("error = $e")
        println("$a >= 0")
    end
end
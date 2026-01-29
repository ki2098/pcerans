ptr = Ptr{Int}(0)
try
    unsafe_store!(ptr, 0)
    # throw(ErrorException("not good..."))
catch e
    bt = catch_backtrace()
    s = sprint(Base.showerror, e, bt)
    print("$s\n")
end

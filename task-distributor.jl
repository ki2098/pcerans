function split_view(x, n)
    l = length(x)
    base = div(l, n)
    extra = mod(l, n)
    views = Vector(undef, n)
    i = 1
    for k in 1:n
        if k <= extra
            len = base + 1
        else
            len = base
        end
        views[k] = @view x[i:i+len-1]
        i += len
    end
    return views
end
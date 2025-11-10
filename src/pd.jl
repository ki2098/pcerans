using CUDA
using JSON
using Intervals

function intersec_len(interv1, interv2)
    intersec = intersect(interv1, interv2)
    return last(intersec) - first(intersec)
end

function prepare_dfunc(wt_info_json, x, y, dx, dy, sz)
    diameter = wt_info_json["diameter"]
    model_diameter = 1.5*diameter
    thick = wt_info_json["thick"]
    C = wt_info_json["C"]
    wt_array_json = wt_info_json["coordinates"]
    wt_count = length(wt_array_json)
    println("PD MODEL INFO")
    println("\tC = $C")
    println("\tturbine thickness = $thick")
    println("\tturbine diameter = $diameter")
    println("\tmodeled turbine diameter = $model_diameter")
    for t = 1:wt_count
        wt_xy = wt_array_json[t]
        println("\tturbine $t at ($(wt_xy[1]), $(wt_xy[2]))")
    end
    dfunc = zeros(sz...)
    for j = 1:sz[2], i = 1:sz[1]
        cell_x_range = Interval(x[i] - 0.5*dx, x[i] + 0.5*dx)
        cell_y_range = Interval(y[j] - 0.5*dy, y[j] + 0.5*dy)
        for t = 1:wt_count
            wt_xy = wt_array_json[t]
            wtx = wt_xy[1]
            wty = wt_xy[2]
            
            wt_x_range = Interval(wtx - 0.5*thick, wtx + 0.5*thick)
            wt_y_range = Interval(wty - 0.5*model_diameter, wty + 0.5*model_diameter)
            cover = intersec_len(cell_x_range, wt_x_range)*intersec_len(cell_y_range, wt_y_range)
            if cover > 0
                sigma = model_diameter/6.0
                dist = abs(y[j] - wty)
                gauss_kernel = exp(- 0.5*(dist/sigma)^2)
                occupancy = cover/(dx*dy)
                dfunc[i, j] = C*gauss_kernel*occupancy
            end
        end
    end
    dfunc_d = CuArray(dfunc)
    return dfunc_d
end
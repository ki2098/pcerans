using LinearAlgebra
using CUDA

function kernel_pressure_eq_A!(A, dx, dy, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        Ae = Aw = 1/(dx^2)
        An = As = 1/(dy^2)
        Ac = - (Ae + Aw + An + As)
        A[i, j, 1] = Ac
        A[i, j, 2] = Ae
        A[i, j, 3] = Aw
        A[i, j, 4] = An
        A[i, j, 5] = As
    end
    nothing
end

function gpu_pressure_eq_A(dx, dy, sz, gc, nthread)
    nblock = (cld(sz[1], nthread[1]), cld(sz[2], nthread[2]))
    A = CUDA.zeros(Float64, sz..., 5)
    A[:, :, 1] .= 1
    @cuda threads=nthread blocks=nblock kernel_pressure_eq_A!(
        A, dx, dy, sz, gc
    )
    maxdiag = maximum(abs.(A[:, :, 1]))
    A ./= maxdiag
    return A, maxdiag
end

function gpu_init_pressure_eq(dx, dy, sz, gc, nthread)
    A, maxdiag = gpu_pressure_eq_A(dx, dy, sz, gc, nthread)
    b = CUDA.zeros(Float64, sz...)
    r = CUDA.zeros(Float64, sz...)
    return A, b, r, maxdiag
end

function cell_residual(A, x, b, i, j)
    xc = x[i, j]
    xe = x[i + 1, j]
    xw = x[i - 1, j]
    xn = x[i, j + 1]
    xs = x[i, j - 1]
    Ac = A[i, j, 1]
    Ae = A[i, j, 2]
    Aw = A[i, j, 3]
    An = A[i, j, 4]
    As = A[i, j, 5]
    r  = b[i, j] - (Ac*xc + Ae*xe + Aw*xw + An*xn + As*xs)
    return r
end

function kernel_residual!(A, x, b, r, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        r[i, j] = cell_residual(A, x, b, i, j)
    end
    nothing
end

function kernel_colored_sor_sweep!(A, x, b, ω, sz, gc, c)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        if (i + j)%2 == c
            x[i, j] += ω*cell_residual(A, x, b, i, j)/A[i, j, 1]
        end
    end
    nothing
end

function gpu_sor!(A, x, b, r, ω, sz, gc, maxerr, maxit, nthread)
    nblock = (cld(sz[1], nthread[1]), cld(sz[2], nthread[2]))
    it = 0
    errmag = 0
    ncell_effective = (sz[1] - 2*gc)*(sz[2] - 2*gc)
    while true
        @cuda threads=nthread blocks=nblock kernel_colored_sor_sweep!(
            A, x, b, ω, sz, gc, 0
        )
        @cuda threads=nthread blocks=nblock kernel_colored_sor_sweep!(
            A, x, b, ω, sz, gc, 1
        )
        @cuda threads=nthread blocks=nblock kernel_residual!(
            A, x, b, r, sz, gc
        )
        errmag = norm(r)
        errmag = errmag / sqrt(ncell_effective)
        it += 1
        if errmag <= maxerr || it >= maxit
            break
        end
        if errmag > 1 || isnan(errmag)
            error("linear solver failed to converge, |r|/N = $errmag")
        end
    end
    return it, errmag
end
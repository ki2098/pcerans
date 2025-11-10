using LinearAlgebra
using CUDA

function utopia_convection(fww, fw, fc, fe, fee, uW, uc, uE, dx)
    xE = uE*(fw - 27*fc + 27*fe - fee)/(24*dx)
    xW = uW*(fww - 27*fw + 27*fc - fe)/(24*dx)
    xA = abs(uc)*(fww - 4*fw + 6*fc - 4*fe + fee)/(12*dx)
    return 0.5*(xE + xW) + xA
end

function cell_convection(f, u, v, uu, vv, dx, dy, i, j)
    fc  = f[i, j]
    fe  = f[i+1, j]
    fee = f[i+2, j]
    fw  = f[i-1, j]
    fww = f[i-2, j]
    fn  = f[i, j+1]
    fnn = f[i, j+2]
    fs  = f[i, j-1]
    fss = f[i, j-2]
    uc  = u[i  , j]
    uE = uu[i  , j]
    uW = uu[i-1, j]
    vc  = v[i, j  ]
    vN = vv[i, j  ]
    vS = vv[i, j-1]
    udfdx = utopia_convection(fww, fw, fc, fe, fee, uW, uc, uE, dx)
    vdfdy = utopia_convection(fss, fs, fc, fn, fnn, vS, vc, vN, dy)
    return udfdx + vdfdy
end

function cell_diffusion(f, nu, dx, dy, i, j)
    fc = f[i, j]
    fe = f[i + 1, j]
    fw = f[i - 1, j]
    fn = f[i, j + 1]
    fs = f[i, j - 1]
    d2fdx2 = (fe - 2*fc + fw)/(dx^2)
    d2fdy2 = (fn - 2*fc + fs)/(dy^2)
    return nu*(d2fdx2 + d2fdy2)
end

function kernel_pseudo_U!(u, v, ut, vt, uu, vv, nut, dfunc, dx, dy, dt, nu, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        nu_eff = nu + nut[i, j]
        u_conv = cell_convection(ut, ut, vt, uu, vv, dx, dy, i, j)
        u_diff = cell_diffusion(ut, nu_eff, dx, dy, i, j)
        v_conv = cell_convection(vt, ut, vt, uu, vv, dx, dy, i, j)
        v_diff = cell_diffusion(vt, nu_eff, dx, dy, i, j)
        uc = ut[i, j]
        vc = vt[i, j]
        UMag = sqrt(uc^2 + vc^2)
        fx = dfunc[i, j]*UMag*uc
        fy = dfunc[i, j]*UMag*vc
        u[i, j] = uc + dt*(- u_conv + u_diff - fx)
        v[i, j] = vc + dt*(- v_conv + v_diff - fy)
    end
    nothing
end

function kernel_interpolate_UU!(u, v, uu, vv, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        uu[i, j] = (u[i, j] + u[i+1, j])/2
    end
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc-1
        vv[i, j] = (v[i, j] + v[i, j+1])/2
    end
    nothing
end

function cell_divU(uu, vv, dx, dy, i, j)
    uE = uu[i    , j]
    uW = uu[i - 1, j]
    vN = vv[i, j    ]
    vS = vv[i, j - 1]
    return (uE - uW)/dx + (vN - vS)/dy
end

function kernel_pressure_eq_b!(uu, vv, b, dx, dy, dt, maxdiag, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        b[i, j] = cell_divU(uu, vv, dx, dy, i, j)/(dt*maxdiag)
    end
    nothing
end

function kernel_project_p!(u, v, uu, vv, p, dx, dy, dt, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        dpdx = (p[i+1, j] - p[i-1, j])/(2*dx)
        dpdy = (p[i, j+1] - p[i, j-1])/(2*dy)
        u[i, j] -= dt*dpdx
        v[i, j] -= dt*dpdy
    end
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        dpdx = (p[i+1, j] - p[i, j])/dx
        uu[i, j] -= dt*dpdx
    end 
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc-1
        dpdy = (p[i, j+1] - p[i, j])/dy
        vv[i, j] -= dt*dpdy
    end
    nothing
end

function kernel_divU!(uu, vv, divU, dx, dy, sz, gc)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    j = (blockIdx().y - 1)*blockDim().y + threadIdx().y
    if gc < i <= sz[1]-gc && gc < j <= sz[2]-gc
        divU[i, j] = cell_divU(uu, vv, dx, dy, i, j)
    end
    nothing
end

function gpu_pseudo_U!(u, v, ut, vt, uu, vv, nut, dfunc, b, dx, dy, dt, nu, maxdiag, sz, gc, nthread)
    nblock = (cld(sz[1], nthread[1]), cld(sz[2], nthread[2]))
    @cuda threads=nthread blocks=nblock kernel_pseudo_U!(
        u, v, ut, vt, uu, vv, nut, dfunc,
        dx, dy, dt, nu, sz, gc
    )
    @cuda threads=nthread blocks=nblock kernel_interpolate_UU!(
        u, v, uu, vv, sz, gc
    )
    @cuda threads=nthread blocks=nblock kernel_pressure_eq_b!(
        uu, vv, b, dx, dy, dt, maxdiag, sz, gc
    )
end

function gpu_project_p!(u, v, uu, vv, p, dx, dy, dt, sz, gc, nthread)
    nblock = (cld(sz[1], nthread[1]), cld(sz[2], nthread[2]))
    @cuda threads=nthread blocks=nblock kernel_project_p!(
        u, v, uu, vv, p, dx, dy, dt, sz, gc
    )
end

function gpu_divU!(uu, vv, divU, dx, dy, sz, gc, nthread)
    nblock = (cld(sz[1], nthread[1]), cld(sz[2], nthread[2]))
    @cuda threads=nthread blocks=nblock kernel_divU!(
        uu, vv, divU, dx, dy, sz, gc
    )
    divmag = norm(divU)
    ncell_effective = (sz[1] - 2*gc)*(sz[2] - 2*gc)
    return divmag/sqrt(ncell_effective)
end

# checked
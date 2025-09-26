# --- Script to read in snapshot, calculate potential and DF ---

module MakeGlobals

# Import necessary packages
using CSV
using DataFrames
using Interpolations
using Statistics

# Physical constants export
export G, M, b, M_bh, N, m_test, m_field, coulomb_log

# Data, extrema export
export m_tab, r_tab, vr_tab, vt_tab, v_tab
export rmin, rmax, vmin, vmax
export vrmin, vrmax, vtmin, vtmax

# Psi export
export find_shell, find_psi_arrays, find_psi_exact, psi_exact
export psi_tab, M_arr, psi_calc, psi_th

# DF export
export get_sampling_DF_rv, bilinear_interp
export tabr_samp, tabv_samp, tabN_rv, tabN_rv_grid, tabF_rv, tabF_rv_grid, tabDFth_rv, tabDFth_rv_grid
export N_1D, F_1D, F_th
export norm_n, norm_f, N_1D_, F_1D_
export read_file, midpoint2D

# Numerical constant export
export dr, dv, int_steps, tol, epsilon


# --- Physical constants ---
const G = 1
const M = 1
const M_bh = 0 #8.500250825e-03
const b = 3*pi/16

# --- Numerical constants ---
const int_steps=500
const tol=1e-8
const epsilon = 0.001

const dr=0.02
const dv=0.02


### Potential calculation ###


# --- Shell finding via bisection---
function find_shell(rr,r_tab)

    #check edge cases
    if rr < r_tab[1]
        return 0
    elseif rr > r_tab[N]
        return N
    elseif rr == r_tab[N]
        return N-1
    end

    #simple bisection
    left, right = 1, N - 1
    while left < right
        mid = div(left + right + 1, 2)
        if r_tab[mid] <= rr
            left = mid
        else
            right = mid - 1
        end
    end

    return left
end

# --- Generate psi for field ---
function find_psi_arrays(r_tab,m_tab)

    # set up array
    psi_arr=zeros(N+1)
    M_arr=zeros(N+1)

    psi_arr[N+1] = 0.0 # Eq 1.7
    M_arr[N] = sum(m_tab) # Eq 1.8
    psi_arr[N] = -G*M_arr[N] / r_tab[N] 

    # iterate over shells
    for k in (N-1):-1:1
        M_arr[k] = M_arr[k+1]-m_tab[k+1] # Eq 1.10
        psi_arr[k] = psi_arr[k+1] - G*M_arr[k]*(1/r_tab[k] - 1/r_tab[k+1]) # Eq 1.9
    end

    # return arrays
    return psi_arr, M_arr
end

# --- Find exact potential (including BH) using arrays ---
function find_psi_exact(rr,r_tab,psi_arr,M_arr)
    k=find_shell(rr,r_tab) # find shell

    if k == 0 # 1st shell
        return psi_arr[1] - G * M_bh / rr

    elseif k == N # last shell
        return -G * (M_arr[N] + M_bh)/rr

    else # between
        rk, rk1 = r_tab[k], r_tab[k+1]
        psi_k, psi_k1 = psi_arr[k], psi_arr[k+1]
        return psi_k+((1/rk - 1/rr)/(1/rk - 1/rk1))*(psi_k1-psi_k) - G*M_bh/rr #Eq 1.15
    end
end

# --- Wrapper function for psi ---
function psi_exact(r_tab,psi_arr,M_arr)
    return rr -> find_psi_exact(rr,r_tab,psi_arr,M_arr)
end


### Distribution function calculation ###


# --- Plummer potential ---
function psi_th(r)
    return -G*M/sqrt(r^2+b^2)
end

# --- Plummer distribution function ---
function DF_E_theoretical(E)
    if (E<=0)
        return 24*sqrt(2)/(7*pi^3) * b^2/(G^5*M^4) * (-E)^(7/2)
    else 
        return 0.0
    end
end

# --- Plummer distribution function ---
function DF_r_v_theoretical(r, v)
    E = psi_th(r) + 0.5*v^2
    return 16*pi^2*r^2*v^2*DF_E_theoretical(E)
end


# --- Sampled distribution function ---
function get_DF_lin(tabr, tabv, tabm)

    # linear bin centers
    tabr_samp = collect(0.0:dr:rmax-dr) .+ dr/2 
    tabv_samp = collect(0.0:dv:vmax-dv) .+ dv/2

    # number of bins 
    nbr = length(tabr_samp)
    nbv = length(tabv_samp)

    # storage arrays
    tabF_rv = zeros(Float64, nbr, nbv) # matrix
    tabF_rv_grid = zeros(Float64, nbr * nbv) # 1D

    tabN_rv = zeros(Float64, nbr, nbv)
    tabN_rv_grid = zeros(Float64, nbr * nbv)

    tabDFth_rv = zeros(Float64, nbr, nbv)
    tabDFth_rv_grid = zeros(Float64, nbr * nbv)

    # emperical DF
    Threads.@threads for i=1:N

        r = tabr[i]
        v = tabv[i]
        mass = tabm[i]
        
        # determine indices (bottom left)
        ir = floor(Int64, r/dr) + 1
        iv = floor(Int64, v/dv) + 1

        # check within cell
        if (ir <= nbr) && (iv <= nbv) && (ir > 0) && (iv > 0)

            # 1D index
            index = iv + (ir-1)*nbv
            
            tabN_rv[ir, iv] += mass/(dr*dv) # Eq 1.27
            tabN_rv_grid[index] += mass/(dr*dv)

            tabF_rv[ir, iv] += (mass^2)/(dr*dv) # Eq 1.28
            tabF_rv_grid[index] += (mass^2)/(dr*dv)
        end
    end

    # theoretical DF 
    Threads.@threads for index=1:nbr*nbv
        
        # bin indices
        ir = floor(Int64, (index-1)/nbv) + 1
        iv = index - (ir-1)*nbv

        # determine r, v; DF
        r = tabr_samp[ir]
        v = tabv_samp[iv]

        df = DF_r_v_theoretical(r, v)
        
        # save to arrays 
        tabDFth_rv[ir, iv] = df
        tabDFth_rv_grid[index] = df
    end

    return tabr_samp, tabv_samp,
           tabN_rv, tabN_rv_grid,
           tabF_rv, tabF_rv_grid,
           tabDFth_rv, tabDFth_rv_grid 
end

function get_DF_log(tabr,tabv,tabm)

    # log edges
    logr_edge = collect(log(rmin):dlogr:log(rmax))
    logv_edge = collect(log(vmin):dlogv:log(vmax))

    # determine bin centers -> log space (!!)
    logr_samp = (logr_edge[1:end-1] .+ logr_edge[2:end]) ./ 2
    logv_samp = (logv_edge[1:end-1] .+ logv_edge[2:end]) ./ 2

    # number of bins
    nbr = length(logr_samp)
    nbv = length(logv_samp)

    # storage arrays
    tabF_rv = zeros(Float64, nbr, nbv) # matrix
    tabF_rv_grid = zeros(Float64, nbr * nbv) # 1D

    tabN_rv = zeros(Float64, nbr, nbv)
    tabN_rv_grid = zeros(Float64, nbr * nbv)

    tabDFth_rv = zeros(Float64, nbr, nbv)
    tabDFth_rv_grid = zeros(Float64, nbr * nbv)

    # emperical DF

    Threads.@threads for i=1:N
        r = tabr[i]
        v = tabv[i]
        mass = tabm[i]

        r_log = log(r)
        v_log = log(v)

        # determine indices (bottom left)
        ir = searchsortedlast(logr_edge,r_log)
        iv = searchsortedlast(logv_edge,v_log)

        if (ir <= nbr) && (iv <= nbv) && (ir > 0) && (iv > 0)

            # 1D index
            index = iv + (ir-1)*nbv
        
            # eq 1.34, 1.35
            tabN_rv[ir, iv] += mass/(dlogr*dlogv) # Eq 1.40
            tabN_rv_grid[index] += mass/(dlogr*dlogv)

            tabF_rv[ir, iv] += (mass^2)/(dlogr*dlogv) # Eq 1.41
            tabF_rv_grid[index] += (mass^2)/(dlogr*dlogv)
        end
    end 

    # theoretical DF
    Threads.@threads for index=1:nbr*nbv

        # bin indices
        ir = floor(Int64, (index-1)/nbv) + 1
        iv = index - (ir-1)*nbv

        # determine r, v; DF
        lnr = logr_samp[ir]
        lnv = logv_samp[iv]

        r = exp(lnr)
        v = exp(lnv)

        # jacobian correction -> F(logr,logv) = r * v * F(r,v)
        df = r * v * DF_r_v_theoretical(r,v)

        tabDFth_rv[ir,iv] = df
        tabDFth_rv_grid[index] = df
    end

    return logr_samp, logv_samp,
           tabN_rv, tabN_rv_grid,
           tabF_rv, tabF_rv_grid,
           tabDFth_rv, tabDFth_rv_grid
end


# --- Bilinear interpolation (weighted means) --- 
function bilinear_interp(rr, vv, rgrid, vgrid, fgrid)

    # calculate grid length 
    nr = length(rgrid)
    nv = length(vgrid)

    # find index of left/bottom center: rgrid[i] <= rr < rgrid[i+1]
    ir = searchsortedlast(rgrid, rr)
    iv = searchsortedlast(vgrid, vv)

    # bounds check 
    if ir < 1 || iv < 1 || ir >= nr || iv >= nv
        return 0.0
    end

    @inbounds begin # safe since above check

        # use adjacent centers as interpolation points:
        x1 = rgrid[ir];   x2 = rgrid[ir+1]
        y1 = vgrid[iv];   y2 = vgrid[iv+1]

        f11 = fgrid[ir, iv]
        f12 = fgrid[ir, iv+1]
        f21 = fgrid[ir+1, iv]
        f22 = fgrid[ir+1, iv+1]

        denom = (x2 - x1) * (y2 - y1)

        return ( (x2 - rr)*(y2 - vv)*f11 +
                 (x2 - rr)*(vv - y1)*f12 +
                 (rr - x1)*(y2 - vv)*f21 +
                 (rr - x1)*(vv - y1)*f22 ) / denom
    end
end

# --- 2D Integration for normalization 
function midpoint2D(func, rmin, rmax, vmin, vmax; Nr=1000, Nv=1000)
    dr = (rmax - rmin) / Nr
    dv = (vmax - vmin) / Nv
    integral = 0.0

    for i in 0:Nr-1
        r0 = rmin + i*dr
        r1 = r0 + dr
        rmid = (r0 + r1) / 2

        for j in 0:Nv-1
            v0 = vmin + j*dv
            v1 = v0 + dv
            vmid = (v0 + v1) / 2

            integral += func(rmid, vmid) * dr * dv
        end
    end

    return integral
end

# --- Read in data, compute psi and DFs ---
function read_file(filename)

    # read in data, make DataFrame
    df=CSV.read(filename,DataFrame)

    # set up arrays, calculate velocity
    #global id_tab, m_tab, r_tab, vr_tab, vt_tab, startype_tab, binflag = eachcol(df[:,1:7])
    global r_tab, vr_tab, vt_tab = eachcol(df[:,1:3])
    global v_tab = sqrt.(vr_tab .^ 2 + vt_tab .^2)

    # determine array length
    global N = length(r_tab)
    global m_tab=fill(1/N,N)

    # find (r,v) bounds
    global rmin, rmax = extrema(r_tab)
    global vmin, vmax = extrema(v_tab)
    global vrmin, vrmax = extrema(vr_tab)
    global vtmin, vtmax = extrema(vt_tab)

    # calculate coulomb logarithm, test mass
    global coulomb_log = log(0.11*N) # M_bh / mean(m_tab)
    global m_test = 1/N # mean(m_tab)
    global m_field = m_test

    # find potential via Henon
    global psi_tab, M_arr = find_psi_arrays(r_tab,m_tab)
    global psi_calc = psi_exact(r_tab,psi_tab,M_arr)

    # compute DataFrame arrays
    global tabr_samp, tabv_samp, tabN_rv, tabN_rv_grid, tabF_rv, tabF_rv_grid, tabDFth_rv, tabDFth_rv_grid = get_DF_lin(r_tab,v_tab,m_tab)

    # compute distribution functions
    global N_1D_ = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabN_rv)
    global F_1D_ = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabF_rv)
    global F_th = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabDFth_rv)

    # normalize
    global norm_n = midpoint2D(N_1D_,rmin,rmax,vmin,vmax)
    global norm_f = midpoint2D(F_1D_,rmin,rmax,vmin,vmax)

    global N_1D = (r,v) -> (1/norm_n)*N_1D_(r,v)
    global F_1D = (r,v) -> (1/(norm_f*N))*F_1D_(r,v)
    

end

end

# --- Distribution Function Sampling --- #

struct SampDF
    tabr_samp::Vector{Float64}
    tabv_samp::Vector{Float64}
    tabn_rv::Matrix{Float64}
    tabf_rv::Matrix{Float64}
    N::Function
    F::Function
end

## Calculate (logr, logv) DF, convolve with Gaussian kernel ##
function kde_rv_df(dat::SnapDat)

    N = length(dat.r)

    # Number of bins 
    nbr = floor(Int64, (log(dat.rmax) - log(dat.rmin)) / dlogr)
    nbv = floor(Int64, (log(dat.vmax) - log(dat.vmin)) / dlogv)

    # Bin centers (log space)
    logr_samp = log(dat.rmin) .+ (collect(0:nbr-1) .+ 0.5) .* dlogr
    logv_samp = log(dat.vmin) .+ (collect(0:nbv-1) .+ 0.5) .* dlogv

    # Storage matrices 
    tabn_rv = zeros(Float64, nbr, nbv)
    tabf_rv = zeros(Float64, nbr, nbv)

    Threads.@threads for i in 1:N # Parallelise over all objects

        # Characterize object 
        r, v, mass = dat.r[i], dat.v[i], dat.m[i]

        # Locate (r,v) bin in log space 
        ir = floor(Int64, (log(r) - log(dat.rmin)) / dlogr) + 1
        iv = floor(Int64, (log(v) - log(dat.vmin)) / dlogv) + 1

        # Only sum if particle falls inside bin / grid 
        if (ir in 1:nbr) && (iv in 1:nbv) && (ir > 0) && (iv > 0)
            @inbounds tabn_rv[ir,iv] += mass / (dlogr * dlogv)
            @inbounds tabf_rv[ir,iv] += mass^2 / (dlogr * dlogv)
        end
    end

    ## KDE Smoothing ##
    sigma_r, sigma_v = 2.0, 2.0 # Kernel width (# bins)
    half_w = 4

    dr_offsets = -half_w:half_w # Offsets for kernel 
    dv_offsets = -half_w:half_w

    # Compute 2D Gaussian kernel; normalize to 1
    kernel = [exp(-0.5 * ((dr/sigma_r)^2 + (dv/sigma_v)^2))
              for dr in dr_offsets, dv in dv_offsets]
    kernel ./= sum(kernel)

    # Convolve DF with kernel 
    tabn_kde = conv2d(tabn_rv, kernel)
    tabf_kde = conv2d(tabf_rv, kernel)

    # Interpolated functions for coeffs
    mass_dens = (r,v) -> bilinear_interp(r, v, logr_samp, logv_samp, tabn_kde)
    sq_mass_dens = (r,v) -> bilinear_interp(r, v, logr_samp, logv_samp, tabf_kde)

    return SampDF(logr_samp, logv_samp, tabn_kde, tabf_kde, mass_dens, sq_mass_dens)
end

## Simple 2D convolution for KDE ##
function conv2d(A::Matrix, K::Matrix)
    
    n_row, n_col = size(A)
    half_row, half_col = (size(K,1)-1)÷2, (size(K,2)-1)÷2

    # Storage matrix 
    mat = zeros(n_row, n_col)

    for row in 1:n_row, col in 1:n_col
        acc = 0.0

        for drow in -half_row:half_row, dcol in -half_col:half_col
            ii, jj = row + drow, col + dcol

            # Skip if outside of matrix 
            if (1 <= ii <= n_row) && (1 <= jj <= n_col)
                # Correction for 1-based indexing 
                @inbounds acc += K[drow+half_row+1, dcol+half_col+1] * A[ii, jj]
            end
        end

        mat[row, col] = acc
    end

    return mat
end

 ## Sample DF in (logr,logv) ##
function log_rv_df(dat::SnapDat)

    N = length(dat.r)

    # Number of bins
    nbr = floor(Int64, (log(dat.rmax) - log(dat.rmin)) / dlogr)
    nbv = floor(Int64, (log(dat.vmax) - log(dat.vmin)) / dlogv)

    # Log bin centers
    logr_samp = log(dat.rmin) .+ (collect(0:nbr-1) .+ 0.5) .* dlogr
    logv_samp = log(dat.vmin) .+ (collect(0:nbv-1) .+ 0.5) .* dlogv

    # Initialize matrices
    tabn_rv = zeros(Float64, nbr, nbv)
    tabf_rv = zeros(Float64, nbr, nbv)

    Threads.@threads for i in 1:N # Iterate over all objects

        # Find object r, v, mass
        r = dat.r[i]
        v = dat.v[i]
        mass = dat.m[i]

        # (r,v) <-> (logr,logv)
        r_log = log(r)
        v_log = log(v)

        # Find bin indices
        ir = floor(Int64, (r_log - log(dat.rmin)) / dlogr) + 1
        iv = floor(Int64, (v_log - log(dat.vmin)) / dlogv) + 1

        # Check within bin, sum
        if (ir <= nbr) && (iv <= nbv) && (ir > 0) && (iv > 0)
            @inbounds tabn_rv[ir,iv] += mass/(dlogr*dlogv) # Eq 1.38
            @inbounds tabf_rv[ir,iv] += (mass^2)/(dlogr*dlogv) # Eq 1.39
        end
    end

    N_samp = (r,v) -> bilinear_interp(r,v,logr_samp,logv_samp,tabn_rv)
    F_samp = (r,v) -> bilinear_interp(r,v,logr_samp,logv_samp,tabf_rv)

    # Return DF (log) struct
    return SampDF(logr_samp,logv_samp,tabn_rv,tabf_rv,N_samp,F_samp)
end

## Sample DF in (r,v) ##
function linear_rv_df(dat::SnapDat)

    N = length(dat.r)

    # Initialize arrays from 0 -> (rmax-dr); Bin centers
    tabr_samp = collect(0.0:dr:dat.rmax-dr) .+ dr/2
    tabv_samp = collect(0.0:dv:dat.vmax-dv) .+ dv/2

    # Number of bins
    nbr, nbv = length(tabr_samp), length(tabv_samp)

    # Initialize matrices
    tabn_rv = zeros(Float64, nbr, nbv)
    tabf_rv = zeros(Float64, nbr, nbv)

    Threads.@threads for i=1:N # Iterate over all objects

        # r = (ir-1)*dr + deltar
        # v = (iv-1)*dv + deltav
        # index-1 = iv-1 + (ir-1)*nbv

        # Find object r, v, mass
        r = dat.r[i]
        v = dat.v[i]
        mass = dat.m[i]

        # Find bin indices
        ir = floor(Int64, r/dr) + 1
        iv = floor(Int64, v/dv) + 1

        # Check within bin, sum
        if (ir <= nbr) && (iv <= nbv) && (ir > 0) && (iv > 0)
            @inbounds tabn_rv[ir,iv] += mass/(dr*dv) # Eq 1.27
            @inbounds tabf_rv[ir,iv] += (mass^2)/(dr*dv) # Eq 1.28
        end
    end

    # Interpolate DF
    N_samp = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabn_rv)
    F_samp = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabf_rv)

    # Return DF (linear) struct
    return SampDF(tabr_samp,tabv_samp,tabn_rv,tabf_rv,N_samp,F_samp)
end

# Compute rho(r) from N(logr, logv) # 
function rho_calc(sdf::SampDF, dat::SnapDat)
    
    # Log bin centers, number of bins
    logr_samp = sdf.tabr_samp
    nbr = length(logr_samp)

    # Initialize density array
    tab_rho = zeros(Float64, nbr)

    Threads.@threads for ir in 1:nbr

        r = exp(logr_samp[ir]) # r from bin centers

        # f(r,v) = f_tilde(logr,logv) / (r*v)
        tab_rho[ir] = 4π * r * midpoint(v -> sdf.N(log(r), log(v)) * v, dat.vmin, dat.vmax, n_steps)
    end

    # Lin interp rho(r) in log space
    rho_samp = function(r)
        logr = log(r)
        ir = searchsortedlast(logr_samp, logr)
        ir = clamp(ir, 1, nbr - 1)
        t = (logr - logr_samp[ir]) / (logr_samp[ir+1] - logr_samp[ir])
        return tab_rho[ir] + t * (tab_rho[ir+1] - tab_rho[ir])
    end

    return rho_samp
end

## Sample DF in (E,j) ##
function linear_jE_df(dat::SnapDat)

    # Calculate energy, momentum 
    Etab1 = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
    Ltab1 = dat.r .* dat.vt

    Etab = Etab1[Etab1 .< 0] # Filter unbound orbits
    Ltab = Ltab1[Etab1 .< 0]

    jtab = zeros(length(Etab))
    
    Threads.@threads for i in eachindex(Etab)
        jtab[i] = Ltab[i] / last(find_Lc(Etab[i]))
    end

    # Determine bounds 
    E_min, E_max = extrema(Etab)
    j_min, j_max = extrema(jtab)

    # Number of bins
    nbE = floor(Int64, (E_max - E_min) / dE)
    nbj = floor(Int64, (j_max - j_min) / dj)

    # Bin centers
    E_samp = E_min .+ (collect(0:nbE-1) .+ 0.5) .* dE
    j_samp = j_min .+ (collect(0:nbj-1) .+ 0.5) .* dj

    # Initialize matrices
    tabn_jE = zeros(Float64, nbE, nbj)
    tabf_jE = zeros(Float64, nbE, nbj)

    Threads.@threads for i in eachindex(Etab) # Loop over all bound objects 
    
    # Find object E, j, mass
        E = Etab[i]
        j = jtab[i]
        mass = dat.m[i]

    # Find bin indices
        iE = floor(Int64, (E - E_min) / dE) + 1
        ij = floor(Int64, (j - j_min) / dj) + 1

    # Check within bin, sum
        if (iE <= nbE) && (ij <= nbj) && (iE > 0) && (ij > 0)
            @inbounds tabn_jE[iE,ij] += mass/(dE*dj)
            @inbounds tabf_jE[iE,ij] += (mass^2)/(dE*dj)
        end
    end

    # Interpolated functions 
    N_samp = (E,j) -> bilinear_interp(E,j,E_samp,j_samp,tabn_jE)
    F_samp = (E,j) -> bilinear_interp(E,j,E_samp,j_samp,tabf_jE)

    # Return DF (linear) struct
    return SampDF(E_samp,j_samp,tabn_jE,tabf_jE,N_samp,F_samp)
end

## Sample DF in [log(-E),log(j)] ##
function log_jE_df(dat::SnapDat)

    # Calculate energy, momentum 
    Etab1 = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
    Ltab1 = dat.r .* dat.vt

    Etab = Etab1[Etab1 .< 0] # Filter unbound orbits
    Ltab = Ltab1[Etab1 .< 0]

    jtab = zeros(length(Etab))

    Threads.@threads for i in eachindex(Etab)
        jtab[i] = Ltab[i] / last(find_Lc(Etab[i]))
    end

    # Convert into log space, determine bounds 
    logE_tab = log.(.-Etab)
    logj_tab = log.(jtab)

    logE_min, logE_max = extrema(logE_tab)
    logj_min, logj_max = extrema(logj_tab[isfinite.(logj_tab)]) # Protect against j = 0

    # Number of bins
    nbE = floor(Int64, (logE_max - logE_min) / dlogE)
    nbj = floor(Int64, (logj_max - logj_min) / dlogj)

    # Log bin centers
    logE_samp = logE_min .+ (collect(0:nbE-1) .+ 0.5) .* dlogE
    logj_samp = logj_min .+ (collect(0:nbj-1) .+ 0.5) .* dlogj

    # Initialize matrices
    tabn_jE = zeros(Float64, nbE, nbj)
    tabf_jE = zeros(Float64, nbE, nbj)

    Threads.@threads for i in eachindex(Etab) # Loop over all bound objects 
        
        # Find object log(-E), logj, mass
        logE = logE_tab[i]
        logj = logj_tab[i]
        mass = dat.m[i]
        
        # Find bin indices
        iE = floor(Int64, (logE - logE_min) / dlogE) + 1
        ij = floor(Int64, (logj - logj_min) / dlogj) + 1
        
        # Check within bin, sum
        if (iE <= nbE) && (ij <= nbj) && (iE > 0) && (ij > 0)
            @inbounds tabn_jE[iE,ij] += mass/(dlogE*dlogj)
            @inbounds tabf_jE[iE,ij] += (mass^2)/(dlogE*dlogj)
        end
    end

    # Interpolated functions 
    N_samp = (E,j) -> bilinear_interp(E,j,logE_samp,logj_samp,tabn_jE)
    F_samp = (E,j) -> bilinear_interp(E,j,logE_samp,logj_samp,tabf_jE)

    # Return DF (log) struct
    return SampDF(logE_samp,logj_samp,tabn_jE,tabf_jE,N_samp,F_samp)
end






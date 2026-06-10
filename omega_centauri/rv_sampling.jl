# --- A Branch of Distribution Function Sampling with Orbit Sampling in R,V --- #

struct SampDF
    tabr_samp::Vector{Float64}
    tabv_samp::Vector{Float64}
    tabn_rv::Matrix{Float64}
    tabf_rv::Matrix{Float64}
    N::Function
    F::Function
end


using Random

## Orbit sampling in (r, v)
##
## Input:
##     r        = current radius
##     vr       = current radial velocity
##     vt       = current tangential velocity
##     N_orb    = number of orbit samples
##     psi_calc = exact potential function, e.g.
##                psi_calc = psi_exact(dat.r, psi_tab, M_tab)
##
## Output:
##     (rs = rs, vs = vs)
##
## where:
##     rs = sampled radii along the orbit
##     vs = corresponding total speed magnitudes

function orbit_sample_rv(
    r::Float64,
    vr::Float64,
    vt::Float64,
    N_orb::Integer,
    psi_calc::Function;
    rng_seed::Union{Nothing, Integer} = nothing,
    n_grid::Integer = 1000,
    tol::Float64 = 1e-15,
    expand_factor::Float64 = 2.0,
    max_expand::Integer = 100,
)

    # Simple Checks
    
    if N_orb < 1
        throw(ArgumentError("N_orb must be at least 1."))
    end

    if n_grid < 10
        throw(ArgumentError("n_grid must be at least 10."))
    end

    if !(isfinite(r) && isfinite(vr) && isfinite(vt) && r > 0.0)
        throw(ArgumentError(
            "Invalid input to orbit_sample_rv: " *
            "r = $r, vr = $vr, vt = $vt. " *
            "Expected finite values and r > 0."
        ))
    end

    if expand_factor <= 1.0
        throw(ArgumentError("expand_factor must be larger than 1."))
    end

    # RNG
    rng = rng_seed === nothing ? Random.default_rng() : MersenneTwister(rng_seed)

    # Orbital constants
    v2_now = vr^2 + vt^2
    E = 0.5 * v2_now + psi_calc(r)
    J = r * abs(vt)

    if !(isfinite(E) && isfinite(J))
        throw(ArgumentError(
            "Invalid orbital constants: E = $E, J = $J."
        ))
    end

    if J <= 0.0
        throw(ArgumentError(
            "Angular momentum J = $J is not positive. " 
        ))
    end
#=
    if E >= 0.0
        throw(ArgumentError(
            "Orbit appears unbound: E = $E >= 0. " 
        ))
    end
=#
    if E >= 0.0
        v_now = sqrt(v2_now)
        return (rs = fill(r, N_orb), vs = fill(v_now, N_orb))
    end

    Jc = last(find_Lc(E))
    #println("E = $E, J = $J, J/Jc = $(J / Jc)")
    
    # Radial velocity squared
    radial_vr2_raw(x) = 2.0 * (E - psi_calc(x)) - J^2 / x^2

    radial_vr2(x) = begin
        y = radial_vr2_raw(x)

        if !isfinite(y)
            return y
        end

        # Treat tiny numerical negatives/positives near a turning point as zero.
        if abs(y) < tol
            return 0.0
        end

        return y
    end

    vr2_now = radial_vr2(r)

    if !(isfinite(vr2_now)) || vr2_now < 0.0
        throw(ArgumentError(
            "Current point is not inside the allowed radial region: " *
            "r = $r, vr = $vr, vt = $vt, E = $E, J = $J, " *
            "radial_vr2(r) = $vr2_now."
        ))
    end

    # bisection root finder
    function bisect_root(f, a, b)
        fa = f(a)
        fb = f(b)

        if !isfinite(fa) || !isfinite(fb)
            throw(ArgumentError(
                "Non-finite function value in root bracket: " *
                "a = $a, f(a) = $fa, b = $b, f(b) = $fb."
            ))
        end

        if abs(fa) <= tol
            return a
        end

        if abs(fb) <= tol
            return b
        end

        if fa * fb > 0.0
            throw(ArgumentError(
                "Root is not bracketed: " *
                "a = $a, f(a) = $fa, b = $b, f(b) = $fb."
            ))
        end

        left, right = a, b
        f_left, f_right = fa, fb

        for _ in 1:200
            mid = 0.5 * (left + right)
            f_mid = f(mid)

            if !isfinite(f_mid)
                throw(ArgumentError(
                    "Non-finite function value during bisection: " *
                    "mid = $mid, f(mid) = $f_mid."
                ))
            end

            if abs(f_mid) <= tol || abs(right - left) <= tol * max(1.0, abs(mid))
                return mid
            end

            if f_left * f_mid > 0.0
                left = mid
                f_left = f_mid
            else
                right = mid
                f_right = f_mid
            end
        end

        return 0.5 * (left + right)
    end

    # Find pericenter rmin
    
    function find_inner_turning_point()
        x_hi = r
        f_hi = radial_vr2(x_hi)

        for _ in 1:max_expand
            x_lo = x_hi / expand_factor
            f_lo = radial_vr2(x_lo)

            if !isfinite(f_lo)
                throw(ArgumentError(
                    "Non-finite radial_vr2 while searching inward: " *
                    "x_lo = $x_lo, f_lo = $f_lo."
                ))
            end

            # We have moved from allowed region to forbidden region.
            if f_lo <= 0.0
                # If the current point itself is the inner turning point.
                if abs(f_hi) <= tol && x_hi == r
                    return r
                end

                return bisect_root(radial_vr2, x_lo, x_hi)
            end

            x_hi = x_lo
            f_hi = f_lo
        end

        throw(ArgumentError(
            "Could not find inner turning point for orbit: " *
            "r = $r, vr = $vr, vt = $vt, E = $E, J = $J."
        ))
    end

    # Find apocenter rmax
    
    function find_outer_turning_point()
        x_lo = r
        f_lo = radial_vr2(x_lo)

        for _ in 1:max_expand
            x_hi = x_lo * expand_factor
            f_hi = radial_vr2(x_hi)

            if !isfinite(f_hi)
                throw(ArgumentError(
                    "Non-finite radial_vr2 while searching outward: " *
                    "x_hi = $x_hi, f_hi = $f_hi."
                ))
            end

            # We have moved from allowed region to forbidden region.
            if f_hi <= 0.0
                # If the current point itself is the outer turning point.
                if abs(f_lo) <= tol && x_lo == r
                    return r
                end

                return bisect_root(radial_vr2, x_lo, x_hi)
            end

            x_lo = x_hi
            f_lo = f_hi
        end

        throw(ArgumentError(
            "Could not find outer turning point for orbit: " *
            "r = $r, vr = $vr, vt = $vt, E = $E, J = $J."
        ))
    end

    rmin = find_inner_turning_point()
    rmax = find_outer_turning_point()

    if !(isfinite(rmin) && isfinite(rmax) && rmin > 0.0 && rmax > 0.0)
        throw(ArgumentError(
            "Invalid turning points: rmin = $rmin, rmax = $rmax."
        ))
    end

    if rmax < rmin
        throw(ArgumentError(
            "Invalid orbit: rmax < rmin. rmin = $rmin, rmax = $rmax."
        ))
    end

    # Near-circular or degenerate orbit: just return the current point.
    if abs(rmax - rmin) <= tol * max(1.0, r)
        v_now = sqrt(v2_now)
        return (rs = fill(r, N_orb), vs = fill(v_now, N_orb))
    end


    # time-weighted radial CDF
    
    r_grid_full = collect(range(rmin, rmax; length = n_grid + 2))

    # Exclude exact turning points because vr = 0 there.
    r_grid = r_grid_full[2:end-1]

    vr2_grid = radial_vr2.(r_grid)

    for j in eachindex(vr2_grid)
        if !(isfinite(vr2_grid[j]) && vr2_grid[j] > 0.0)
            throw(ArgumentError(
                "Invalid radial velocity squared on sampling grid: " *
                "j = $j, r = $(r_grid[j]), vr2 = $(vr2_grid[j]), " *
                "rmin = $rmin, rmax = $rmax, E = $E, J = $J."
            ))
        end
    end

    weights = 1.0 ./ sqrt.(vr2_grid)

    # Cumulative trapezoidal integral.
    cdf = zeros(Float64, length(r_grid))

    for j in 2:length(r_grid)
        dr_local = r_grid[j] - r_grid[j-1]
        cdf[j] = cdf[j-1] + 0.5 * (weights[j-1] + weights[j]) * dr_local
    end

    total = cdf[end]

    if !(isfinite(total) && total > 0.0)
        throw(ArgumentError(
            "Invalid radial CDF normalization: total = $total."
        ))
    end

    cdf ./= total

    # Inverse-CDF sample radii
    
    rs = Vector{Float64}(undef, N_orb)

    for k in 1:N_orb
        u = rand(rng)

        idx = searchsortedfirst(cdf, u)

        if idx <= 1
            rs[k] = r_grid[1]

        elseif idx > length(r_grid)
            rs[k] = r_grid[end]

        else
            c0 = cdf[idx-1]
            c1 = cdf[idx]
            r0 = r_grid[idx-1]
            r1 = r_grid[idx]

            if c1 == c0
                rs[k] = r0
            else
                t = (u - c0) / (c1 - c0)
                rs[k] = r0 + t * (r1 - r0)
            end
        end
    end

    # v at r
    
    vs = Vector{Float64}(undef, N_orb)

    for k in 1:N_orb
        v2_sample = 2.0 * (E - psi_calc(rs[k]))

        if !(isfinite(v2_sample) && v2_sample > 0.0)
            throw(ArgumentError(
                "Invalid sampled total speed: " *
                "k = $k, rs = $(rs[k]), v2_sample = $v2_sample, " *
                "E = $E, J = $J."
            ))
        end

        vs[k] = sqrt(v2_sample)
    end

    return (rs = rs, vs = vs)
end



function orbit_sampled_kde_rv_df(
    dat::SnapDat,
    N_orb::Integer = 1;
    rng_seed::Union{Nothing, Integer} = nothing,
)

    if N_orb < 1
        throw(ArgumentError("N_orb must be at least 1."))
    end

    N_obj = length(dat.r)

    # Build potential function from current snapshot
    psi_tab, psi_tot_tab, M_tab = find_psi_arrays(dat.r, dat.m)
    psi_calc = psi_exact(dat.r, psi_tab, M_tab)

    # Number of bins 
    nbr = floor(Int64, (log(dat.rmax) - log(dat.rmin)) / dlogr)
    nbv = floor(Int64, (log(dat.vmax) - log(dat.vmin)) / dlogv)

    # Bin centers (log space)
    logr_samp = log(dat.rmin) .+ (collect(0:nbr-1) .+ 0.5) .* dlogr
    logv_samp = log(dat.vmin) .+ (collect(0:nbv-1) .+ 0.5) .* dlogv

    # Storage matrices 
    tabn_rv = zeros(Float64, nbr, nbv)
    tabf_rv = zeros(Float64, nbr, nbv)

    # Each orbit sample gets weight 1/N_orb
    w = 1.0 / N_orb


    orbit_counter = Threads.Atomic{Int}(0)
    t_start = time()
    
    Threads.@threads for i in 1:N_obj
    
        r    = dat.r[i]
        vr   = dat.vr[i]
        vt   = dat.vt[i]
        mass = dat.m[i]
    
        particle_seed = rng_seed === nothing ? nothing : rng_seed + 1000003 * i
    
        samples = orbit_sample_rv(
            r,
            vr,
            vt,
            N_orb,
            psi_calc;
            rng_seed = particle_seed,
        )
    
        done = Threads.atomic_add!(orbit_counter, 1) + 1
        if done % 10000 == 0
            elapsed = time() - t_start
            rate = done / elapsed
            eta = (N_obj - done) / rate
            println("Orbit sampled: $done / $N_obj, elapsed = $(round(elapsed/60, digits=2)) min, rate = $(round(rate, digits=2)) orbits/s, ETA = $(round(eta/60, digits=2)) min")
        end
    
        # Bin each orbit-sampled point
        for k in eachindex(samples.rs, samples.vs)

            rs = samples.rs[k]
            vs = samples.vs[k]

            # Orbit samples should always be valid before taking logs
            if !(isfinite(rs) && isfinite(vs) && rs > 0.0 && vs > 0.0)
                throw(ArgumentError(
                    "i = $i, sample k = $k: " *
                    "rs = $rs, vs = $vs. " 
                ))
            end
            
            # Locate sampled (r,v) bin in log space 
            ir = floor(Int64, (log(rs) - log(dat.rmin)) / dlogr) + 1
            iv = floor(Int64, (log(vs) - log(dat.vmin)) / dlogv) + 1

            # Only sum if sampled particle falls inside bin / grid 
            if (1 <= ir <= nbr) && (1 <= iv <= nbv)
                @inbounds tabn_rv[ir,iv] += mass * w / (dlogr * dlogv)
                @inbounds tabf_rv[ir,iv] += mass^2 * w / (dlogr * dlogv)
            end
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






# --- Distribution Function Sampling --- #

struct SampDF
    tabr_samp::Vector{Float64}
    tabv_samp::Vector{Float64}
    tabn_rv::Matrix{Float64}
    tabf_rv::Matrix{Float64}
    N::Function
    F::Function
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








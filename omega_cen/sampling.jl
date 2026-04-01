# --- Distribution Function Sampling --- #
# --- Working (Tested 1/24/2026) --- #

## Struct to store DF data ##
struct SampDF
    tabr_samp::Vector{Float64}
    tabv_samp::Vector{Float64}
    tabn_rv::Matrix{Float64}
    tabf_rv::Matrix{Float64}
    tabth_rv::Matrix{Float64}
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
    tabth_rv = zeros(Float64, nbr, nbv)

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

    Threads.@threads for index=1:nbr*nbv # Iterate over all bins

        # Find bin indices
        ir = floor(Int64, (index-1)/nbv) + 1
        iv = index - (ir-1)*nbv

        # Find r,v; theoretical DF
        r = tabr_samp[ir]
        v = tabv_samp[iv]
        df = plum_df_rv(r,v)

        # Save to array
        @inbounds tabth_rv[ir,iv] = df
    end

    # Interpolate DF
    N_samp = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabn_rv)
    F_samp = (r,v) -> bilinear_interp(r,v,tabr_samp,tabv_samp,tabf_rv)

    # Return DF (linear) struct
    return SampDF(tabr_samp,tabv_samp,tabn_rv,tabf_rv,tabth_rv,N_samp,F_samp)
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
    tabth_rv = zeros(Float64, nbr, nbv)

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

    Threads.@threads for index in 1:(nbr*nbv) # Iterate over all bins

        # Find bin indices
        ir = floor(Int64, (index - 1) / nbv) + 1
        iv = index - (ir - 1) * nbv

        # Find (r,v); DF
        r, v = exp(logr_samp[ir]), exp(logv_samp[iv])

        # Jacobian correction -> F(logr,logv) = r * v * F(r,v)
        df = r * v * plum_df_rv(r,v) # Eq 1.35

        # Save to array
        @inbounds tabth_rv[ir,iv] = df
    end

    N_samp = (r,v) -> bilinear_interp(r,v,logr_samp,logv_samp,tabn_rv)
    F_samp = (r,v) -> bilinear_interp(r,v,logr_samp,logv_samp,tabf_rv)

    # Return DF (log) struct
    return SampDF(logr_samp,logv_samp,tabn_rv,tabf_rv,tabth_rv,N_samp,F_samp)
end












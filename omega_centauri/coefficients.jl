# --- Diffusion coefficient calculations --- #

## --- Precomputed cumulative velocity-space integrals of the DF --- ##
#
# The four integrals in vel_coef are cumulative-in-logv integrals of the
# (fixed) DF at a given logr, evaluated at the object's speed vel(r):
#
#   calc1 = ∫_{vmin}^{vel} F   dlogv
#   calc2 = ∫_{vmin}^{vel} N   dlogv
#   calc3 = (1/vel^3) ∫_{vmin}^{vel} v^2 F dlogv
#   calc4 = ∫_{vel}^{vmax} (F/v) dlogv
#
# Since DF is bilinear on a fixed (logr, logv) grid, precompute the
# cumulative logv-integrals column-by-column (per logr node) ONCE. Then
# vel_coef reduces from O(n_steps) quadrature per call to O(1) lookups.
struct VelIntTables
    logr_nodes::Vector{Float64}
    logv_nodes::Vector{Float64}
    dlogv::Float64
    cou_log::Float64
    tabF::Matrix{Float64}       # DF.F samples this (squared mass)
    tabN::Matrix{Float64}       # DF.N samples this (mass)
    cumF::Matrix{Float64}       # ∫ F   dlogv  (to each logv node)
    cumN::Matrix{Float64}       # ∫ N   dlogv
    cumV2F::Matrix{Float64}     # ∫ v^2 F dlogv
    cumFinv::Matrix{Float64}    # ∫ F/v dlogv
end

# Populated by generate_coeffs before the grid loop; read (only) by vel_coef.
const VEL_INT_TABLES = Ref{Union{Nothing, VelIntTables}}(nothing)

## Build cumulative velocity-integral tables from a (fixed) DF ##
function build_vel_int_tables(df; sub::Integer = 1)

    logr_nodes = collect(df.tabr_samp)
    logv_nodes = collect(df.tabv_samp)
    tabF = df.tabf_rv           # bilinear-sampled by DF.F
    tabN = df.tabn_rv           # bilinear-sampled by DF.N

    n_r = length(logr_nodes)
    n_v = length(logv_nodes)

    if n_r < 2 || n_v < 2
        error("build_vel_int_tables: DF grid too small (n_r = $n_r, n_v = $n_v).")
    end

    dlv = logv_nodes[2] - logv_nodes[1]

    cumF    = zeros(Float64, n_r, n_v)
    cumN    = zeros(Float64, n_r, n_v)
    cumV2F  = zeros(Float64, n_r, n_v)
    cumFinv = zeros(Float64, n_r, n_v)

    # Trapezoidal cumulative integral along logv, at each logr node.
    # `sub` sub-samples each DF bin to tighten the exp-weighted (v^2 F, F/v)
    # integrals; sub = 1 integrates on the native DF nodes.
    Threads.@threads for i in 1:n_r
        @inbounds for k in 2:n_v
            F_lo = tabF[i, k-1]; F_hi = tabF[i, k]
            N_lo = tabN[i, k-1]; N_hi = tabN[i, k]

            lv_lo = logv_nodes[k-1]
            dv    = (logv_nodes[k] - lv_lo) / sub

            aF = aN = aV2F = aFinv = 0.0

            # Trapezoid over `sub` sub-intervals of bin [k-1, k].
            for s in 1:sub
                lvL = lv_lo + (s - 1) * dv
                lvR = lv_lo + s * dv
                tL = (s - 1) / sub
                tR = s / sub

                FL = F_lo + tL * (F_hi - F_lo)
                FR = F_lo + tR * (F_hi - F_lo)
                NL = N_lo + tL * (N_hi - N_lo)
                NR = N_lo + tR * (N_hi - N_lo)

                vL = exp(lvL); vR = exp(lvR)
                half = 0.5 * dv

                aF    += half * (FL + FR)
                aN    += half * (NL + NR)
                aV2F  += half * (vL^2 * FL + vR^2 * FR)
                aFinv += half * (FL / vL + FR / vR)
            end

            cumF[i, k]    = cumF[i, k-1]    + aF
            cumN[i, k]    = cumN[i, k-1]    + aN
            cumV2F[i, k]  = cumV2F[i, k-1]  + aV2F
            cumFinv[i, k] = cumFinv[i, k-1] + aFinv
        end
    end

    cou_log = log(M_bh / mean(dat.m))

    return VelIntTables(logr_nodes, logv_nodes, dlv, cou_log,
                        tabF, tabN, cumF, cumN, cumV2F, cumFinv)
end

## Cumulative integrals up to log(vel) for a single logr column ##
@inline function _col_cum_up(tb::VelIntTables, i::Int, log_vel::Float64)

    logv_nodes = tb.logv_nodes
    n_v = length(logv_nodes)

    lv1 = logv_nodes[1]
    lvN = logv_nodes[n_v]

    # Below the grid the DF is zero; above it, the integral is complete.
    if log_vel <= lv1
        return 0.0, 0.0, 0.0, 0.0
    elseif log_vel >= lvN
        @inbounds return tb.cumF[i, n_v], tb.cumN[i, n_v], tb.cumV2F[i, n_v], tb.cumFinv[i, n_v]
    end

    # Locate bin k with logv_nodes[k] <= log_vel < logv_nodes[k+1].
    k = floor(Int, (log_vel - lv1) / tb.dlogv) + 1
    k = clamp(k, 1, n_v - 1)

    @inbounds begin
        vk  = logv_nodes[k]
        dvp = log_vel - vk
        t   = dvp / tb.dlogv

        # DF columns are linear in logv between nodes.
        Fk = tb.tabF[i, k]; Fv = Fk + t * (tb.tabF[i, k+1] - Fk)
        Nk = tb.tabN[i, k]; Nv = Nk + t * (tb.tabN[i, k+1] - Nk)

        ev_k = exp(vk)
        ev_v = exp(log_vel)

        # Add the partial-bin trapezoid from node k up to log_vel.
        F_up    = tb.cumF[i, k]    + 0.5 * (Fk + Fv) * dvp
        N_up    = tb.cumN[i, k]    + 0.5 * (Nk + Nv) * dvp
        V2F_up  = tb.cumV2F[i, k]  + 0.5 * (ev_k^2 * Fk + ev_v^2 * Fv) * dvp
        Finv_up = tb.cumFinv[i, k] + 0.5 * (Fk / ev_k + Fv / ev_v) * dvp
    end

    return F_up, N_up, V2F_up, Finv_up
end

## Evaluate the four velocity integrals at (logr, vel) via table lookup ##
function vel_integrals(tb::VelIntTables, logr::Float64, vel::Float64)

    logr_nodes = tb.logr_nodes
    n_r = length(logr_nodes)
    n_v = length(tb.logv_nodes)

    # Outside the logr grid the DF is zero (matches bilinear_interp).
    if logr < logr_nodes[1] || logr >= logr_nodes[n_r]
        return 0.0, 0.0, 0.0, 0.0
    end

    ir = searchsortedlast(logr_nodes, logr)
    if ir < 1 || ir >= n_r
        return 0.0, 0.0, 0.0, 0.0
    end

    @inbounds begin
        r1 = logr_nodes[ir]; r2 = logr_nodes[ir+1]
        w2 = (logr - r1) / (r2 - r1)
        w1 = 1.0 - w2

        log_vel = log(vel)

        F1, N1, V1, Fi1 = _col_cum_up(tb, ir,     log_vel)
        F2, N2, V2, Fi2 = _col_cum_up(tb, ir + 1, log_vel)

        F_up    = w1 * F1  + w2 * F2
        N_up    = w1 * N1  + w2 * N2
        V2F_up  = w1 * V1  + w2 * V2
        Finv_up = w1 * Fi1 + w2 * Fi2

        # calc4 integrates from vel up to vmax = (full) - (up to vel).
        Finv_tot = w1 * tb.cumFinv[ir, n_v] + w2 * tb.cumFinv[ir + 1, n_v]
    end

    calc1 = F_up
    calc2 = N_up
    calc3 = V2F_up / vel^3
    calc4 = Finv_tot - Finv_up

    return calc1, calc2, calc3, calc4
end

## Calculate velocity diff coeffs using log sampling ##
function vel_coef(r, E, L)

    tb = VEL_INT_TABLES[]
    if tb === nothing
        error("vel_coef: velocity-integral tables not built. Call " *
              "generate_coeffs first (or set VEL_INT_TABLES[] = build_vel_int_tables(DF)).")
    end

    # Object total velocity
    vel = sqrt(find_vr(r, E, L)^2 + find_vt(r, E, L)^2)
    logr = log(r) # Convert into log space

    # Cumulative velocity-space integrals (precomputed -> O(1) lookups)
    calc1, calc2, calc3, calc4 = vel_integrals(tb, logr, vel)

    # Constants (cou_log precomputed once in the tables)
    c1 = (-G^2 * tb.cou_log) / (vel^2 * r^3)
    c2 = (2 * G^2 * tb.cou_log) / (3 * r^3)

    # Drift coefficient -> Split by mass dependence
    vp1 = c1 * calc1 # Eq 1.52
    vp2 = c1 * calc2 # Mass dependant

    # Diffusion coefficients -> No mass dependence
    vp_sq = c2 * (calc3 + calc4) # Eq 1.53
    vt_sq = c2 * ((3 / vel) * calc1 - calc3 + 2 * calc4) # Eq 1.54

    return vp1, vp2, vp_sq, vt_sq
end

## Reference (pre-optimization) velocity diff coeffs -- kept for A/B validation ##
function vel_coef_ref(r, E, L)

    # Object total velocity
    vel = sqrt(find_vr(r, E, L)^2 + find_vt(r, E, L)^2)

    logr = log(r)
    log_vel = log(vel)

    cou_log = log(M_bh / mean(dat.m))

    c1 = (-G^2 * cou_log) / (vel^2 * r^3)
    c2 = (2 * G^2 * cou_log) / (3 * r^3)

    int1(logv) = DF.F(logr, logv)
    int2(logv) = DF.N(logr, logv)
    int3(logv) = (v = exp(logv); (v^2 / vel^3) * DF.F(logr, logv))
    int4(logv) = (v = exp(logv); DF.F(logr, logv) / v)

    calc1 = midpoint(x -> int1(x), log(dat.vmin), log_vel, n_steps)
    calc2 = midpoint(x -> int2(x), log(dat.vmin), log_vel, n_steps)
    calc3 = midpoint(x -> int3(x), log(dat.vmin), log_vel, n_steps)
    calc4 = midpoint(x -> int4(x), log_vel, log(dat.vmax), n_steps)

    vp1 = c1 * calc1
    vp2 = c1 * calc2
    vp_sq = c2 * (calc3 + calc4)
    vt_sq = c2 * ((3 / vel) * calc1 - calc3 + 2 * calc4)

    return vp1, vp2, vp_sq, vt_sq
end

## Calculate local (E,L) diff coeffs ##
function loc_coef(r, E, L)

    if r < 0
        @warn "loc_coef received negative r; using abs(r) for calculation" E L r
        r = abs(r)
    elseif r == 0
        throw(DomainError(r, "loc_coef received r = 0, cannot take log(r)."))
    end

    # Object total velocity
    vel = sqrt(find_vr(r, E, L)^2 + find_vt(r, E, L)^2)

    logr = log(r) # Convert into log space
    log_vel = log(vel)

    # Vel diff coefs
    vp1, vp2, vp_sq, vt_sq = vel_coef(r, E, L)

    # 1st order coefs split, 2nd order unchanged
    delE1 = 0.5 * (vp_sq + vt_sq) + vel * vp1 # Eq 1.55
    delE2 = vel * vp2 # Mass dependent 

    delL1 = (L / vel) * vp1 + (r^2) / (4 * L) * vt_sq # Eq 1.57
    delL2 = (L / vel) * vp2 # Mass dependent

    delE_sq = vel^2 * vp_sq # Eq 1.56
    delL_sq = (L^2 / vel^2) * vp_sq + 0.5 * (r^2 - L^2 / vel^2) * vt_sq # Eq 1.58
    delEL = L * vp_sq # Eq 1.59

    return delE1, delE2, delL1, delL2, delE_sq, delL_sq, delEL
end


## Calculate orbit averaged (E,L) diff coeffs ##
function avg_coef_el(E,L)

    if E > 0
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end
    
    # Calculate orbit turning points, period
    rp, ra = find_rp_ra(E,L)
    T = find_period(E,L)
    h = pi/n_steps

    # Initialize variables 
    DE1 = DE2 = DL1 = DL2 = DEE = DLL = DEL = 0.0

    # Simultaneous midpoint integration over orbit #
    for i in 1:n_steps
        
        # Anomaly conversion
        tt = -pi/2 + (i - 0.5) * h
        r_hat = (ra - rp) / 2 # Eq 1.28
        t_hat = (ra + rp) / 2 # Eq 1.29
        rr = r_hat * sin(tt) + t_hat # Eq 1.27

        pref = (2h * r_hat * cos(tt)) / (T * find_vr(rr,E,L)) # Prefactor

        # Compute local diff coeffs 
        dE1, dE2, dL1, dL2, dEE, dLL, dEL = loc_coef(rr,E,L)

        DE1 += pref * dE1 # Eq 1.65
        DE2 += pref * dE2 # Mass dependent

        DL1 += pref * dL1 # Eq 1.67
        DL2 += pref * dL2 # Mass dependent 

        DEE += pref * dEE # Eq 1.66
        DLL += pref * dLL # Eq 1.68
        DEL += pref * dEL # Eq 1.69
    end

    return DE1, DE2, DL1, DL2, DEE, DLL, DEL
end

## Compute orbit averaged coeffs in (E,j) ##
function avg_coef_ej(E, j)
    # Find Jc, derivatives; J
    Jc = find_Lc(E)[2]
    dJc = dLc_dE(E)
    d2Jc = d2Lc_dE2(E)
    J = j * Jc

    # Compute avg diff coeffs in (E,L)
    DE1_L, DE2_L, DL1_L, DL2_L, DEE_L, DLL_L, DEL_L = avg_coef_el(E, J)

    # Energy coeffs preserved under transformation
    DE1, DE2, DEE = DE1_L, DE2_L, DEE_L

    # Convert j coeffs -> Eq 1.73, 1.74
    Dj1_1 = (DL1_L/Jc) - (J/(Jc^2))*dJc*DE1_L - (1/(Jc^2))*dJc*DEL_L
    Dj1_2 = -(J/2)*((d2Jc/(Jc^2)) - (2*(dJc^2))/(Jc^3))*DEE_L

    Dj1   = Dj1_1 + Dj1_2
    Dj2   = (DL2_L/Jc) - (J/(Jc^2))*dJc*DE2_L # Mass dependent 

    # Diffusion -> No Mass dependence
    Djj   = ((J*dJc)/(Jc^2))^2 * DEE_L - 2*(J*dJc/(Jc^3))*DEL_L + (DLL_L/(Jc^2)) # Eq 1.76
    DEj   = -(J/(Jc^2))*dJc*DEE_L + (DEL_L/Jc) # Eq 1.75

    return DE1, DE2, Dj1, Dj2, DEE, Djj, DEj
end

## Compute diff coeff grid (E,j) for sampling ##
function coef_grid(E_tab,j_tab,res)

    # Initialize matrices
    DE1_tab = zeros(res,res)
    DE2_tab = zeros(res,res)
    Dj1_tab = zeros(res,res)
    Dj2_tab = zeros(res,res)
    DEE_tab = zeros(res,res)
    Djj_tab = zeros(res,res)
    DEj_tab = zeros(res,res)

    # Loop over energy, momentum
    Threads.@threads for i in 1:res
        E = E_tab[i]
        for k in 1:res
            j = j_tab[k]
            DE1, DE2, Dj1, Dj2, DEE, Djj, DEj = avg_coef_ej(E,j)
            @inbounds begin
                DE1_tab[k,i] = DE1
                DE2_tab[k,i] = DE2
                Dj1_tab[k,i] = Dj1
                Dj2_tab[k,i] = Dj2
                DEE_tab[k,i] = DEE
                Djj_tab[k,i] = Djj
                DEj_tab[k,i] = DEj
            end
        end
    end
    return DE1_tab, DE2_tab, Dj1_tab, Dj2_tab, DEE_tab, Djj_tab, DEj_tab
end


# ======================================================================================
# DIAGNOSTIC: fully decompose the first-order J drift Dj1 into the INDIVIDUAL terms of its
# (E,L)->(E,j) transform (Eq 1.73/1.74):
#   Dj1_1 = Dj1_1_1 + Dj1_1_2 + Dj1_1_3   (three terms)
#   Dj1_2 = Dj1_2_1 + Dj1_2_2             (two terms)
#   Dj1   = Dj1_1 + Dj1_2
# Lc(E) derivatives use finite differences (the cubic-spline comparison showed the spline
# makes no difference, so it is dropped here). Not used by production (generate_coeffs);
# driven only by run_dj1_decomp.jl.
# ======================================================================================
function avg_coef_ej_decomp(E, j)
    Jc = find_Lc(E)[2]
    J  = j * Jc

    # Raw orbit-averaged (E,L) coefficients (the expensive part).
    DE1_L, DE2_L, DL1_L, DL2_L, DEE_L, DLL_L, DEL_L = avg_coef_el(E, J)

    # Lc(E) derivatives (finite difference).
    dJc  = dLc_dE_fd(E)
    d2Jc = d2Lc_dE2_fd(E)

    # Dj1_1 (Eq 1.73), three terms:
    Dj1_1_1 = DL1_L / Jc                                 # angular-momentum drift, DL1/Jc
    Dj1_1_2 = -(J/(Jc^2)) * dJc * DE1_L                  # E-drift through the moving ceiling Lc(E)
    Dj1_1_3 = -(1/(Jc^2)) * dJc * DEL_L                  # E-L cross diffusion via the j tilt
    Dj1_1   = Dj1_1_1 + Dj1_1_2 + Dj1_1_3

    # Dj1_2 (Eq 1.74), two terms (the Ito/curvature drift, split by Lc-derivative order):
    Dj1_2_1 = -(J/2) * (d2Jc/(Jc^2)) * DEE_L             # d2Lc/dE2 term
    Dj1_2_2 = -(J/2) * (-(2*(dJc^2))/(Jc^3)) * DEE_L     # (dLc/dE)^2 term = +(J*dJc^2/Jc^3)*DEE_L
    Dj1_2   = Dj1_2_1 + Dj1_2_2

    Dj1 = Dj1_1 + Dj1_2

    return (Jc = Jc, J = J,
            DE1_L = DE1_L, DL1_L = DL1_L, DEE_L = DEE_L, DEL_L = DEL_L,
            dJc = dJc, d2Jc = d2Jc,
            Dj1_1_1 = Dj1_1_1, Dj1_1_2 = Dj1_1_2, Dj1_1_3 = Dj1_1_3, Dj1_1 = Dj1_1,
            Dj1_2_1 = Dj1_2_1, Dj1_2_2 = Dj1_2_2, Dj1_2 = Dj1_2,
            Dj1 = Dj1)
end

# Fields emitted per grid point (matrices are [j, E] = [k, i], matching coef_grid).
const DJ1_DECOMP_FIELDS = (:Jc, :J, :DE1_L, :DL1_L, :DEE_L, :DEL_L, :dJc, :d2Jc,
                           :Dj1_1_1, :Dj1_1_2, :Dj1_1_3, :Dj1_1,
                           :Dj1_2_1, :Dj1_2_2, :Dj1_2, :Dj1)

function coef_grid_decomp(E_tab, j_tab, res)
    out = Dict{Symbol,Matrix{Float64}}(f => zeros(res, res) for f in DJ1_DECOMP_FIELDS)
    done = Threads.Atomic{Int}(0)
    Threads.@threads for i in 1:res            # threads write disjoint E-columns -> no race
        E = E_tab[i]
        for k in 1:res
            r = avg_coef_ej_decomp(E, j_tab[k])
            for f in DJ1_DECOMP_FIELDS
                @inbounds out[f][k, i] = getproperty(r, f)
            end
        end
        d = Threads.atomic_add!(done, 1) + 1
        d % 10 == 0 && (println("  coef_grid_decomp: $d / $res E-columns done"); flush(stdout))
    end
    return out
end

# Set up the SAME grid as generate_coeffs (from the active DF/potential globals) and run the
# term decomposition. Returns (E_tab, j_tab, grids).
function generate_dj1_decomp(res)
    VEL_INT_TABLES[] = build_vel_int_tables(DF)
    EE = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
    emin, emax = extrema(EE[(EE .< 0) .& (dat.startype .> -100)])

    j_min = 1e-6; j_max = 0.999
    j_edges = 10 .^ range(log10(j_min), log10(j_max), length = res + 1)
    j_tab   = sqrt.(j_edges[1:end-1] .* j_edges[2:end])
    E_abs_edges = 10 .^ range(log10(abs(emax)), log10(abs(emin)), length = res + 1)
    E_tab   = reverse(-sqrt.(E_abs_edges[1:end-1] .* E_abs_edges[2:end]))

    println("generate_dj1_decomp: res = $res, E in [$(E_tab[1]), $(E_tab[end])], " *
            "j in [$(j_tab[1]), $(j_tab[end])]")
    flush(stdout)
    grids = coef_grid_decomp(E_tab, j_tab, res)
    return (E_tab = E_tab, j_tab = j_tab, emin = emin, emax = emax, grids = grids)
end


struct DiffusionCoeffs

    # Orbit-averaged diff coeff grids
    DE1_tab::Matrix{Float64}
    DE2_tab::Matrix{Float64} # Mass dependent 
    Dj1_tab::Matrix{Float64}
    Dj2_tab::Matrix{Float64} # Mass dependent
    DEE_tab::Matrix{Float64}
    Djj_tab::Matrix{Float64}
    DEj_tab::Matrix{Float64}

    # Calculated energy and bounds (quantile)
    E_tab::AbstractVector{Float64}
    j_tab::AbstractVector{Float64}
    emin::Float64
    emax::Float64

    # Interpolated functions 
    DE1_samp::Function
    DE2_samp::Function # Mass dependent 
    Dj1_samp::Function
    Dj2_samp::Function # Mass dependent
    DEE_samp::Function
    Djj_samp::Function
    DEj_samp::Function
end

## Generate diff coeffs, interpolated functions ##
function generate_coeffs(res)

    # Precompute cumulative velocity-space integrals of the DF once,
    # so vel_coef does O(1) lookups instead of O(n_steps) quadrature.
    VEL_INT_TABLES[] = build_vel_int_tables(DF)

    # Calculate energy (& limits) for snapshot
    EE = 0.5 .* (dat.vr.^2 .+ dat.vt.^2) .+ psi_calc.(dat.r)
    #emin, emax = extrema(EE[(EE .< 0) .& (dat.startype .== 14)])
    emin, emax = extrema(EE[(EE .< 0) .& (dat.startype .> -100)])

    j_min = 1e-6
    j_max = 0.999
    
    j_edges = 10 .^ range(log10(j_min), log10(j_max), length=res + 1)
    
    # Geometric centers of the log bins
    j_tab = sqrt.(j_edges[1:end-1] .* j_edges[2:end])
    
    E_abs_edges = 10 .^ range(log10(abs(emax)), log10(abs(emin)), length=res + 1)
    
    E_edges = reverse(-E_abs_edges)
    E_tab = reverse(-sqrt.(E_abs_edges[1:end-1] .* E_abs_edges[2:end]))
        
    (DE1_tab, DE2_tab, Dj1_tab, Dj2_tab, 
         DEE_tab, Djj_tab, DEj_tab) = coef_grid(E_tab,j_tab,res)

    # Bilinear interpolation
    interp(Z) = (j,E) -> bilinear_interp_clamped(j, E, j_tab, E_tab, Z)
    
    return DiffusionCoeffs(
        DE1_tab, DE2_tab, Dj1_tab, Dj2_tab,
        DEE_tab, Djj_tab, DEj_tab,
        E_tab, j_tab, emin, emax,
        interp(DE1_tab), interp(DE2_tab),
        interp(Dj1_tab), interp(Dj2_tab),
        interp(DEE_tab), interp(Djj_tab),
        interp(DEj_tab))
end

## Save orbit-averaged diff coeff struct to HDF5 ##
function save_coeffs(name::String, dc::DiffusionCoeffs)
    h5open(name, "w") do file
        file["E_tab"] = collect(dc.E_tab)
        file["j_tab"] = collect(dc.j_tab)
        file["emin"] = dc.emin
        file["emax"] = dc.emax
        file["DE1_tab"] = dc.DE1_tab
        file["DE2_tab"] = dc.DE2_tab
        file["Dj1_tab"] = dc.Dj1_tab
        file["Dj2_tab"] = dc.Dj2_tab
        file["DEE_tab"] = dc.DEE_tab
        file["Djj_tab"] = dc.Djj_tab
        file["DEj_tab"] = dc.DEj_tab
    end
end

## Load diff coeff HDF5 ##
function load_coeffs(name::String)
    h5open(name, "r") do file
        # Read in coeff data
        E_tab = read(file["E_tab"])
        j_tab = read(file["j_tab"])
        emin = read(file["emin"])
        emax = read(file["emax"])
        DE1_tab = read(file["DE1_tab"])
        DE2_tab = read(file["DE2_tab"])
        Dj1_tab = read(file["Dj1_tab"])
        Dj2_tab = read(file["Dj2_tab"])
        DEE_tab = read(file["DEE_tab"])
        Djj_tab = read(file["Djj_tab"])
        DEj_tab = read(file["DEj_tab"])

        # Bilinear interpolation
        interp(Z) = (j,E) -> bilinear_interp_clamped(j, E, j_tab, E_tab, Z)

        return DiffusionCoeffs(
            DE1_tab, DE2_tab, Dj1_tab, Dj2_tab,
            DEE_tab, Djj_tab, DEj_tab,
            E_tab, j_tab, emin, emax,
            interp(DE1_tab), interp(DE2_tab),
            interp(Dj1_tab), interp(Dj2_tab),
            interp(DEE_tab), interp(Djj_tab),
            interp(DEj_tab))
    end
end

## Rebuild a DiffusionCoeffs' sampler closures with a chosen out-of-grid behavior, WITHOUT
## changing the grid. Used to select the MC extrapolation method for a query outside the
## coefficient grid (the physical extrapolation rule is unknown, so we bracket it). Note the
## J-axis ALWAYS clamps (a very small / near-radial j uses the nearest j-edge value, so those
## walkers keep diffusing toward the loss cone); only the E-axis differs by method:
##   boundary = :zero  -> 0 when E is outside the grid, clamp in j         (method 1)
##   boundary = :clamp -> nearest grid-edge value in both E and j          (method 3;
##                        this is also the generate_coeffs / load_coeffs default)
## Method 2 instead EXTENDS the grid downward via extrapolate_coeffs (power law), so it is
## not handled here.
function set_extrap_boundary(dc::DiffusionCoeffs, boundary::Symbol)
    boundary in (:zero, :clamp) || error("set_extrap_boundary: boundary must be :zero or :clamp")
    itp = boundary === :zero ? bilinear_interp_clampj_zeroE : bilinear_interp_clamped
    interp(Z) = (j, E) -> itp(j, E, dc.j_tab, dc.E_tab, Z)
    return DiffusionCoeffs(
        dc.DE1_tab, dc.DE2_tab, dc.Dj1_tab, dc.Dj2_tab,
        dc.DEE_tab, dc.Djj_tab, dc.DEj_tab,
        dc.E_tab, dc.j_tab, dc.emin, dc.emax,
        interp(dc.DE1_tab), interp(dc.DE2_tab),
        interp(dc.Dj1_tab), interp(dc.Dj2_tab),
        interp(dc.DEE_tab), interp(dc.Djj_tab),
        interp(dc.DEj_tab))
end

# ============================================================================
# Analytic ground-truth diffusion coefficients for the Keplerian Bahcall-Wolf cusp.
#
# Injects the SMOOTH analytic DF f(eps) ∝ eps^(gamma-3/2) and a PURE Keplerian potential
# (psi = -G M_bh/r) into the EXISTING coefficient machinery (build_vel_int_tables + coef_grid),
# giving noise-free coefficients on an (E,j) grid that can be extended below the lowest sampled
# particle. This is the ground truth against which the sampled+extrapolated coefficients are
# compared. Uses the untruncated cusp DF (the closed form is the infinite cusp); the r-grid is
# extended past [r_min, r_infl] so orbits that dip below r_min still see field stars.
#
# Requires main.jl (coefficient/potential/orbit fns) + analytic_sampling.jl (AnalyticCusp).
# ============================================================================

using Printf

## Beta function B(a,b) = ∫_0^1 x^(a-1)(1-x)^(b-1) dx (midpoint; avoids a SpecialFunctions dep). ##
function _beta(a::Real, b::Real; n::Integer = 40000)
    s = 0.0; h = 1.0 / n
    @inbounds for i in 1:n
        x = (i - 0.5) * h
        s += x^(a - 1) * (1 - x)^(b - 1)
    end
    return s * h
end

## Pure Keplerian potential in the pipeline globals. Zero stellar mass -> find_psi_arrays gives
## psi_stellar = 0 and M_enc = 0, so psi_tot = -G M_bh/r everywhere and find_Lc/find_rp_ra/
## find_period reproduce Kepler. `M_bh` (code) MUST already be set. A fine log-r grid covers the
## cusp with margin so tabulated look-ups are accurate. Activates it as the orbit potential.
function setup_kepler_potential!(model::AnalyticCusp; npts::Integer = 6000, margin::Real = 100.0)
    r_lo   = (model.r_min  / model.lengthunitpc) / margin      # code length, below r_min
    r_hi   = (model.r_infl / model.lengthunitpc) * 5.0
    r_grid = 10.0 .^ range(log10(r_lo), log10(r_hi), length = npts)
    psi_tab, psi_tot, M_tab = find_psi_arrays(r_grid, zeros(length(r_grid)))   # uses global M_bh
    psi_c = psi_exact(r_grid, psi_tab, M_tab)
    ORBIT_PSI[] = (r = r_grid, psi_tab = psi_tab, psi_tot = psi_tot, M_tab = M_tab, psi_calc = psi_c)
    activate_orbit_psi!()
    return nothing
end

## Analytic DF as a SampDF (grid of N = mass density, F = mass^2 density in (log r, log v)).
## N_log = K r^3 v^3 eps^p, F_log = m_star * N_log, with eps = M_bh/r - v^2/2 (code, G=1) and
## K = (2 sqrt(2) pi / B(p+1,3/2)) * A_code * M_bh^(-gamma). The r-grid is extended below r_min
## (untruncated cusp) so deep orbits see field stars; v-grid reaches v_esc at the inner edge.
function analytic_df(model::AnalyticCusp; nr::Integer = 600, nv::Integer = 600, margin::Real = 100.0)
    lpc = model.lengthunitpc; mu = model.massunitmsun
    M_bh_c = model.M_bh / mu
    p  = model.gamma - 1.5
    B  = _beta(p + 1, 1.5)
    A_code = model.rho_coeff * lpc^(3.0 - model.gamma) / mu           # rho = A_code r^-gamma (code)
    K  = (2 * sqrt(2) * pi / B) * A_code * M_bh_c^(-model.gamma)

    r_lo = (model.r_min / lpc) / margin
    r_hi = (model.r_infl / lpc) * 3.0
    v_hi = sqrt(2 * M_bh_c / r_lo)                                    # ~ v_esc at inner edge
    v_lo = v_hi / 1.0e5

    logr = collect(range(log(r_lo), log(r_hi), length = nr))
    logv = collect(range(log(v_lo), log(v_hi), length = nv))
    m_c  = model.m_star / mu
    tabn = zeros(nr, nv); tabf = zeros(nr, nv)
    @inbounds for i in 1:nr
        r = exp(logr[i])
        for k in 1:nv
            v = exp(logv[k])
            eps = M_bh_c / r - 0.5 * v^2
            if eps > 0.0
                Nl = K * r^3 * v^3 * eps^p
                tabn[i, k] = Nl
                tabf[i, k] = m_c * Nl
            end
        end
    end
    Nfun = (lr, lv) -> bilinear_interp(lr, lv, logr, logv, tabn)
    Ffun = (lr, lv) -> bilinear_interp(lr, lv, logr, logv, tabf)
    return SampDF(logr, logv, tabn, tabf, Nfun, Ffun)
end

## Density recovered from the analytic DF grid: rho(r) = ∫ N_log dlogv / (4 pi r^3). Used to
## VALIDATE the normalization against A_code r^-gamma.
function df_density(df::SampDF, r::Real)
    lr = log(r)
    lv = df.tabv_samp
    s = 0.0
    @inbounds for k in 2:length(lv)
        s += 0.5 * (df.N(lr, lv[k-1]) + df.N(lr, lv[k])) * (lv[k] - lv[k-1])
    end
    return s / (4 * pi * r^3)
end

## Build the analytic ground-truth coefficients on a deep (E,j) grid and save as a
## DiffusionCoeffs .h5 (same format as the sampled ones). E spans [E_deep, E_shallow] with
## E_deep = circular energy at r_min (deeper than the lowest sampled particle) and E_shallow
## from the sampled emax; j is the usual log grid. Returns the DiffusionCoeffs.
function generate_analytic_coeffs(model::AnalyticCusp, snapshot_file::AbstractString,
                                  conv_path::AbstractString, out_path::AbstractString;
                                  res::Integer = 200, r_deep_pc::Real = model.r_min)
    global M_bh, dat, DF

    set_model_constants!(conv_path)                    # units (c, r_conv, massunitmsun, ...)
    M_bh = model.M_bh / massunitmsun                    # code units

    key = snapshot_keys_by_time(snapshot_file; sep = model.r_infl >= 0 ? 0.1 : 0.1)[1]
    dat = read_file_h5(snapshot_file, key)              # for cou_log = log(M_bh/mean(dat.m))

    # Potential + DF grids must reach below the deepest orbit -> widen to cover r_deep_pc with a
    # decade of margin (df_margin = 100 when r_deep_pc = r_min, larger when going deeper).
    df_margin = max(100.0, 10.0 * model.r_min / r_deep_pc)
    setup_kepler_potential!(model; margin = df_margin)  # pure-Kepler potential globals
    df = analytic_df(model; margin = df_margin)
    DF = df
    VEL_INT_TABLES[] = build_vel_int_tables(df)

    # (E, j) grid: shallow end from the sampled particles; deep end = circular energy at r_deep_pc,
    # DECOUPLED from the sampling floor r_min (set r_deep_pc < r_min to extend deeper).
    EE = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
    emax_s = maximum(EE[EE .< 0])                       # least-bound sampled particle
    r_deep_code = r_deep_pc / model.lengthunitpc
    E_deep = -M_bh / (2 * r_deep_code)                  # circular energy at r_deep (deep bound)

    j_min = 1.0e-6; j_max = 0.999
    j_edges = 10.0 .^ range(log10(j_min), log10(j_max), length = res + 1)
    j_tab   = sqrt.(j_edges[1:end-1] .* j_edges[2:end])
    E_abs_edges = 10.0 .^ range(log10(abs(emax_s)), log10(abs(E_deep)), length = res + 1)
    E_tab   = reverse(-sqrt.(E_abs_edges[1:end-1] .* E_abs_edges[2:end]))

    println("generate_analytic_coeffs: res = $res, E in [$(E_tab[1]), $(E_tab[end])] " *
            "(sampled emax = $emax_s, r_deep = $r_deep_pc pc -> E_deep = $E_deep), j in [$(j_tab[1]), $(j_tab[end])]")
    flush(stdout)

    DE1, DE2, Dj1, Dj2, DEE, Djj, DEj = coef_grid(E_tab, j_tab, res)

    interp(Z) = (j, E) -> bilinear_interp_clamped(j, E, j_tab, E_tab, Z)
    dc = DiffusionCoeffs(DE1, DE2, Dj1, Dj2, DEE, Djj, DEj,
                         E_tab, j_tab, E_tab[1], E_tab[end],
                         interp(DE1), interp(DE2), interp(Dj1), interp(Dj2),
                         interp(DEE), interp(Djj), interp(DEj))
    save_coeffs(out_path, dc)
    println("generate_analytic_coeffs: wrote $out_path")
    flush(stdout)
    return dc
end


## Monte Carlo on the analytic model: seed one walker per BH from the snapshot, evolve in the
## ANALYTIC coefficients + pure-Keplerian potential (the same psi the coefficients were built in),
## and integrate for `max_t` Myr. Same output layout as automate_mc (per-snapshot group with the
## final states + a traj/<k> subgroup, including Egw, for each loss-cone capture) so the existing
## analysis reads it unchanged. No extrapolation method: the analytic grid already reaches r_g, and
## load_coeffs clamps at the edges. Heavy (threaded walker loop) -> submit as a job.
function automate_analytic_mc(snapshot_file::AbstractString, conv_path::AbstractString,
                              coeff_file::AbstractString, out_file::AbstractString;
                              M_bh_msun::Real, max_t::Real = 100.0,
                              max_steps::Integer = 100000, F_safe::Real = 10,
                              compute_rp_ra::Bool = false)
    global M_bh, dat

    set_model_constants!(conv_path)                      # units (c, r_conv, massunitmsun, ...)
    M_bh = M_bh_msun / massunitmsun                       # code units (IMBH not stored in the snapshot)
    set_F_safe!(F_safe)

    key   = snapshot_keys_by_time(snapshot_file; sep = 0.1)[1]
    t_gyr = _snapshot_time(key); tag = _time_tag(t_gyr)
    dat   = read_file_h5(snapshot_file, key)

    # Pure-Keplerian potential (matches the analytic coefficients): zero stellar mass -> psi = -M_bh/r.
    # The grid brackets the orbits, from r_g/10 (deep, near the loss cone) to 5x the outermost particle.
    r_grid = 10.0 .^ range(log10((M_bh / c^2) / 10), log10(5 * dat.rmax), length = 6000)
    kpsi, ktot, kM = find_psi_arrays(r_grid, zeros(length(r_grid)))
    ORBIT_PSI[] = (r = r_grid, psi_tab = kpsi, psi_tot = ktot, M_tab = kM,
                   psi_calc = psi_exact(r_grid, kpsi, kM))
    activate_orbit_psi!()

    coef = load_coeffs(coeff_file)                        # analytic coefficients (clamp at grid edges)
    ics  = black_hole_ics(dat)
    N    = length(ics.indices)
    println("automate_analytic_mc: N = $N walkers, max_t = $max_t Myr, F_safe = $F_safe, " *
            "M_bh = $M_bh (code), coeff = $(basename(coeff_file))")
    flush(stdout)

    Ef  = Vector{Float64}(undef, N); jf  = Vector{Float64}(undef, N)
    rpf = Vector{Float64}(undef, N); raf = Vector{Float64}(undef, N)
    reasons = Vector{Int}(undef, N)
    traj_t  = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)
    traj_E  = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)
    traj_j  = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)
    traj_rp = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)
    traj_ra = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)
    traj_Egw = Vector{Union{Nothing,Vector{Float64}}}(nothing, N)

    Threads.@threads :dynamic for k in 1:N
        ts, es, js, rps, ras, reason, egws =
            run_mc_rp_ra(ics.E0[k], ics.J0[k], coef, ics.m0[k];
                         max_step = max_steps, max_t = max_t, compute_rp_ra = compute_rp_ra)
        Ef[k] = es[end]; jf[k] = js[end]; rpf[k] = rps[end]; raf[k] = ras[end]
        reasons[k] = reason
        if reason == 1
            traj_t[k]  = ts;  traj_E[k]  = es;  traj_j[k] = js
            traj_rp[k] = rps; traj_ra[k] = ras; traj_Egw[k] = egws
        end
    end

    h5open(out_file, "w") do of
        attrs(of)["snapshot_file"] = snapshot_file
        attrs(of)["coeff_file"]    = coeff_file
        attrs(of)["max_t"]         = float(max_t)
        attrs(of)["max_steps"]     = max_steps
        attrs(of)["F_safe"]        = float(F_safe)
        attrs(of)["M_bh"]          = M_bh
        attrs(of)["nthreads"]      = Threads.nthreads()

        g = create_group(of, "$(tag)Gyr")
        attrs(g)["time"] = t_gyr; attrs(g)["M_bh"] = M_bh; attrs(g)["N"] = N
        g["indices"] = ics.indices
        g["E0s"] = ics.E0; g["J0s"] = ics.J0; g["m0s"] = ics.m0
        g["Ef"]  = Ef;     g["jf"]  = jf
        g["rpf"] = rpf;    g["raf"] = raf
        g["reason_to_end"] = reasons

        tg = create_group(g, "traj"); ncap = 0
        for k in 1:N
            traj_t[k] === nothing && continue
            ncap += 1
            wk = create_group(tg, string(k))
            attrs(wk)["orig_index"] = ics.indices[k]
            attrs(wk)["m"]  = ics.m0[k]
            attrs(wk)["E0"] = ics.E0[k]
            attrs(wk)["J0"] = ics.J0[k]
            wk["t"]   = traj_t[k];  wk["E"]  = traj_E[k];  wk["j"] = traj_j[k]
            wk["rp"]  = traj_rp[k]; wk["ra"] = traj_ra[k]; wk["Egw"] = traj_Egw[k]
        end
        attrs(tg)["n_captured"] = ncap
        println("    $N walkers integrated; $ncap captured (reason 1)")
        flush(stdout)
    end

    println("automate_analytic_mc: complete -> $out_file")
    flush(stdout)
    return out_file
end

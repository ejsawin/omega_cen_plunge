# --- Monte Carlo Functions --- #

function mc_step(E,j,coef,m_obj)

    # Characterize orbit 
    rc, Lc = find_Lc(E)
    L = j * Lc

    rp, ra = find_rp_ra(E, L)
    
    if orbit_failed(rp, ra)
        # Return whatever your code uses to mark a terminated/rejected orbit.
        # Do not calculate a, e, period, or GW rates.
        return -100.0, -100.0, 1000, 1000, false, 0.0
    end
    
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    
    # Extra protection against bad elements.
    if !isfinite(a) || !isfinite(e) || a <= 0 || e < 0 || e >= 1
        return -100.0, -100.0, 1000, 1000, false, 0.0
    end
    
    period = find_period(E, L, rp, ra)   # reuse rp, ra from above (skips a 2nd find_rp_ra)

    if !isfinite(period) || period <= 0
        return -100.0, -100.0, 1000, 1000, false, 0.0
    end
    # (period, a, e are already set above -- the duplicate find_period/a/e recompute was
    #  removed; find_period does root-finding/integration, so this is one fewer per step.)

    # Orbit averaged diff coeffs (per unit time!)
    dE_NR = coef.DE1_samp(j,E) + m_obj * coef.DE2_samp(j,E)
    dj_NR = coef.Dj1_samp(j,E) + m_obj * coef.Dj2_samp(j,E)
    dEE_NR = coef.DEE_samp(j,E)
    djj_NR = coef.Djj_samp(j,E)
    dEj_NR = coef.DEj_samp(j,E)

    # GW Emission (per unit time)
    dE_GW, dL_GW = gw_rates(a, e, m_obj)
    dj_GW = dL_GW / Lc
    
    # Combined NR + GW
    dE_tot = dE_NR + dE_GW
    dj_tot = dj_NR + dj_GW

    # Step Size. F_safe is a typed global (default 10, settable via automate_mc for
    # convergence tests); read here instead of a hardcoded local. N = max(1, N_safe) orbits.

    # Paper's two timescales: GW on E, and L-diffusion (in j here). t_jj is clamped
    # >= 0 so noisy interpolation can't give a negative timescale (-> Inf = unconstrained).
    t_gwE = abs(E / dE_GW)           # GW on E
    t_jj  = j^2 / max(djj_NR, 0.0)   # L-diffusion (in j)

    # Ignore degenerate L-diffusion timescales (t_jj = 0, Inf, or NaN -- e.g. djj_NR <= 0
    # off-grid, or j -> 0). The diffusion step is ill-defined there, so bail with the
    # sentinel; the furball / driver then drops this sample rather than guessing a step.
    # EXCEPTION (GW_ONLY, no-diffusion mode): djj_NR = 0 everywhere is intentional, so a
    # t_jj = Inf means "L-diffusion does not constrain the step" (GW sets it) rather than a
    # degenerate off-grid sample -- keep the walker instead of dropping it.
    if (!isfinite(t_jj) || t_jj <= 0.0) && !(GW_ONLY && isinf(t_jj))
        return -100.0, -100.0, 1000, 1000, false, 0.0
    end

    # Single step limiter (NO plunge branch): GW (t_gwE) and L-diffusion (t_jj) both cap
    # it. At least one full orbit, no upper cap -- the diffusion coeffs are orbit-averaged,
    # so a sub-orbital step would misuse them. When t_jj = Inf (GW_ONLY) L-diffusion drops
    # out of the min and the GW timescale alone sets the step.
    plunge = t_gwE < period   # unused here; kept only for the commented-out pinhole block
    jj_lim = isfinite(t_jj) ? sqrt(t_jj / period) : Inf
    N_safe = (1/F_safe) * min(t_gwE / period, jj_lim)
    dt = isfinite(N_safe) ? max(N_safe, 1.0) * period : period

    # Deterministic GW energy drained this step: the GW E-rate (dE_GW, per unit time and
    # specific in m_obj) times the actual step dt. Unlike the diffusion piece there is no
    # stochastic term -- it is exactly dE_GW*dt regardless of the random kicks, and dt is
    # the same dt the state update uses. The driver accumulates this into a per-step,
    # cumulative "GW energy dissipated" trajectory (dE_GW < 0, so each contribution is
    # negative; |cumsum| is the specific energy radiated -- no mass division, per request).
    dEgw_step = dE_GW * dt

    denom = sqrt(abs(dEE_NR * djj_NR))

    g1, g2 = randn(), randn()

    if denom > 0 # Diffusion present
        rho = clamp(dEj_NR / denom, -1.0, 1.0) # Ensure safe rho
        g1t = g1
        g2t = g1 * rho + g2 * sqrt(max(0.0, 1 - rho^2))
    else # No diffusion
        rho = 0.0
        g1t = g1
        g2t = g2
    end

    # Update (E,L)
    #dE = dE_tot * dt + g1t * sqrt(dEE_NR * dt)
    #dj = dj_tot * dt + g2t * sqrt(djj_NR * dt)

    #E_step = E + dE
    #j_step = j + dj

    dE_new = dE_NR * dt + g1t * sqrt(max(dEE_NR, 0.0) * dt)
    dj_new = dj_NR * dt + g2t * sqrt(max(djj_NR, 0.0) * dt)

    E_new = E + dE_new
    j_new = j + dj_new

    # Diffused out of the bound region -> escaped. Return now (captured = false): for
    # E >= 0, find_Lc gives Lc -> 0 and j = L/Lc -> +-Inf, which must NOT be read as a
    # loss-cone capture. The driver logs this as escaped/unbound.
    if E_new >= 0.0
        return E_new, j_new, dt, period, false, dEgw_step
    end

    rc_new, Lc_new = find_Lc(E_new)
    L_new = j_new * Lc_new

    # === TEMP (pinhole merger-count test) ===================================
    # The pinhole-aware loss-cone treatment below -- probabilistic once-per-orbit
    # pericenter gating (check_lc), the plunge per-pericenter GW burst, and the
    # pre-burst capture test -- is commented out. In its place (further down): GW
    # integrated CONTINUOUSLY and the loss cone tested DIRECTLY every step, with no
    # dt/period probability. This overcounts relative to the per-orbit treatment and
    # lets walkers "fly away"; both are intentional for this test. To restore the
    # proper treatment, un-comment the #= ... =# block and delete the naive block.
    #=
    # Reached pericenter this step? Once per orbit on average.
    check_lc = rand() < dt / period

    if plunge
        if check_lc
            # Pre-burst loss-cone test: capture is decided on the APPROACHING orbit
            # (post-diffusion, pre-GW), so it's recorded at entry rather than at the
            # deep post-burst E.
            if reflect_j(j_new) <= loss_cone(E_new)
                return E_new, reflect_j(j_new), dt, period, true
            end
            # Not captured: radiate one full orbit's GW as it passes pericenter
            # (E[burst] = (dt/P)*(dE_GW*P) = dE_GW*dt, the right mean).
            E_step = E_new + dE_GW * period
            L_step = L_new + dL_GW * period
        else
            # Coasting far from pericenter: no GW, no loss-cone test.
            E_step = E_new
            L_step = L_new
        end
    else
        # GW is a slow perturbation: integrate continuously.
        E_step = E_new + dE_GW * dt
        L_step = L_new + dL_GW * dt
    end

    rc_step, Lc_step = find_Lc(E_step)
    j_step = L_step / Lc_step

    # Normal-regime loss cone: once per orbit on the current state (continuous GW is a
    # small perturbation here, so pre/post are equivalent).
    if !plunge && check_lc && reflect_j(j_step) <= loss_cone(E_step)
        return E_step, j_step, dt, period, true
    end

    return E_step, j_step, dt, period, false
    =#

    # --- Naive treatment: continuous GW + direct loss-cone check every step ---
    E_step = E_new + dEgw_step        # dEgw_step == dE_GW * dt (the deterministic GW piece)
    L_step = L_new + dL_GW * dt

    rc_step, Lc_step = find_Lc(E_step)
    j_step = L_step / Lc_step

    # Disregard pathological j BEFORE the loss-cone test: j outside the single-reflection
    # range [-1, 2] (or non-finite) is unphysical -- a huge one-orbit diffusion kick --
    # and reflect_j would ALIAS it into the loss cone (e.g. j = 9.9 -> 2 - 9.9 = -7.9,
    # which is <= j_lc, a spurious "capture"). Return it NOT captured with the raw j so the
    # driver's pathological-j guard logs it distinctly ("pathological j = ...") and ends
    # the walker as reason -100 -- rather than aliasing it into a merger, OR (as the plain
    # sentinel did) mislabeling this fly-away as an "uncharacterizable orbit".
    if !isfinite(j_step) || j_step < -1.0 || j_step > 2.0
        return E_step, j_step, dt, period, false, dEgw_step
    end

    # Direct loss-cone capture: no dt/period probability, tested every step. Return the
    # REFLECTED j (j_step can be a small negative, e.g. -1.8e-4; reflect_j maps it to its
    # physical positive value) so the driver's log, the stored j, and its find_rp_ra
    # (L = j*Lc) all see j > 0 instead of a spurious negative.
    if reflect_j(j_step) <= loss_cone(E_step)
        return E_step, reflect_j(j_step), dt, period, true, dEgw_step
    end

    return E_step, j_step, dt, period, false, dEgw_step

end

## Set the MC step-size safety factor F_safe (the typed global mc_step reads). Use this
## rather than `global F_safe = ...` from a scope that also has an F_safe local/kwarg --
## there the `global` declaration makes every F_safe reference the global and the intended
## value is silently dropped. The argument here is named `x`, so there is no such clash.
set_F_safe!(x::Real) = (global F_safe = float(x); return nothing)

## Toggle GW-only (no-diffusion) mode the mc_step hot loop reads (see GW_ONLY in
## constants.jl). Same rationale as set_F_safe! for using a helper: a bare
## `global GW_ONLY = ...` from a scope with a local of that name silently no-ops.
set_gw_only!(x::Bool) = (global GW_ONLY = x; return nothing)

## Reflecting boundary: single reflection at j = 0 and j = 1 (for j in [-1, 2]) ##
function reflect_j(j)
    if j < 0.0
        return -j       # reflect at j = 0
    elseif j > 1.0
        return 2.0 - j  # reflect at j = 1
    else
        return j
    end
end


function run_mc(E0, j0, coef, m_obj; max_step=100, max_t=Inf)

    # max_t is in Myr; the internal clock t is in code (Henon) units, so convert once.
    max_t_code = max_t / t_conv

    # Initialize storage arrays
    E_stor = Vector{Float64}(undef, max_step + 1)
    j_stor = Vector{Float64}(undef, max_step + 1)
    t_stor = Vector{Float64}(undef, max_step + 1)

    E_stor[1], j_stor[1], t_stor[1] = E0, j0, 0.0

    # Initial values
    E, j, t = E0, j0, 0.0
    n = 1

    for step in 1:max_step
        j_lc = loss_cone(E)

        if j <= j_lc
            println("Entered loss cone (j = $j, j_lc = $j_lc, M = $m_obj)")
            flush(stdout)
            break
        elseif E >= 0.0 || t >= max_t_code
            break
        end

        E, j, dt = mc_step(E,j,coef,m_obj)

        # Pathological step (Inf/NaN or j far outside [0, 1]) -> end this
        # walker; the batch continues with the next black hole.
        if !isfinite(j) || j < -1.0 || j > 2.0
            println(stderr, "WARNING: pathological j = $j (E = $E) at step $step; ending walker early.")
            flush(stderr)
            break
        end

        j = reflect_j(j) # Reflecting boundary at j = 0 and j = 1
        t += dt
        n += 1

        E_stor[n], j_stor[n], t_stor[n] = E, j, t
        println("E = $E, j = $j, t = $(t_conv * t)")
    end

    return t_conv .* t_stor[1:n], E_stor[1:n], j_stor[1:n]
end


# compute_rp_ra: if true (default, for manual single-snapshot testing) compute rp/ra every
# step for every walker. If false (automate_mc batch runs) skip the per-step find_Lc +
# find_rp_ra -- the dominant cost -- and only backfill rp/ra for walkers that end in the
# loss cone (~1% of walkers), which are the only trajectories saved. The final rp/ra of a
# captured walker is set at the crossing regardless, so its summary rpf/raf stays correct;
# non-captured walkers get NaN rp/ra (discarded downstream).
function run_mc_rp_ra(E0, j0, coef, m_obj; max_step=100, max_t=Inf, compute_rp_ra=true)

    # max_t is in Myr; the internal clock t is in code (Henon) units, so convert once.
    max_t_code = max_t / t_conv

    # Initialize storage arrays
    E_stor = Vector{Float64}(undef, max_step + 1)
    j_stor = Vector{Float64}(undef, max_step + 1)
    t_stor = Vector{Float64}(undef, max_step + 1)
    
    rp_stor = Vector{Float64}(undef, max_step + 1)
    ra_stor = Vector{Float64}(undef, max_step + 1)

    # Cumulative GW energy dissipated (specific, per unit m_obj), one entry per stored step.
    # Egw_stor[k] = sum over steps 1..k of dE_GW*dt (deterministic GW E-contribution); starts
    # at 0 at the IC. Negative (GW removes orbital energy); |value| grows with radiated energy.
    Egw_stor = Vector{Float64}(undef, max_step + 1)

    rc0, Lc0 = find_Lc(E0)
    rp_stor[1], ra_stor[1] = find_rp_ra(E0, j0*Lc0)


    E_stor[1], j_stor[1], t_stor[1] = E0, j0, 0.0
    Egw_stor[1] = 0.0

    # Initial values
    E, j, t = E0, j0, 0.0
    Egw_cum = 0.0
    n = 1
    reason = 3   # default: reached max_step

    # Initial-condition loss-cone check (a walker could start inside).
    j_lc0 = loss_cone(E)
    if j <= j_lc0
        reason = 1
        println("Entered loss cone (j = $j, j_lc = $j_lc0, M = $m_obj, Egw = $Egw_cum)")
        flush(stdout)
        return t_conv .* t_stor[1:n], E_stor[1:n], j_stor[1:n],
               rp_stor[1:n], ra_stor[1:n], reason, Egw_stor[1:n]
    end

    for step in 1:max_step

        if E <= -1000000000
            reason = -100
            println(stderr, "REASON -100: below potential minimum (E = $E) at step $step; ending walker early.")
            flush(stderr)
            break
        elseif E >= 0.0
            reason = -100
            println(stderr, "REASON -100: escaped/unbound at loop top (E = $E) at step $step; ending walker early.")
            flush(stderr)
            break
        elseif t >= max_t_code
            reason = 2
            break
        end

        E_new, j_new, dt, P_orb, captured, dEgw_step = mc_step(E,j,coef,m_obj)

        # mc_step could not characterize the orbit (find_rp_ra / find_period failed);
        # it returns the sentinel (-100, -100, 1000). Report the real (E, j).
        if E_new == -100.0 && j_new == -100.0
            reason = -100
            println(stderr, "REASON -100: orbit uncharacterizable at step $step (E = $E, j = $j); ending walker early.")
            flush(stderr)
            break
        end

        # Captured: mc_step tested the loss cone at pericenter (pre-burst). Store the
        # crossing point (E_new, j_new) as the FINAL step before stopping -- otherwise
        # es/js/rps/ras end one step short of the loss cone and plots never reach it.
        if captured
            reason = 1
            n += 1
            E_stor[n], j_stor[n], t_stor[n] = E_new, j_new, t + dt
            Egw_cum += dEgw_step
            Egw_stor[n] = Egw_cum
            _, Lc_c = find_Lc(E_new)
            rp_stor[n], ra_stor[n] = find_rp_ra(E_new, j_new * Lc_c)
            println("Entered loss cone (j = $j_new, E = $E_new, M = $m_obj, Egw = $Egw_cum) at step $step")
            flush(stdout)
            break
        end

        # Pathological step (Inf/NaN or j far outside [0, 1]) -> end this walker.
        if !isfinite(j_new) || j_new < -1.0 || j_new > 2.0
            reason = -100
            println(stderr, "REASON -100: pathological j = $j_new (E = $E_new) at step $step; ending walker early.")
            flush(stderr)
            break
        end

        # Completing this step would pass max_t: discard it and time out.
        if t + dt > max_t_code
            reason = 2
            break
        end

        # Diffused out of the bound region (E >= 0) -> object escaped.
        if E_new >= 0.0
            reason = -100
            println(stderr, "REASON -100: escaped/unbound after step (E = $E_new) at step $step; ending walker early.")
            flush(stderr)
            break
        end

        E = E_new
        # Adaptive j floor: below ~2e-8 the eccentricity rounds to e = 1.0 and the orbit
        # becomes uncharacterizable, and a j well inside the loss cone is a plunge where
        # the period is unchanged with j at fixed E. Clamp at 10% of j_lc(E) (deep inside
        # the cone, but adapts to E instead of a fixed 1e-6); never below 1e-7 so the
        # orbit stays characterizable. The loss-cone capture is decided in mc_step.
        j = max(reflect_j(j_new), 0.1 * loss_cone(E_new), 1.0e-7)
        t += dt
        Egw_cum += dEgw_step
        n += 1

        E_stor[n], j_stor[n], t_stor[n] = E, j, t
        Egw_stor[n] = Egw_cum
        # Per-step rp/ra (find_Lc + find_rp_ra) is the dominant cost but only matters for
        # walkers that reach the loss cone. With compute_rp_ra=false skip it here (store NaN)
        # and backfill the captured trajectory after the loop.
        if compute_rp_ra
            rc, Lc = find_Lc(E)
            rp, ra = find_rp_ra(E, j*Lc)
            rp_stor[n], ra_stor[n] = rp, ra
        else
            rp_stor[n], ra_stor[n] = NaN, NaN
        end
        #println("E = $E, j = $j, t = $t")
    end

    # Reached the step cap without capture (1), timeout (2), or an error (-100): reason stays
    # at its default 3. Log it so these show up in the error log like the other exit reasons.
    if reason == 3
        println(stderr, "REASON 3: reached max_step ($max_step) without capture or timeout " *
                        "(E = $E, j = $j, M = $m_obj); ending walker.")
        flush(stderr)
    end

    # If per-step rp/ra were skipped but this walker ended in the loss cone, backfill its
    # full trajectory now (captured walks are ~1% and short, so this is cheap). This
    # reproduces exactly what the per-step path would have stored. The IC-start-in-cone
    # path returns earlier with rp_stor[1] already set at the top, so it is unaffected.
    if !compute_rp_ra && reason == 1
        @inbounds for i in 1:n
            _, Lc_i = find_Lc(E_stor[i])
            rp_stor[i], ra_stor[i] = find_rp_ra(E_stor[i], j_stor[i] * Lc_i)
        end
    end

    return t_conv .* t_stor[1:n], E_stor[1:n], j_stor[1:n], rp_stor[1:n], ra_stor[1:n], reason, Egw_stor[1:n]
end


## --- Mean-drift "flow field" diagnostic (Fokker-Planck streamlines) --- ##
#
# Trace the ENSEMBLE-MEAN trajectory of a walker in (E, j). At each point, run `n_furball`
# independent single mc_step realizations from the SAME (E, j), collapse that cloud (the
# "furball") to its centroid, and step deterministically to the centroid; then repeat from
# there. Averaging cancels the diffusive spread and leaves only the systematic drift, so
# the path is a streamline of the mean flow -- a clean way to see whether the (loaded)
# diffusion coefficients push walkers sensibly toward the loss cone or do something
# pathological (stall, reverse, run away in E).
#
# All mc_step physics is kept (drift + GW + reflection + loss-cone capture). The centroid
# includes captured and escaped members, so the mean flows into the loss cone naturally as
# the capture fraction turns on near the boundary. Because every furball member starts from
# the same (E, j), they share one deterministic timestep dt, so the recorded time is the
# drift time. There is no physical time limit -- termination is loss cone / escape / max_step.
#
# Returns per-step arrays E, j, t[Myr], rp, ra, captured-fraction, and a termination reason
# (1 = loss cone, 2 = unused, -100 = escaped / uncharacterizable, 3 = reached max_step).
function mc_mean_flow(E0, j0, coef, m_obj;
                      n_furball::Integer = 1000, max_step::Integer = 5000)

    E_stor    = Float64[]
    j_stor    = Float64[]
    t_stor    = Float64[]      # cumulative drift time [Myr]
    rp_stor   = Float64[]
    ra_stor   = Float64[]
    fcap_stor = Float64[]

    E, j, t = E0, j0, 0.0
    reason  = 3                # default: reached max_step

    # Record the starting point.
    _, Lc0 = find_Lc(E)
    rp0, ra0 = find_rp_ra(E, clamp(j, 1.0e-7, 1.0) * Lc0)
    push!(E_stor, E); push!(j_stor, j); push!(t_stor, t_conv * t)
    push!(rp_stor, rp0); push!(ra_stor, ra0); push!(fcap_stor, 0.0)

    for step in 1:max_step

        # Terminate on the mean state.
        j_lc = loss_cone(E)
        if isfinite(j_lc) && j <= j_lc
            reason = 1
            break
        elseif E >= 0.0 || E <= -1.0e7
            reason = -100
            break
        end

        # Furball: n_furball single mc_steps from (E, j), collapsed to the centroid.
        sumE = 0.0; sumj = 0.0; sumdt = 0.0
        nvalid = 0; ncap = 0
        for _ in 1:n_furball
            E_new, j_new, dt, period, captured = mc_step(E, j, coef, m_obj)

            (E_new == -100.0 && j_new == -100.0) && continue          # sentinel orbit
            (isfinite(E_new) && isfinite(j_new) && isfinite(dt)) || continue

            jr = reflect_j(j_new)                                     # captured j already reflected
            (jr < -1.0 || jr > 2.0) && continue                       # pathological

            sumE += E_new; sumj += jr; sumdt += dt
            nvalid += 1
            captured && (ncap += 1)
        end

        if nvalid == 0
            reason = -100
            break
        end

        E = sumE / nvalid
        j = sumj / nvalid
        t += sumdt / nvalid        # dt is identical across members -> this is that dt

        _, Lc = find_Lc(E)
        rp, ra = find_rp_ra(E, clamp(j, 1.0e-7, 1.0) * Lc)
        push!(E_stor, E); push!(j_stor, j); push!(t_stor, t_conv * t)
        push!(rp_stor, rp); push!(ra_stor, ra); push!(fcap_stor, ncap / n_furball)
    end

    return E_stor, j_stor, t_stor, rp_stor, ra_stor, fcap_stor, reason
end


## --- One furball "step": mean displacement (drift vector) at a fixed (E, j) --- ##
#
# Run `n_furball` independent single mc_steps from the SAME (E, j) and return the MEAN
# displacement (dE, dj) -- the local drift vector of the Fokker-Planck flow -- plus the
# captured fraction. Averaging cancels the diffusive spread, leaving the systematic
# direction. Used to build a gradient/quiver field over an (E, j) grid: each grid point
# gets exactly one furball, with NO time integration (see run_mean_flow.jl). Captured
# members (which jump to the loss cone) are included, so the vector turns into the cone
# near the boundary. Returns NaN displacement if no member produced a usable step.
function mc_furball_step(E, j, coef, m_obj; n_furball::Integer = 1000)

    # Rule 1: a grid point that STARTS inside (or on) the loss cone is already a plunge --
    # don't compute a vector there at all. (loss_cone is NaN where undefined, e.g. E >= 0;
    # treat that as skip too.)
    j_lc = loss_cone(E)
    if !isfinite(j_lc) || j <= j_lc
        return (dE = NaN, dj = NaN, fcap = NaN, n = 0)
    end

    # Rule 2: a grid point OUTSIDE the loss cone runs the full furball and counts EVERY
    # random walk -- including members that diffuse INTO the loss cone (captured=true, whose
    # returned j is the loss-cone value) -- so the mean vector reflects the true flux toward
    # the boundary. Only pathological fly-aways (j outside [-1, 2]) are dropped.
    sumdE = 0.0; sumdj = 0.0
    nvalid = 0; ncap = 0
    for _ in 1:n_furball
        E_new, j_new, dt, period, captured = mc_step(E, j, coef, m_obj)
        (E_new == -100.0 && j_new == -100.0) && continue          # sentinel orbit
        (isfinite(E_new) && isfinite(j_new)) || continue
        jr = reflect_j(j_new)                                     # captured j already reflected
        (jr < -1.0 || jr > 2.0) && continue                       # pathological fly-away
        sumdE += E_new - E
        sumdj += jr - j
        nvalid += 1
        captured && (ncap += 1)
    end
    nvalid == 0 && return (dE = NaN, dj = NaN, fcap = NaN, n = 0)
    return (dE = sumdE / nvalid, dj = sumdj / nvalid, fcap = ncap / n_furball, n = nvalid)
end


## Semimajor axis a = (rp+ra)/2 at (E, j); NaN if the orbit can't be characterized. ##
function _semimajor(E, j)
    (isfinite(E) && E < 0.0 && isfinite(j) && j > 0.0) || return NaN
    _, Lc = find_Lc(E)
    (isfinite(Lc) && Lc > 0.0) || return NaN
    rp, ra = find_rp_ra(E, j * Lc)
    orbit_failed(rp, ra) && return NaN
    return (rp + ra) / 2
end


## --- Deterministic drift vector at (E, j): the Fokker-Planck first moment (a RATE) --- ##
#
# The instantaneous mean drift (dE/dt, dj/dt) = diffusion drift (Dj1/DE1 coeffs, already in
# (E,j) space) + GW drift (Peters, transformed from (E,L) to j-space via the Jacobian of
# j = L/Jc(E): dj/dt = (dL/dt - j (dJc/dE) dE/dt)/Jc). This is the EXACT drift field -- no
# stochastic step, no timestep, no capture/overshoot exclusion -- so it is the right
# quantity for a gradient/quiver plot, and GW cleanly dominates near the loss cone where
# the Peters rate blows up. Rule 1: NaN inside the loss cone. Returns the total drift plus
# the GW-only and diffusion-only pieces so their relative size can be checked directly.
function mc_drift(E, j, coef, m_obj)
    nanret = (dE = NaN, dj = NaN, dE_gw = NaN, dj_gw = NaN, dE_nr = NaN, dj_nr = NaN,
              a = NaN, da = NaN, da_gw = NaN, da_nr = NaN)

    j_lc = loss_cone(E)
    (!isfinite(j_lc) || j <= j_lc) && return nanret          # Rule 1: inside the cone

    rc, Lc = find_Lc(E)
    L = j * Lc
    rp, ra = find_rp_ra(E, L)
    orbit_failed(rp, ra) && return nanret
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    (!isfinite(a) || !isfinite(e) || a <= 0 || e < 0 || e >= 1) && return nanret

    # Diffusion drift (orbit-averaged, already transformed into (E,j) space).
    dE_nr = coef.DE1_samp(j, E) + m_obj * coef.DE2_samp(j, E)
    dj_nr = coef.Dj1_samp(j, E) + m_obj * coef.Dj2_samp(j, E)

    # GW (Peters) drift, transformed (E, L) -> j-space.
    dE_gw, dL_gw = gw_rates(a, e, m_obj)
    dj_gw = (dL_gw - j * dLc_dE(E) * dE_gw) / Lc

    # (E, j) -> (a, j) transform: a = (rp+ra)/2, so da/dt = (da/dE) dE/dt + (da/dj) dj/dt.
    # Central finite differences with tiny relative steps -- the perturbed points stay
    # outside the cone (j > j_lc), so they characterize; if one fails, that partial (hence
    # da) becomes NaN, which the plot skips. The grid stays (E, j); a/da are just for
    # plotting in (a, j).
    hE = 1.0e-4 * abs(E)
    hj = 1.0e-4 * j
    dAdE = (_semimajor(E + hE, j) - _semimajor(E - hE, j)) / (2 * hE)
    dAdj = (_semimajor(E, j + hj) - _semimajor(E, j - hj)) / (2 * hj)

    da    = dAdE * (dE_nr + dE_gw) + dAdj * (dj_nr + dj_gw)
    da_gw = dAdE * dE_gw + dAdj * dj_gw
    da_nr = dAdE * dE_nr + dAdj * dj_nr

    return (dE = dE_nr + dE_gw, dj = dj_nr + dj_gw,
            dE_gw = dE_gw, dj_gw = dj_gw, dE_nr = dE_nr, dj_nr = dj_nr,
            a = a, da = da, da_gw = da_gw, da_nr = da_nr)
end


## Loss cone prescription from Qunbar / Stone ##
function loss_cone(E)
    J_icbo = (4 * G * M_bh) / c
    Jc = last(find_Lc(E))

    # Degenerate circular angular momentum -> loss cone undefined here. Returning
    # NaN makes j <= j_lc false, so the walker fails loudly (-100) downstream
    # instead of being recorded as a spurious capture via j_lc = Inf.
    if !isfinite(Jc) || Jc <= 0.0
        println(stderr, "WARNING: find_Lc gave Jc = $Jc at E = $E; loss cone undefined here.")
        return NaN
    end

    return J_icbo / Jc
end


## Loss-cone boundary j_lc(E) on a log-spaced |E| grid (for plotting the curve) ##
function loss_cone_curve(; Emin_abs = 0.01, Emax_abs = 1.0e7, nbins = 300)
    E_abs = 10.0 .^ range(log10(Emin_abs), log10(Emax_abs), length = nbins)
    lcE = -E_abs                          # bound energies (negative): -0.1 .. -1e6
    lcJ = [loss_cone(E) for E in lcE]     # normalized loss-cone j via find_Lc
    return lcE, lcJ
end


## (E, L, j) for a phase-space point using the current global potential ##
function ej_from_phase(r, vr, vt)
    E = 0.5 * (vr^2 + vt^2) + psi_calc(r)
    L = r * vt

    if E >= 0.0
        return E, L, NaN            # unbound orbit: j undefined
    end

    Lc = last(find_Lc(E))
    return E, L, L / Lc
end


## Initial (E, j, mass) for every bound compact objectin a snapshot ##
function black_hole_ics(dat; startypes = (14))

    inds = Int[]
    E0s  = Float64[]
    J0s  = Float64[]
    m0s  = Float64[]

    for i in eachindex(dat.startype)
        dat.startype[i] in startypes || continue    # target types only

        E, L, j = ej_from_phase(dat.r[i], dat.vr[i], dat.vt[i])
        isfinite(j) || continue                      # skip unbound

        push!(inds, i)
        push!(E0s, E)
        push!(J0s, reflect_j(j))                     # keep j0 within [0, 1]
        push!(m0s, dat.m[i])                         # each BH uses its own mass
    end

    return (indices = inds, E0 = E0s, J0 = J0s, m0 = m0s)
end


# ============================================================================
# METHOD 4: physically-motivated deep-E extrapolation of the diffusion coefficients.
#
# Peer to methods 1 (:zero), 2 (extrapolate_coeffs, power-law), 3 (:clamp). Instead of
# guessing the coefficients below the innermost sampled particle, we exploit that this region
# is deep inside the influence radius -> the potential is ~Keplerian and a relaxed cusp is a
# power law. We fit the cusp from the snapshot and replace the deep tail with ANALYTIC
# coefficients built from that fit (the same machinery validated by the analytic-sampling test).
#
# Two power-law slopes are fit (both feed the diffusion coefficients, see coefficients.jl):
#   gamma_rho  from the MASS density   rho(r)   = A_rho r^-gamma_rho   -> DF.N (friction: DE2/Dj2)
#   gamma_F    from the MASS^2 density  rho2(r) = A_F   r^-gamma_F     -> DF.F (diffusion: DEE/Djj/DEj, DE1/Dj1)
# In a single-mass cusp these coincide (F = m*N); under mass segregation gamma_F (BH-dominated)
# is steeper and is the one that actually drives loss-cone diffusion.
# ============================================================================

## Beta B(a,b) via midpoint quadrature (no SpecialFunctions dep; distinct name from _beta). ##
function _beta_pl(a::Real, b::Real; n::Integer = 40000)
    s = 0.0; h = 1.0 / n
    @inbounds for i in 1:n
        x = (i - 0.5) * h
        s += x^(a - 1) * (1 - x)^(b - 1)
    end
    return s * h
end

## Fit a power-law density rho_w(r) = A r^(-gamma) to weighted particles over [r_lo, r_hi].
## `weights` = masses for the mass density (first moment), masses.^2 for the mass^2 density
## (second moment). Bins particles in log r, forms the shell density per bin, and does a
## count-weighted log-log least-squares fit (same estimator as the analytic-sampling density
## plot). Returns (gamma, A) in the SAME (code) units as r and the weights. ##
function fit_powerlaw_density(rs::AbstractVector, weights::AbstractVector,
                              r_lo::Real, r_hi::Real; nbins::Integer = 18)
    r_lo > 0 && r_hi > r_lo || error("fit_powerlaw_density: need 0 < r_lo < r_hi (got $r_lo, $r_hi)")
    edges = 10.0 .^ range(log10(r_lo), log10(r_hi), length = nbins + 1)
    vol   = (4pi / 3) .* (edges[2:end] .^ 3 .- edges[1:end-1] .^ 3)   # shell volumes
    wsum  = zeros(nbins); cnt = zeros(Int, nbins)
    @inbounds for i in eachindex(rs)
        r = rs[i]
        (r < r_lo || r >= r_hi) && continue
        b = searchsortedlast(edges, r)
        (1 <= b <= nbins) || continue
        wsum[b] += weights[i]; cnt[b] += 1
    end
    good = (cnt .> 0) .& (wsum .> 0)
    count(good) >= 3 || error("fit_powerlaw_density: <3 populated bins in [$r_lo, $r_hi]; " *
                              "widen the range or lower nbins (got $(count(good)))")
    rc   = sqrt.(edges[1:end-1] .* edges[2:end])[good]        # geometric bin centers
    rho  = (wsum ./ vol)[good]
    w    = Float64.(cnt[good])                                 # count-weighted fit
    lx = log10.(rc); ly = log10.(rho)
    xb = sum(w .* lx) / sum(w); yb = sum(w .* ly) / sum(w)
    slope = sum(w .* (lx .- xb) .* (ly .- yb)) / sum(w .* (lx .- xb) .^ 2)
    gamma = -slope
    A = 10.0 ^ (yb - slope * xb)
    return gamma, A
end

## Two-moment analytic DF (as a SampDF) for a Keplerian power-law cusp: tabn from (A_rho,
## gamma_rho), tabf INDEPENDENTLY from (A_F, gamma_F). Each is the closed-form Eddington DF of
## a power law in a Kepler potential: N_log = K r^3 v^3 eps^(gamma-3/2), eps = M_bh/r - v^2/2
## (code, G=1), K = (2 sqrt(2) pi / B(gamma-1/2, 3/2)) * A * M_bh^(-gamma). Grid spans [r_lo,
## r_hi] with v up to ~v_esc(r_lo). ##
function build_powerlaw_df_two_moment(gamma_rho::Real, A_rho::Real, gamma_F::Real, A_F::Real,
                                      M_bh_code::Real, r_lo::Real, r_hi::Real;
                                      nr::Integer = 600, nv::Integer = 600)
    p_rho = gamma_rho - 1.5; p_F = gamma_F - 1.5
    K_rho = (2 * sqrt(2) * pi / _beta_pl(p_rho + 1, 1.5)) * A_rho * M_bh_code ^ (-gamma_rho)
    K_F   = (2 * sqrt(2) * pi / _beta_pl(p_F   + 1, 1.5)) * A_F   * M_bh_code ^ (-gamma_F)

    v_hi = sqrt(2 * M_bh_code / r_lo); v_lo = v_hi / 1.0e5
    logr = collect(range(log(r_lo), log(r_hi), length = nr))
    logv = collect(range(log(v_lo), log(v_hi), length = nv))
    tabn = zeros(nr, nv); tabf = zeros(nr, nv)
    @inbounds for i in 1:nr
        r = exp(logr[i])
        for k in 1:nv
            v = exp(logv[k])
            eps = M_bh_code / r - 0.5 * v^2
            if eps > 0.0
                rv3 = r^3 * v^3
                tabn[i, k] = K_rho * rv3 * eps^p_rho
                tabf[i, k] = K_F   * rv3 * eps^p_F
            end
        end
    end
    Nfun = (lr, lv) -> bilinear_interp(lr, lv, logr, logv, tabn)
    Ffun = (lr, lv) -> bilinear_interp(lr, lv, logr, logv, tabf)
    return SampDF(logr, logv, tabn, tabf, Nfun, Ffun)
end

## Method 4 extrapolation: keep the sampled coefficient rows at E >= E_fit_deep (identical to
## method 2's kept region) and replace the deep tail (down to E_end) with analytic coefficients.
## Builds the SAME n_E_total E-axis as extrapolate_coeffs, so a method-2 and a method-4 grid
## differ ONLY in the tail. Requires the SAMPLED orbit potential to be active on entry (to
## compute nothing here -- E_fit_deep is passed in, computed by the caller in that potential).
## The analytic tail is evaluated in a pure-Kepler potential (rebuilt internally down to r_g/10)
## with the two-moment fitted DF; the caller's potential + VEL_INT_TABLES are restored on exit.
## `dat_s` = the snapshot (for the gamma fits + r_infl); `c_code` sets r_g and the default E_end.
## Returns the method-4 DiffusionCoeffs. Fitted slopes are printed (verbose) for inspection.
function extrapolate_coeffs_method4(coef::DiffusionCoeffs, dat_s, M_bh_code::Real, c_code::Real,
                                    E_fit_deep::Real;
                                    E_end::Real = -c_code^2 / 2, n_E_total::Integer = 300,
                                    r_fit_hi_frac::Real = 0.1, r_fit_lo::Union{Nothing,Real} = nothing,
                                    nr::Integer = 600, nv::Integer = 600, fit_nbins::Integer = 18,
                                    verbose::Bool = true)

    E_tab = collect(float.(coef.E_tab)); j_tab = collect(float.(coef.j_tab))
    n_j = length(j_tab)
    E_end < E_fit_deep || error("extrapolate_coeffs_method4: E_end ($E_end) must be below E_fit_deep ($E_fit_deep)")

    # Kept (reliable) rows + deep tail energies -- IDENTICAL construction to extrapolate_coeffs.
    keep_mask = E_tab .>= E_fit_deep
    E_keep = E_tab[keep_mask]; n_keep = length(E_keep)
    n_keep > 0 || error("extrapolate_coeffs_method4: no grid rows at/above E_fit_deep = $E_fit_deep")
    n_tail = n_E_total - n_keep
    n_tail >= 1 || error("extrapolate_coeffs_method4: n_E_total ($n_E_total) <= kept rows ($n_keep)")
    E_tail = -exp10.(range(log10(abs(E_end)), log10(abs(E_keep[1])), length = n_tail + 1))[1:end-1]
    E_new  = vcat(E_tail, E_keep)

    # Influence radius (where enclosed STELLAR mass = M_bh) + two-slope cusp fit over [r_lo, 0.1 r_infl].
    _, _, M_tab_s = find_psi_arrays(dat_s.r, dat_s.m)
    idx    = clamp(searchsortedfirst(view(M_tab_s, 1:length(dat_s.r)), M_bh_code), 1, length(dat_s.r))
    r_infl = dat_s.r[idx]
    r_hi   = r_fit_hi_frac * r_infl
    r_lo   = r_fit_lo === nothing ? minimum(dat_s.r) : r_fit_lo
    gamma_rho, A_rho = fit_powerlaw_density(dat_s.r, dat_s.m,       r_lo, r_hi; nbins = fit_nbins)
    gamma_F,  A_F    = fit_powerlaw_density(dat_s.r, dat_s.m .^ 2,  r_lo, r_hi; nbins = fit_nbins)

    r_g = M_bh_code / c_code^2
    if verbose
        println("  [method4] r_infl = $r_infl (code),  fit range [$r_lo, $r_hi] = [min r, $(r_fit_hi_frac) r_infl]")
        println("  [method4] mass density    : gamma_rho = $(round(gamma_rho, digits=4)),  A_rho = $A_rho")
        println("  [method4] mass^2 density   : gamma_F   = $(round(gamma_F,  digits=4)),  A_F   = $A_F")
        println("  [method4] r_g = $r_g,  E_end = $E_end  (= -c^2/2, circular energy at r_g; c = $c_code)")
        println("  [method4] keep $n_keep sampled rows (E >= $E_fit_deep); analytic tail = $n_tail rows down to E_end")
        flush(stdout)
    end

    # --- analytic tail in a pure-Kepler potential; restore the caller's state afterwards ---
    saved_psi = ORBIT_PSI[]; saved_vit = VEL_INT_TABLES[]
    # r-range spanning the tail orbits: from below r_g out past the shallowest tail apocenter
    # (~2*a for the least-bound tail energy; use a generous 10*a = 5*M_bh/|E_keep[1]|).
    r_span_lo = r_g / 10
    r_span_hi = 5 * M_bh_code / abs(E_keep[1])
    tail_mats = try
        rgrid = 10.0 .^ range(log10(r_span_lo), log10(r_span_hi), length = 6000)
        kpsi, ktot, kM = find_psi_arrays(rgrid, zeros(length(rgrid)))
        ORBIT_PSI[] = (r = rgrid, psi_tab = kpsi, psi_tot = ktot, M_tab = kM,
                       psi_calc = psi_exact(rgrid, kpsi, kM))
        activate_orbit_psi!()

        df = build_powerlaw_df_two_moment(gamma_rho, A_rho, gamma_F, A_F, M_bh_code,
                                          r_span_lo, r_span_hi; nr = nr, nv = nv)
        VEL_INT_TABLES[] = build_vel_int_tables(df)   # cou_log from global M_bh/dat (matches sampled)

        DE1t = zeros(n_j, n_tail); DE2t = zeros(n_j, n_tail); Dj1t = zeros(n_j, n_tail)
        Dj2t = zeros(n_j, n_tail); DEEt = zeros(n_j, n_tail); Djjt = zeros(n_j, n_tail)
        DEjt = zeros(n_j, n_tail)
        Threads.@threads for i in 1:n_tail
            E = E_tail[i]
            for k in 1:n_j
                d1, d2, jj1, jj2, ee, jj, ej = avg_coef_ej(E, j_tab[k])
                @inbounds begin
                    DE1t[k, i] = d1; DE2t[k, i] = d2; Dj1t[k, i] = jj1; Dj2t[k, i] = jj2
                    DEEt[k, i] = ee; Djjt[k, i] = jj; DEjt[k, i] = ej
                end
            end
        end
        (DE1t, DE2t, Dj1t, Dj2t, DEEt, Djjt, DEjt)
    finally
        if saved_psi !== nothing
            ORBIT_PSI[] = saved_psi; activate_orbit_psi!()
        end
        VEL_INT_TABLES[] = saved_vit
    end

    # Assemble: analytic tail (cols 1:n_tail) + kept sampled rows (cols n_tail+1:end).
    sampled = (coef.DE1_tab, coef.DE2_tab, coef.Dj1_tab, coef.Dj2_tab,
               coef.DEE_tab, coef.Djj_tab, coef.DEj_tab)
    out = Matrix{Float64}[]
    for (Msamp, Mtail) in zip(sampled, tail_mats)
        Mn = Matrix{Float64}(undef, n_j, n_E_total)
        Mn[:, 1:n_tail]         = Mtail
        Mn[:, (n_tail + 1):end] = Msamp[:, keep_mask]
        push!(out, Mn)
    end

    DE1, DE2, Dj1, Dj2, DEE, Djj, DEj = out
    interp(Z) = (j, E) -> bilinear_interp_clamped(j, E, j_tab, E_new, Z)
    return DiffusionCoeffs(DE1, DE2, Dj1, Dj2, DEE, Djj, DEj,
        E_new, j_tab, E_new[1], coef.emax,
        interp(DE1), interp(DE2), interp(Dj1), interp(Dj2),
        interp(DEE), interp(Djj), interp(DEj))
end

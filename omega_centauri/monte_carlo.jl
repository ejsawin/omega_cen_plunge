# --- Monte Carlo Functions --- #

function mc_step(E,j,coef,m_obj)

    # Characterize orbit 
    rc, Lc = find_Lc(E)
    L = j * Lc

    rp, ra = find_rp_ra(E, L)
    
    if orbit_failed(rp, ra)
        # Return whatever your code uses to mark a terminated/rejected orbit.
        # Do not calculate a, e, period, or GW rates.
        return -100.0, -100.0, 1000, 1000, false
    end
    
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    
    # Extra protection against bad elements.
    if !isfinite(a) || !isfinite(e) || a <= 0 || e < 0 || e >= 1
        return -100.0, -100.0, 1000, 1000, false
    end
    
    period = find_period(E, L)
    
    if !isfinite(period) || period <= 0
        return -100.0, -100.0, 1000, 1000, false
    end
    
    period = find_period(E,L)
    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)

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

    # Step Size
    F_safe = 10 # step covers N = max(1, N_safe) orbits (paper default; tested 3-30)

    # Paper's two timescales: GW on E, and L-diffusion (in j here). t_jj is clamped
    # >= 0 so noisy interpolation can't give a negative timescale (-> Inf = unconstrained).
    t_gwE = abs(E / dE_GW)           # GW on E
    t_jj  = j^2 / max(djj_NR, 0.0)   # L-diffusion (in j)

    # Two regimes. When the GW timescale drops below one orbit (t_gwE < period),
    # dE_GW*period > |E|: the orbit-averaged Peters rate is past its validity, because
    # GW is really a discrete kick delivered once per orbit at pericenter, not a smooth
    # drift. There we step on the diffusion-limited dt (GW excluded) and deliver the
    # full per-orbit GW burst only when the walker actually reaches pericenter this
    # step (probability check_lc = dt/period); the driver tests the loss cone on that
    # same event. Otherwise (GW slow, many orbits) integrate GW continuously.
    #
    # Paper: N_safe = (1/F_safe) min(E/dE_GW, sqrt(L^2/(<dL^2>_o P))). The linear GW
    # term is dropped in the plunge regime, where GW is delivered as a per-pericenter
    # burst rather than limiting the step.
    plunge = t_gwE < period

    if plunge
        N_safe = (1/F_safe) * sqrt(t_jj / period)
        # Cap at one orbit: off-grid (djj = 0) gives N_safe = Inf, and dt/period must
        # stay <= 1 to be a valid pericenter-reach probability.
        dt = min(N_safe, 1.0) * period
    else
        N_safe = (1/F_safe) * min(t_gwE / period, sqrt(t_jj / period))
        dt = N_safe * period
    end

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
        return E_new, j_new, dt, period, false
    end

    rc_new, Lc_new = find_Lc(E_new)
    L_new = j_new * Lc_new

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

end

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
            println("Entered loss cone (j = $j, j_lc = $j_lc)")
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


function run_mc_rp_ra(E0, j0, coef, m_obj; max_step=100, max_t=Inf)

    # max_t is in Myr; the internal clock t is in code (Henon) units, so convert once.
    max_t_code = max_t / t_conv

    # Initialize storage arrays
    E_stor = Vector{Float64}(undef, max_step + 1)
    j_stor = Vector{Float64}(undef, max_step + 1)
    t_stor = Vector{Float64}(undef, max_step + 1)
    
    rp_stor = Vector{Float64}(undef, max_step + 1)
    ra_stor = Vector{Float64}(undef, max_step + 1)

    rc0, Lc0 = find_Lc(E0)
    rp_stor[1], ra_stor[1] = find_rp_ra(E0, j0*Lc0)

    
    E_stor[1], j_stor[1], t_stor[1] = E0, j0, 0.0

    # Initial values
    E, j, t = E0, j0, 0.0
    n = 1
    reason = 3   # default: reached max_step

    # Initial-condition loss-cone check (a walker could start inside).
    j_lc0 = loss_cone(E)
    if j <= j_lc0
        reason = 1
        println("Entered loss cone (j = $j, j_lc = $j_lc0)")
        return t_conv .* t_stor[1:n], E_stor[1:n], j_stor[1:n],
               rp_stor[1:n], ra_stor[1:n], reason
    end

    for step in 1:max_step

        if E <= -10000000
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

        E_new, j_new, dt, P_orb, captured = mc_step(E,j,coef,m_obj)

        # mc_step could not characterize the orbit (find_rp_ra / find_period failed);
        # it returns the sentinel (-100, -100, 1000). Report the real (E, j).
        if E_new == -100.0 && j_new == -100.0
            reason = -100
            println(stderr, "REASON -100: orbit uncharacterizable at step $step (E = $E, j = $j); ending walker early.")
            flush(stderr)
            break
        end

        # Captured: mc_step tested the loss cone at pericenter (pre-burst).
        if captured
            reason = 1
            println("Entered loss cone (j = $j_new, E = $E_new) at step $step")
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
        n += 1

        E_stor[n], j_stor[n], t_stor[n] = E, j, t
        rc, Lc = find_Lc(E)

        rp, ra = find_rp_ra(E, j*Lc)
        rp_stor[n], ra_stor[n] = rp, ra
        #println("E = $E, j = $j, t = $t, rp = $rp, ra = $ra")
    end

    return t_conv .* t_stor[1:n], E_stor[1:n], j_stor[1:n], rp_stor[1:n], ra_stor[1:n], reason
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

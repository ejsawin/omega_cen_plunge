# --- Monte Carlo Functions --- #

function mc_step(E,j,coef,m_obj)

    # Characterize orbit 
    rc, Lc = find_Lc(E)
    L = j * Lc

    rp, ra = find_rp_ra(E, L)
    
    if orbit_failed(rp, ra)
        # Return whatever your code uses to mark a terminated/rejected orbit.
        # Do not calculate a, e, period, or GW rates.
        return -100.0, -100.0, 1000
    end
    
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)
    
    # Extra protection against bad elements.
    if !isfinite(a) || !isfinite(e) || a <= 0 || e < 0 || e >= 1
        return -100.0, -100.0, 1000
    end
    
    period = find_period(E, L)
    
    if !isfinite(period) || period <= 0
        return -100.0, -100.0, 1000
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

    # Adaptive timestep (Qunbar & Stone)
    F_safe = 50

    gw_term = E / dE_GW 
#    nr_term = sqrt((j^2) / (djj_NR * period))
    nr_term = sqrt((j^2) / (djj_NR))


    
    N_safe = (1/F_safe) * min(gw_term, nr_term) # Choose smaller timestep
#    println("GW: $gw_term, NR: $nr_term")
#    N_orb = max(1.0, N_safe)
#    dt = N_orb * period 
    dt = max(N_safe,period)

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

    dE_new = dE_NR * dt + g1t * sqrt(dEE_NR * dt)
    dj_new = dj_NR * dt + g2t * sqrt(djj_NR * dt)
    
    E_new = E + dE_new
    j_new = j + dj_new
    
    rc_new, Lc_new = find_Lc(E_new)
    L_new = j_new * Lc_new
    
    E_step = E_new + dE_GW * dt
    L_step = L_new + dL_GW * dt
    
    rc_step, Lc_step = find_Lc(E_step)
    j_step = L_step / Lc_step
    
    return E_step, j_step, dt
    
end

function run_mc(E0, j0, coef, m_obj; max_step=100, max_t=Inf)

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

        if E <= -10000000
            println("Below potential minimum: (E = $E)")
            break
        elseif j <= j_lc
            println("Entered loss cone (j = $j, j_lc = $j_lc)")
            break
        elseif j >= 1.0 || E >= 0.0 || t >= max_t
            break
        end

        E, j, dt = mc_step(E,j,coef,m_obj)
        j = clamp(j, 0.0, 1.0) # Ensure physical j
        t += dt
        n += 1

        E_stor[n], j_stor[n], t_stor[n] = E, j, t
        println("E = $E, j = $j, t = $t")
    end

    return 0.07453 .* t_stor[1:n], E_stor[1:n], j_stor[1:n]
end


function run_mc_rp_ra(E0, j0, coef, m_obj; max_step=100, max_t=Inf)

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

    for step in 1:max_step
        j_lc = loss_cone(E)

        if E <= -10000000
            println("Below potential minimum: (E = $E)")
            break
        elseif j <= j_lc
            println("Entered loss cone (j = $j, j_lc = $j_lc)")
            break
        elseif j >= 1.0 || E >= 0.0 || t >= max_t
            break
        end

        E, j, dt = mc_step(E,j,coef,m_obj)
        j = clamp(j, 0.0, 1.0) # Ensure physical j
        t += dt
        n += 1

        E_stor[n], j_stor[n], t_stor[n] = E, j, t
        rc, Lc = find_Lc(E)
        
        rp, ra = find_rp_ra(E, j*Lc)
        rp_stor[n], ra_stor[n] = rp, ra
        println("E = $E, j = $j, t = $t, rp = $rp, ra = $ra")
    end

    return 0.07453 .* t_stor[1:n], E_stor[1:n], j_stor[1:n], rp_stor[1:n], ra_stor[1:n]
end


## Loss cone prescription from Qunbar / Stone ##
function loss_cone(E)
    J_icbo = (4 * G * M_bh) / c
    Jc = last(find_Lc(E))
    return J_icbo / Jc
end

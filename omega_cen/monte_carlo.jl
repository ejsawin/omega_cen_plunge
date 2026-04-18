# --- Monte Carlo Functions --- #

function mc_step(E,j,coef,m_obj)

    # Characterize orbit 
    rc, Lc = find_Lc(E)
    L = j * Lc

    rp, ra = find_rp_ra(E,L)
    period = find_period(E,L)
    a = (ra + rp) / 2
    e = (ra - rp) / (ra + rp)

    # Orbit averaged diff coeffs
    dE_NR = coef.DE1_samp(j,E) + m_obj * coef.DE2_samp(j,E)
    dj_NR = coef.Dj1_samp(j,E) + m_obj * coef.Dj2_samp(j,E)
    dEE_NR = coef.DEE_samp(j,E)
    djj_NR = coef.Djj_samp(j,E)
    dEj_NR = coef.DEj_samp(j,E)

    # GW Emission
    dE_GW, dL_GW = gw_emission(a, e, m_obj)
    dj_GW = dL_GW / Lc

    # Combined NR + GW
    dE_tot = dE_NR + dE_GW
    dj_tot = dj_NR + dj_GW

    # Adaptive timestep (Qunbar & Stone)
    F_safe = 10

    gw_term = E / dE_GW
    nr_term = sqrt((j^2) / (djj_NR * period))

    N_safe = (1/F_safe) * min(gw_term, nr_term) # Choose smaller timestep
    N_orb = max(1.0, N_safe)
    dt = N_orb * period 

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
    dE = dE_tot * N_orb + g1t * sqrt(dEE_NR * N_orb)
    dL = (dj_tot * N_orb + g2t * sqrt(djj_NR * N_orb)) * Lc

    E_step = E + dE
    L_step = L + dL

    # Compute Lc at new E, calculate j
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

        if E <= minimum(Etab)
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

## Loss cone prescription from Qunbar / Stone ##
function loss_cone(E)
    J_icbo = (4 * G * M_bh) / c
    Jc = last(find_Lc(E))
    return J_icbo / Jc
end
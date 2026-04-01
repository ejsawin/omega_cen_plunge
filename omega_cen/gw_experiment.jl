function coef_gw_comp(E,j)
    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit 
    rp, ra = find_rp_ra(E,L)
    period = find_period(E,L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    # Compute GW emissions -> Per orbit and rate
    DE, DL = gw_emission(a,e,m_obj)

    # Compute orbit-averaged diffusion coefficients
    d_E,d_EE,d_L,d_LL,d_EL = find_avg_coef_ej(E,j,m_obj)

    # Ratio of diff coefs with GWs
    e_rat = DE / d_E
    l_rat = DL / d_L

    return DE, (DL/Lc), d_E, d_L, e_rat, l_rat, a, e
end

function gw_step(E, j)
    F_safe = 10

    rc, Lc = find_Lc(E)
    L = j * Lc
    rp, ra = find_rp_ra(E, L)
    period = find_period(E, L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    DE, DL = gw_emission(a, e, m_obj)

    N_safe = (1/F_safe) * abs(E / DE)
    N_time = max(1, N_safe)
    dt = N_time * period

    E_new = E + DE * N_time
    L_new = L + DL * N_time

    # New Lc
    _, Lc_new = find_Lc(E_new)
    j_new = L_new / Lc_new

    return E_new, j_new, dt, a, e
end

function run_mc_gw(E0, j0; max_step=100, max_t=Inf)
    E_stor = Array{Float64}(undef, max_step + 1)
    j_stor = Array{Float64}(undef, max_step + 1)
    t_stor = Array{Float64}(undef, max_step + 1)

    E_stor[1] = E0
    j_stor[1] = j0
    t_stor[1] = 0.0

    E, j, t = E0, j0, 0.0
    n = 1

    for step in 1:max_step
        j_lc = loss_cone(E)
        j_co = 1.0 - 1e-7

        if j <= j_lc || j >= j_co || E >= 0.0 || t >= max_t
            break
        end

        E, j, dt, a, e= gw_step(E, j)
        j = clamp(abs(j), 0.0, 1.0) # Ensure physical j

        t += dt
        n += 1

        E_stor[n] = E
        j_stor[n] = j
        t_stor[n] = t
        
        println("-----")
        println("E = $E, j = $j, dt = $dt, a = $a, e = $e")
        println("-----")

        if E <= minimum(Etab)
            println("Below potential minimum at step $step, breaking.")
            break
        end
    end

    return t_stor[1:n], E_stor[1:n], j_stor[1:n]
end

function loss_cone(E)
    J_icbo = (4 * G * M_bh) / c
    Jc = last(find_Lc(E))
    return J_icbo / Jc
end

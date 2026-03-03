# --- Monte Carlo Functions --- #
using QuadGK

function run_mc(E0, j0; max_step=20, max_t=Inf)

    # Initalize storage arrays
    E_stor = Array{Float64}(undef, max_step+1)
    j_stor = Array{Float64}(undef, max_step+1)
    t_stor = Array{Float64}(undef, max_step+1)

    E_stor[1] = E0
    j_stor[1] = j0
    t_stor[1] = 0.0

    E = E0
    j = j0
    t = 0.0
    n = 1

    for step in 1:max_step 
        j_lc = loss_cone(E)
        j_co = 1.0 - 1e-7

        # If inside loss cone or unbound, break
        if j <= j_lc || j >= j_co || E >= 0.0 || t >= max_t
            break
        end

        # Calculate new (E,j), ensure physical j
        E, j, dt = mc_step(E,j)
        j = clamp(abs(j),0.0,1.0)

        #println(E, j, dt)

        # Increment time, index
        t += dt
        n += 1

        # Save to arrays
        E_stor[n] = E
        j_stor[n] = j
        t_stor[n] = t
    end

    return t_stor[1:n], E_stor[1:n], j_stor[1:n]
end

function mc_step(E,j)

    F_safe = 10 # From Qunbar / Stone

    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit
    rp, ra = find_orbit(E,L)
    period = find_period(E,L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    # Compute orbit averaged diffusion coeffs in (E,j)
    dE_NR = coef.DE_samp(j,E)
    dEE_NR = coef.DEE_samp(j,E)
    dj_NR = coef.Dj_samp(j,E)
    djj_NR = coef.Djj_samp(j,E)
    dEj_NR = coef.DEj_samp(j,E)

    # Compute GW emissions (per orbit); dL -> dj
    dE_GW, dL_GW = gw_emission(a,e,m_obj)
    dj_GW = dL_GW / Lc

    # Combined NR + GW 
    dE_tot = dE_NR + dE_GW
    dj_tot = dj_NR + dj_GW

    # Compute adaptive timestep (Qunbar / Stone)
    e_term = E / dE_GW
    l_term = sqrt((j^2) / (djj_NR * period))

    N_safe = (1/F_safe) * min(e_term,l_term)
    N_time = max(1,N_safe)

    dt = N_time * period

    # Noise 
    g1 = randn()
    g2 = randn()

    # Guard against zero diff coeffs
    denom = sqrt(abs(dEE_NR * djj_NR))

    if denom > 0.0
        rho = dEj_NR / denom 
        rho = clamp(rho,-1.0,1.0) # Ensure safe rho
        g1t = g1
        g2t = g1 * rho + g2 * sqrt(max(0.0,1 - rho^2)) # Float rounding
    else # No diffusion -> Uncorrelated
        rho = 0.0
        g1t = g1
        g2t = g2
    end 

    # Random walk changes in (E,j)
    dE = dE_tot * dt + g1t * sqrt(dEE_NR * dt)
    dj = dj_tot * dt + g2t * sqrt(djj_NR * dt)

    return E + dE, j + dj, dt
end

function gw_step(E,j)

    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit (find rp, ra, period, a, e)
    rp, ra = find_orbit(E,L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    # Compute GW emission rate
    dE_GW, dL_GW = gw_rates(a,e,m_obj)
    dj_GW = dL_GW / Lc

    return dE_GW, dj_GW
end

function gw_timescale(E,j)

    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit (find rp, ra, period, a, e)
    rp, ra = find_orbit(E,L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    # Calculate constants
    c0_1 = a * (1-e^2) / (e^(12/19))
    c0_2 = (1 + (121/304)*e^2)^(-870/2299)
    c0 = c0_1 * c0_2

    beta = (64 * G^3 * M_bh * m_obj * (M_bh + m_obj))/(5*c^5)

    # Integrand -> Eq 5.14
    c1 = (12 * c0^4) / (19 * beta)
    tint = x -> c1*(x^(29/19)*(1+(121/304)*x^2)^(1181/2299))/(1-x^2)^(3/2)

    # Calculate t_gw (integral)
    int_calc, err_calc= quadgk(x -> tint(x), 0, e, rtol=1e-8)
    
    # Return rates
    return int_calc
end

## Loss cone prescription from Qunbar / Stone ##
function loss_cone(E)
    J_icbo = (4 * G * M_bh) / c
    Jc = last(find_Lc(E))
    return J_icbo / Jc
end
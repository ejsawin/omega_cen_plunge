# --- GW Emission Calculations --- #

## GW emission per orbit (E, L) ##
function gw_emission(aa,ee,m_obj)

    E1 = -64 * pi * (G^(7/2)) * (M_bh^2) * (m_obj) * (M_bh + m_obj)^(1/2)
    E2 = 5 * (c^5) * (aa^(7/2)) * (1-(ee^2)) ^ (7/2)
    E3 = 1 + (73/24)*(ee^2) + (37/96)*ee^4

    delE_GW = (E1/E2) * E3

    L1 = -64 * pi * (G^3) * (M_bh^2) * (m_obj)
    L2 = 5 * (c^5) * (aa^2) * (1-(ee^2))^2
    L3 = 1 + (7/8)*(ee^2)

    delL_GW = (L1/L2) * L3

    return delE_GW, delL_GW
end

## GW emission rates (E, L) ##
function gw_rates(aa,ee,m_obj)

    de1 = -32 * (G^4) * (M_bh^2) * (m_obj^2) * (M_bh + m_obj) # Eq 5.4 Peters 1964 
    de2 = 5 * (c^5) * (aa^5) * (1 - (ee^2))^(7/2)
    de3 = 1 + (73/24) * (ee^2) + (37/96) * (ee^4)
    dE_dt = (de1 / de2) * de3

    dl1 = -32 * (G^(7/2)) * (M_bh^2) * (m_obj^2) * (M_bh + m_obj)^(1/2) # Eq 5.5 Peters 1964
    dl2 = 5 * (c^5) * (aa^(7/2)) * (1-(ee^2))^2
    dl3 = 1 + (7/8)*(ee^2)
    dL_dt = (dl1/dl2) * dl3

    return dE_dt, dL_dt
end

## Peters GW Timescale (Integral)##
function gw_timescale(E,j)

    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit (find rp, ra, period, a, e)
    rp, ra = find_rp_ra(E,L)
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
    
    # Return t_GW in Myr
    return 0.07453 * int_calc
end

## Peter GW Timescale Approximations ## 
function peter_approx(E,j)
    # Find Lc, L
    rc, Lc = find_Lc(E)
    L = j * Lc

    # Characterize orbit (find rp, ra, period, a, e)
    rp, ra = find_rp_ra(E,L)
    a = (rp + ra) / 2
    e = (ra - rp) / (ra + rp)

    # Calculate constants
    c0_1 = a * (1-e^2) / (e^(12/19))
    c0_2 = (1 + (121/304)*e^2)^(-870/2299)
    c0 = c0_1 * c0_2

    beta = (64 * G^3 * M_bh * m_obj * (M_bh + m_obj))/(5*c^5)

    # Approximate t_gw
    Tc = (a^4)/(4*beta) # e = 0
    small = (c0^4)/(4*beta)*e^(48/19) # e ~ 0
    large = (768/425)*Tc*(1-e^2)^(7/2) # e ~ 1

    return 0.07453 .* (Tc, small, large) # Convert to Myrs
end
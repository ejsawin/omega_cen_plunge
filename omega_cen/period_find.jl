# --- Orbital period finding using anomaly ---

function period_integrand(tt, _E, _L, rp, ra)
    r_hat = (ra - rp)/2 # Eq 2.18
    t_hat = (ra + rp)/2 # Eq 2.19
    r = r_hat * sin(tt) + t_hat # Eq 2.17
    vr = v_r(r, _E, _L)
    return r_hat * cos(tt) / vr # Eq 2.20
end

function period(_E, _L)
    rp, ra = find_rp_ra(_E, _L)
    return 2 * midpoint(tt -> period_integrand(tt, _E, _L, rp, ra), -pi/2, pi/2, int_steps)
end
# --- GW Emission Calculation --- #
# --- Working (Tested 2/10/2026) --- #

## GW emission per orbit (E, L) ##
function gw_emission(aa,ee,m_obj)

    E1 = -64 * pi * (G^(7/2)) * (M_bh^2) * (m_obj^2) * (M_bh + m_obj)^(1/2)
    E2 = 5 * (c^5) * (aa^(7/2)) * (1-(ee^2)) ^ (7/2)
    E3 = 1 + (73/24)*(ee^2) + (37/96)*ee^4

    delE_GW = (E1/E2) * E3

    L1 = -64 * pi * (G^3) * (M_bh^2) * (m_obj^2)
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
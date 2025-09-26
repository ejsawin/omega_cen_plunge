# --- calculate radial, tangential velocity for (r,E,L)

function v_r(rr,EE,LL) # Eq 1.6
    return sqrt(abs(2*(EE-psi_calc(rr))-(LL^2)/(rr^2)))
end

function v_t(rr,EE,LL) # Eq 1.5
    return LL/rr
end
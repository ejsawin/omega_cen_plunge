# ---- Plummmer model -> Analytic formulas --- #

## Plummer potential ##
function plum_psi(r)
    return -G * M / sqrt(r^2 + b^2)
end

## Plummer DF (E) ##
function plum_df_e(E)
    E <= 0 ? 24√2/(7π^3) * b^2/(G^5 * M^4) * (-E)^(7/2) : 0.0 # Orbit must be bound
end

## Plummer DF (r,v) ##
function plum_df_rv(r, v)
    E = plum_psi(r) + 0.5 * v^2
    return 16π^2 * r^2 * v^2 * plum_df_e(E)
end

## Plummer density ##
function plum_rho(r)
    rho_1 = 3 * M / (4 * pi * b^3)
    rho_2 = (1 + r^2 / b^2) ^ (-5/2)

    return rho_1 * rho_2
end

using SpecialFunctions

## Estimate velocity dispersion from surrounding objects ##
function sigma(r,rad)

    idx = find_shell(r,dat.r) # Find target index
    indices = [idx-rad:idx-1; idx+1:idx+rad]

    tabm = dat.m[indices] 
    tabv = dat.v[indices]

    sigma = (1/sqrt(3))*sqrt(sum(tabm .* (tabv .^ 2)) / sum(tabm))

    return sigma
end

## Find circular velocity at given radius ##
function vc(r)

    idx = find_shell(r,dat.r) # Find target index 
    Menc = M_tab[idx] # Enclosed mass 

    return sqrt((G*Menc)/r) # vc^2 / r = GM(r) / r^2
end

## Finite difference derivative (rvc) ##
function drvc_dr(r)
    rvc(r) = r * vc(r)
    return (rvc(r+epsilon)-rvc(r-epsilon))/(2*epsilon) # Finite difference
end

## Calculate dr/dt from NRR ##
function dr_dt(r,M,rad)

    cou_log = 0.11*length(dat.r) # Coulomb logarithm

    X = vc(r) / (sqrt(2) * sigma(r,rad))
    bracket = erf(X) - (2X/sqrt(pi))*exp(-X^2)

    num = -4 * pi * G^2 * M * rho(r) * cou_log * bracket * r
    denom = vc(r)^2 * drvc_dr(r)

    return num / denom
end
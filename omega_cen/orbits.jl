# --- Orbit bound (rp, ra), velocity (vr,vt), and period (T) finding --- #
# --- Working (Tested 1/26/2026) --- #

## 1D bisection algorithm ##
function bisect(f,r1,r2)

    # Function evaluated at left, right
    f1 = f(r1)
    f2 = f(r2)

    if f1 * f2 > 0 # No root in interval
        return nothing
    end

    for iter in 1:maxiter

        rmid = (r1 + r2) / 2 # Midpoint
        fmid = f(rmid)

        # Found root (within tol)
        if abs(fmid) < tol || (r2 - r1) / 2 < tol
            return rmid
        end

        if f1 * fmid <= 0 # Root to left
            r2 = rmid
            f2 = fmid
        else # Root to right
            r1 = rmid
            f1 = fmid
        end
    end

    return (r1 + r2) / 2
end

## Analytic formulas for rp, ra ##

# rp < r1 -> psi_BH + psi[1] + L^2 / (2 * r ^ 2) = E
function analytic_rp(E,L)
    GM = G * M_bh
    psi1 = psi_tab[1]
    return (GM - sqrt(GM^2 - 2*(psi1-E)*L^2)) / (2*(psi1-E))
end

# ra > rN -> -GMtot/r + L^2 / 2r^2 = E
function analytic_ra(E,L)
    GMtot = G * (sum(dat.m) + M_bh)
    return (GMtot + sqrt(GMtot ^ 2 + 2 * E * L^2))/(-2*E)
end

## Find periapsis, apopsis for given E, L ##
function find_orbit(E,L)

    rc, Lc = find_Lc(E) # Calculate rc, Lc

    # Check for nonphysical / (nearly) circular orbits
    epsL = 1e-8 * Lc # Scalable epsilon

    if L > Lc + epsL # Nonphysical
        error("No physical orbit: L > Lc")
        return nothing

    elseif abs(L - Lc) <= epsL # Circular orbit
        error("Circular orbit: L ~ Lc")
        return nothing

    end

    N = length(dat.r)

    # Psi_eff at leftmost two points
    psi_1 = psi_tot_tab[1] + (L^2)/(2*(dat.r[1]^2))
    psi_2 = psi_tot_tab[2] + (L^2)/(2*(dat.r[2]^2))

    # Psi_eff at rightmost two ponts
    psi_N1 = psi_tot_tab[N-1] + (L^2)/(2*(dat.r[N-1]^2))
    psi_N = psi_tot_tab[N] + (L^2)/(2*(dat.r[N]^2))

    # Left, right indices
    il = 1
    ir = N-1

    # Slope of psi_eff at left, right
    del_psi_l = psi_2 - psi_1
    del_psi_r = psi_N - psi_N1

    if ((del_psi_l < 0) && (del_psi_r > 0)) # Minimum in interval [r1, rN]

        while (ir - il) > 1 # Iterate until adjacent indices

            im = div(il+ir,2) # Midpoint

            # Slope of psi_eff at midpoint
            psi_m = psi_tot_tab[im] + (L^2)/(2*(dat.r[im]^2))
            psi_m1 = psi_tot_tab[im + 1] + (L^2)/(2*(dat.r[im + 1]^2))

            del_psi_m = psi_m1 - psi_m

            if (del_psi_m > 0) # Minimum to left
                ir = im
                del_psi_r = del_psi_m

            elseif (del_psi_m < 0) # Minimum to right
                il = im
                del_psi_l = del_psi_m

            else # Derivative ~ 0 -> Found minimum
                il = im
                ir = im + 1
                break
            end
        end

    elseif del_psi_l > 0 # Minimum to left of grid
        il = 1
        ir = 2

    elseif del_psi_r < 0 # Minimum to right of grid
        il = N - 2
        ir = N - 1
    end

    # Turning points -> E = psi_eff
    f(r) = psi_eff(r,L) - E

    # Initialize variables
    rp = nothing
    ra = nothing

    # Find rp
    fr1 = f(dat.r[1]) # f at 1, il
    fr_il = f(dat.r[il])

    if fr1 * fr_il <= 0 # Sign change -> root within bracket
        rp = bisect(f, dat.r[1], dat.r[il]) # Find via bisection
    end

    if rp === nothing # Default to analytic soln if no root
        rp = analytic_rp(E,L)
    end

    # Find ra
    fr_ir = f(dat.r[ir]) # f at ir, N
    frN = f(dat.r[N])

    if fr_ir * frN <= 0 # Sign change -> Root within bracket
        ra = bisect(f, dat.r[ir], dat.r[N])
    end

    if ra === nothing && rp !== nothing && rp < dat.r[1] # If both rp,ra < r[1]
        fr1 = f(dat.r[1])
        frp = f(rp)
        if frp * fr1 <= 0
            ra = bisect(f, rp, dat.r[1]) # Bisect from rp
        end
    end
    
    if ra === nothing # Default to analytic soln if no root
        ra = analytic_ra(E,L)
    end

    # Return orbit bounds
    return rp, ra
end

## Find radial velocity ##
function find_vr(r,E,L)
    return sqrt(abs(2*(E - psi_calc(r)) - L^2 / r^2)) # Eq 1.6
end

## Find tangential velocity ##
function find_vt(r,E,L)
    return L / r # Eq 1.5
end

## Find orbital period via anomaly ##
function period_int(t,E,L,rp,ra)
    r_hat = (ra - rp) / 2 # Eq 1.18
    t_hat = (ra + rp) / 2 # Eq 1.19

    r = r_hat * sin(t) + t_hat # Eq 1.17
    vr = find_vr(r,E,L) # Radial velocity

    return r_hat * cos(t) / vr # Eq 1.20
end

function find_period(E,L)
    rp, ra = find_orbit(E,L)
    return 2 * midpoint(t -> period_int(t,E,L,rp,ra),-pi/2,pi/2,n_steps) # Eq 1.21
end



















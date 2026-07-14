# --- Orbit Finding (Bounds, Period) --- #

## 1D bisection algorithm ##

#=
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

## Analytic turning point solns ##

# rp < r1: psi_BH + psi[1] + (L^2)/(2*r^2) = E
function anal_rp(E,L)
    GM = G * M_bh
    psi1 = psi_tab[1]
    return L^2 / (GM + sqrt((GM)^2 - 2*(psi1-E)*L^2))
end

# ra > rN: -GMtot/r + (L^2)/(2*r^2) = E
function anal_ra(E,L) 
    GMtot = G * (psi_Mtot + M_bh)
    return (GMtot + sqrt(GMtot ^ 2 + 2 * E * L^2))/(-2*E)
end

## Find effective potential minimum for given L ##
function find_rmin(L)

    # Protect against nonphysical / radial orbits 
    if L < 0
        error("find_rmin: Unphysical orbit, L = $L < 0")
    elseif L == 0
        return nothing
    end

    N = length(psi_rtab)

    # Helper functions
    psi_eff_tab(i) = psi_tot_tab[i] + (L^2) / (2 * psi_rtab[i]^2)
    del_psi(i) = psi_eff_tab(i+1) - psi_eff_tab(i)

    # Check if minimum bracketed 
    if !(del_psi(1) < 0 && del_psi(N-1) > 0)
        #@warn "find_rmin: Minimum off grid for L = $L, Anal soln needed"
        return nothing
    end

    ## Index search for minimum ##

    # Left, right indices
    il = 1
    ir = N - 1 

    while (ir - il) > 1 # Until adjacent indices

        im = div(il + ir, 2) # Midpoint

        if del_psi(im) > 0 # Minimum to left
            ir = im
        elseif del_psi(im) < 0 # Minimum to right 
            il = im
        else
            return psi_rtab[im] # Found minimum!
        end
    end

    return psi_rtab[ir] # Return minimum
end

## Calculate orbital turning points via effective potential ##
function find_rp_ra(E, L)
    Emin = -9999999999
    rc, Lc = find_Lc(E)

    # Ensure physical orbit
    if E >= 0
        #error("find_rp_ra: Unbound orbit, E = $E >= 0")
        return -100, -100
    elseif E < Emin
        #error("find_rp_ra: E = $E below potential minimum Emin = $Emin")
        return -100, -100
    elseif L < 0
        error("find_rp_ra: unphysical input, L = $L < 0")
        return -100, -100
    elseif L >= Lc * (1 - lc_tol)
        #error("find_rp_ra: L = $L >= Lc = $Lc, Circular or Nonphysical Orbit")
        return -100, -100
    end

    f(r) = psi_eff(r, L) - E
    N = length(psi_rtab)

    rmin = find_rmin(L)

    # If no minimum found, default to rc
    if !isnothing(rmin) && f(rmin) > 0
        rmin = rc
    end

    if isnothing(rmin) # Minimum off grid 
        psi_eff_tab(i) = psi_tot_tab[i] + (L^2) / (2 * psi_rtab[i]^2)
        del_psi(i) = psi_eff_tab(i + 1) - psi_eff_tab(i)

        if del_psi(1) > 0
            rp = anal_rp(E, L)
            if f(psi_rtab[N]) > 0
                ra = bisect(f, rc, psi_rtab[N])
            else
                ra = anal_ra(E, L)
            end
        elseif del_psi(N - 1) < 0
            ra = anal_ra(E, L)
            if f(psi_rtab[1]) > 0 && f(rc) < 0
                rp = bisect(f, psi_rtab[1], rc)
            else
                rp = anal_rp(E, L)
            end
        end

    else
        # Find rp
        if f(psi_rtab[1]) < 0
            rp = anal_rp(E, L)
        else
            rp = bisect(f, psi_rtab[1], rmin)
        end

        # Find ra
        if f(psi_rtab[N]) < 0
            ra = anal_ra(E, L)
        else
            ra = bisect(f, rmin, psi_rtab[N])
        end
    end
    return rp, ra
end

## Find radial velocity ##
function find_vr(r,E,L)
    return sqrt(abs(2*(E - psi_calc(r)) - L^2 / r^2)) # Eq 1.16
end

## Find tangential velocity ##
function find_vt(r,E,L)
    return L / r # Eq 1.15
end

## Find orbital period via anomaly ##
function period_int(t,E,L,rp,ra)
    r_hat = (ra - rp) / 2 # Eq 1.28
    t_hat = (ra + rp) / 2 # Eq 1.29

    r = r_hat * sin(t) + t_hat # Eq 1.27
    vr = find_vr(r,E,L) # Radial velocity

    return r_hat * cos(t) / vr # Eq 1.30
end

function find_period(E,L)
    rp, ra = find_rp_ra(E,L) # Eq 1.31
    return 2 * midpoint(t -> period_int(t,E,L,rp,ra),-pi/2,pi/2,n_steps)
end

=#

# --- Orbit Finding (Bounds, Period) --- #

const ORBIT_FAIL = (-100.0, -100.0)

@inline orbit_failed(rp, ra) = (
    !isfinite(rp) ||
    !isfinite(ra) ||
    rp <= 0.0 ||
    ra < rp
)

@inline orbit_circular(rp, ra) = (
    isfinite(rp) &&
    isfinite(ra) &&
    rp > 0.0 &&
    rp == ra
)


# -------------------------------------------------------------------
# "Virtual" extension of the potential outside the original r grid.
#
# Below psi_rtab[1], assume the force is purely Keplerian from the BH.
# The contribution from exterior spherical mass is a constant offset.
# -------------------------------------------------------------------

@inline function inner_psi0()
    r1 = psi_rtab[1]
    GM = G * M_bh

    # Makes psi_inner(r1) exactly match psi_tot_tab[1].
    return psi_tot_tab[1] + GM / r1
end


function psi_orbit(r)   

    if !isfinite(r) || r <= 0
        return NaN
    end

    r1 = psi_rtab[1]
    rN = psi_rtab[end]

    if r < r1
        # Inner extension: constant offset + Keplerian BH potential.
        return inner_psi0() - G * M_bh / r

    elseif r > rN
        # Outer extension, consistent with anal_ra.
        return -G * (psi_Mtot + M_bh) / r

    else
        return psi_calc(r)
    end
end


@inline function psi_eff_orbit(r, L)
    return psi_orbit(r) + L^2 / (2 * r^2)
end


# Circular orbit calculation that uses a Keplerian solution only
# when the circular radius is below the first grid point.
function find_Lc_orbit(E)

    GM = G * M_bh
    Ein = E - inner_psi0()

    # Inner Kepler circular orbit:
    # rc = -GM / (2 Ein)
    # Lc = GM / sqrt(-2 Ein)
    if isfinite(Ein) && Ein < 0
        rc_kep = -GM / (2 * Ein)

        if rc_kep < psi_rtab[1]
            Lc_kep = GM / sqrt(-2 * Ein)
            return rc_kep, Lc_kep
        end
    end

    # Otherwise use your original tabulated-potential routine.
    return find_Lc(E)
end


## 1D bisection algorithm ##
function bisect(f, r1, r2; rtol = 1e-12, atol = 0.0)

    f1 = f(r1)
    f2 = f(r2)

    if !isfinite(f1) || !isfinite(f2)
        return nothing
    end

    if f1 == 0.0
        return r1
    elseif f2 == 0.0
        return r2
    elseif signbit(f1) == signbit(f2)
        return nothing
    end

    for _ in 1:maxiter

        rmid = r1 + (r2 - r1) / 2
        fmid = f(rmid)

        if !isfinite(fmid)
            return nothing
        end

        x_tol = atol + rtol * max(abs(r1), abs(r2))

        if fmid == 0.0 || abs(r2 - r1) / 2 <= x_tol
            return rmid
        end

        if signbit(f1) != signbit(fmid)
            r2 = rmid
            f2 = fmid
        else
            r1 = rmid
            f1 = fmid
        end
    end

    return (r1 + r2) / 2
end

## Analytic turning point solutions ##

# Pericenter in the inner Kepler region.
function anal_rp(E, L)

    GM = G * M_bh
    psi0 = inner_psi0()

    # psi0 - GM/r + L^2/(2r^2) = E
    disc2 = GM^2 - 2 * (psi0 - E) * L^2

    if !isfinite(disc2) || disc2 < 0
        return nothing
    end

    return L^2 / (GM + sqrt(disc2))
end


# Apocenter outside the outer grid.
function anal_ra(E, L)

    if !isfinite(E) || !isfinite(L) || E >= 0 || L < 0
        return nothing
    end

    GMtot = G * (psi_Mtot + M_bh)
    disc2 = GMtot^2 + 2 * E * L^2

    if !isfinite(disc2) || disc2 < 0
        return nothing
    end

    ra = (GMtot + sqrt(disc2)) / (-2 * E)

    return isfinite(ra) && ra > 0 ? ra : nothing
end


## Find effective-potential minimum for given L ##
function find_rmin(L)

    if L < 0 || !isfinite(L)
        return nothing
    elseif L == 0
        return nothing
    end

    N = length(psi_rtab)

    psi_eff_tab(i) = psi_tot_tab[i] + L^2 / (2 * psi_rtab[i]^2)
    del_psi(i) = psi_eff_tab(i + 1) - psi_eff_tab(i)

    # Minimum is not bracketed by the tabulated r grid.
    if !(del_psi(1) < 0 && del_psi(N - 1) > 0)
        return nothing
    end

    il = 1
    ir = N - 1

    while (ir - il) > 1

        im = div(il + ir, 2)

        if del_psi(im) > 0
            ir = im
        elseif del_psi(im) < 0
            il = im
        else
            return psi_rtab[im]
        end
    end

    return psi_rtab[ir]
end


## Calculate orbital turning points via effective potential ##
function find_rp_ra(E, L)

    try
        if !isfinite(E) || !isfinite(L) || E >= 0 || L < 0
            return ORBIT_FAIL
        end

        N = length(psi_rtab)

        rc, Lc = find_Lc_orbit(E)
        
        if !isfinite(rc) || !isfinite(Lc) || rc <= 0.0 || Lc <= 0.0
            return ORBIT_FAIL
        end
        
        # Circularize nearly circular orbits.
        #
        # L > Lc is formally unphysical at fixed E, but allow a tiny overshoot
        # because Lc and the tabulated/interpolated potential are numerical.
        L_circ_min = Lc * (1 - lc_tol)
        L_circ_max = Lc * (1 + lc_tol)
        
        if L > L_circ_max
            return ORBIT_FAIL
        elseif L >= L_circ_min
            return rc, rc
        end

        # Important change:
        # f now works below psi_rtab[1] and above psi_rtab[end].
        f(r) = psi_eff_orbit(r, L) - E

        rmin = find_rmin(L)

        # Keep your original logic. rc is now safe even if it is
        # outside the original potential table.
        if !isnothing(rmin) && f(rmin) > 0
            rmin = rc
        end

        rp = nothing
        ra = nothing

        if isnothing(rmin)

            psi_eff_tab(i) = psi_tot_tab[i] + L^2 / (2 * psi_rtab[i]^2)
            del_psi(i) = psi_eff_tab(i + 1) - psi_eff_tab(i)

            # Effective-potential minimum is below the inner r grid.
            if del_psi(1) > 0

                rp = anal_rp(E, L)

                if f(psi_rtab[N]) > 0
                    # rc may be below psi_rtab[1], but f is now valid there.
                    ra = bisect(f, rc, psi_rtab[N])
                else
                    ra = anal_ra(E, L)
                end

            # Effective-potential minimum is beyond the outer r grid.
            elseif del_psi(N - 1) < 0

                ra = anal_ra(E, L)

                if f(psi_rtab[1]) > 0 && f(rc) < 0
                    # rc may be beyond psi_rtab[end], but f is valid there too.
                    rp = bisect(f, psi_rtab[1], rc)
                else
                    rp = anal_rp(E, L)
                end
            end

        else
            # Find rp
            if f(psi_rtab[1]) < 0
                rp = anal_rp(E, L)
            else
                rp = bisect(f, psi_rtab[1], rmin)
            end

            # Find ra
            if f(psi_rtab[N]) < 0
                ra = anal_ra(E, L)
            else
                ra = bisect(f, rmin, psi_rtab[N])
            end
        end

        if isnothing(rp) || isnothing(ra)
            return ORBIT_FAIL
        end

        if !isfinite(rp) || !isfinite(ra) || rp <= 0 || ra <= rp
            return ORBIT_FAIL
        end

        return rp, ra

    catch err
        err isa InterruptException && rethrow()

        # Uncomment while debugging:
        # @warn "find_rp_ra failed" E L exception=(err, catch_backtrace())

        return ORBIT_FAIL
    end
end


## Find radial velocity ##
function find_vr(r, E, L)

    # Important: uses psi_orbit, not psi_calc directly.
    vr2 = 2 * (E - psi_orbit(r)) - L^2 / r^2

    # Retains your original behavior near turning points.
    return sqrt(abs(vr2))
end


## Find tangential velocity ##
function find_vt(r, E, L)
    return L / r
end


## Find orbital period via anomaly ##
function period_int(t, E, L, rp, ra)

    r_hat = (ra - rp) / 2
    t_hat = (ra + rp) / 2

    r = r_hat * sin(t) + t_hat
    vr = find_vr(r, E, L)

    return r_hat * cos(t) / vr
end


function find_period(E, L)

    rp, ra = find_rp_ra(E, L)

    if orbit_failed(rp, ra)
        return -100.0
    end

    # Circularized orbit: rp == ra == rc.
    #
    # Use Lc, not input L, since a near-circular input orbit is being
    # deliberately replaced by the exact circular orbit at energy E.
    if orbit_circular(rp, ra)
        rc, Lc = find_Lc_orbit(E)

        if !isfinite(rc) || !isfinite(Lc) || rc <= 0.0 || Lc <= 0.0
            return -100.0
        end

        return 2 * pi * rc^2 / Lc
    end

    try
        return 2 * midpoint(
            t -> period_int(t, E, L, rp, ra),
            -pi / 2,
            pi / 2,
            n_steps,
        )
    catch
        return -100.0
    end
end

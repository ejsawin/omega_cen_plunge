# --- Orbit Finding (Bounds, Period) --- #


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

## Analytic solns ##
function anal_rp(E,L)
    GM = G * M_bh
    psi1 = psi_tab[1]
    return L^2 / (GM + sqrt((GM)^2 - 2*(psi1-E)*L^2))
end

function anal_ra(E,L) 
    GMtot = G * (sum(dat.m) + M_bh)
    return (GMtot + sqrt(GMtot ^ 2 + 2 * E * L^2))/(-2*E)
end

## Find psi_eff minimum for given L ##
function find_rmin(L)

    # Protect against nonphysical / radial orbits 
    if L < 0
        error("find_rmin: Unphysical orbit, L = $L < 0")
    elseif L == 0
        #@warn "find_rmin: Radial orbit, Anal soln needed"
        return nothing
    end

    N = length(dat.r)

    # Helper functions
    psi_eff_tab(i) = psi_tot_tab[i] + (L^2) / (2 * dat.r[i]^2)
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
            return dat.r[im] # Found minimum!
        end
    end

    return dat.r[ir] # Return minimum
end

## Find orbit turning points ##
function find_rp_ra(E,L)

    Emin = psi_tot_tab[1] # Potential minimum 
    rc, Lc = find_Lc(E)

    # Check if orbit physical
    if E >= 0 # Unbound
        error("find_rp_ra: Unbound orbit, E = $E >= 0")
    elseif E < Emin # Below potential minimum 
        error("find_rp_ra: E = $E below potential minimum Emin = $Emin")
    elseif L < 0
        error("find_rp_ra: unphysical input, L = $L < 0")
    elseif L >= Lc * (1- lc_tol) # L => Lc
        error("find_rp_ra: L = $L >= Lc = $Lc, Circular or Nonphysical Orbit")
    end

    f(r) = psi_eff(r,L) - E # Zeros <-> rp, ra
    N = length(dat.r)

    # Find psi_eff minimum
    rmin = find_rmin(L)

    # If invalid rmin, fall back to rc
    if !isnothing(rmin) && f(rmin) > 0
        rmin = rc
    end
    
    if isnothing(rmin) # Minimum off grid

        psi_eff_tab(i) = psi_tot_tab[i] + (L^2) / (2*dat.r[i]^2)
        del_psi(i) = psi_eff_tab(i+1) - psi_eff_tab(i)

        if del_psi(1) > 0 # Minimum to left of grid -> Anal rp

            rp = anal_rp(E,L)

            if f(dat.r[N]) > 0 # ra on grid -> bisect
                ra = bisect(f, dat.r[1], dat.r[N])
            else # ra off grid -> anal soln
                ra = anal_ra(E,L)
            end

        elseif del_psi(N-1) < 0 # Minimum to right of grid -> Anal ra

            ra = anal_ra(E, L)

            if f(dat.r[1]) > 0 && f(rc) < 0 # rp bracketed by [r1, rc]
                rp = bisect(f, dat.r[1], rc)
            else
                rp = anal_rp(E, L)
            end
        end
    else # Minimum on grid 
        
        # Find rp 
        if f(dat.r[1]) < 0 # rp < r1 -> Anal soln
            rp = anal_rp(E,L)
        else # rp on grid -> Bisect
            rp = bisect(f,dat.r[1],rmin)
        end

        # Find ra
        if f(dat.r[N]) < 0 # ra > rN -> Anal soln
            ra = anal_ra(E,L)
        else # ra on grid  -> Bisect
            ra = bisect(f,rmin,dat.r[N])
        end
    end    

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
    rp, ra = find_rp_ra(E,L)
    return 2 * midpoint(t -> period_int(t,E,L,rp,ra),-pi/2,pi/2,n_steps) # Eq 1.21
end

function find_root(f,r1,r2) # Simple bisection

    f1 = f(r1)
    f2 = f(r2)

    if f1 * f2 > 0
        return nothing
    end

    for iter in 1:maxiter
        rmid = (r1 + r2) / 2
        fmid = f(rmid)

        if abs(fmid) < tol || (r2 - r1) / 2 < tol
            return rmid
        end

        if f1 * fmid <= 0
            r2 = rmid
            f2 = fmid
        else
            r1 = rmid
            f1 = fmid
        end
    end

    return(r1 + r2) / 2
end

function analytic_rp(E,L) # Only for rp < r1
    return (G*M_bh - sqrt((G*M_bh)^ 2 - 2*(psi_tab[1] - E)*(L^2)))/(2*(psi_tab[1]-E))
end

function analytic_ra(E,L) # Only for ra > rN
    return (G*M_tot + sqrt((G*M_tot)^2 + 2*E*(L^2)))/(-2*E)
end


function find_rp_ra(E,L)

    rc, Lc = find_Lc(E)

    f(r)=psi_eff(r,L) - E # zeros <-> orbit bounds

    # look for rp in [rmin, rc], default to analytic if none found
    rp = find_root(f,rmin,rc)

    if rp === nothing
        rp = analytic_rp(E,L)
    end

    # same for ra
    ra=find_root(f,rc,rmax)
        
    if ra === nothing
        ra = analytic_ra(E,L)
    end

    return rp, ra
end
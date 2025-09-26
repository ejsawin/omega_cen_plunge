# --- Find rp, ra --- 

function compute_rp_ra(E, L)
    N = length(r_tab)

    # --- standard method ---
    psi_eff_1  = psi_tab[1] + (L^2) / (2 * (r_tab[1])^2)
    psi_eff_2  = psi_tab[2] + (L^2) / (2 * (r_tab[2])^2)
    psi_eff_N1 = psi_tab[N-1] + (L^2) / (2 * (r_tab[N-1])^2)
    psi_eff_N  = psi_tab[N]   + (L^2) / (2 * (r_tab[N])^2)

    il, ir = 1, N-1
    delta_psi_eff_l = psi_eff_2 - psi_eff_1
    delta_psi_eff_r = psi_eff_N - psi_eff_N1

    if ((delta_psi_eff_l < 0) && (delta_psi_eff_r > 0))
        while (ir - il > 1)
            im = div(il+ir, 2)
            psi_eff_m  = psi_tab[im]   + (L^2) / (2 * (r_tab[im])^2)
            psi_eff_m1 = psi_tab[im+1] + (L^2) / (2 * (r_tab[im+1])^2)
            delta_psi_eff_m = psi_eff_m1 - psi_eff_m
            if delta_psi_eff_m < 0
                il = im
            elseif delta_psi_eff_m > 0
                ir = im
            else
                il, ir = im, im+1
                break
            end
        end
    end

    function bisect(i,j,E)
        psi_eff_l = psi_tab[i] + (L^2)/(2*(r_tab[i]^2)) - E
        psi_eff_r = psi_tab[j] + (L^2)/(2*(r_tab[j]^2)) - E
        while j - i > 1
            im = div(i+j,2)
            fm = psi_tab[im] + (L^2)/(2*(r_tab[im]^2)) - E
            if psi_eff_l*fm <= 0
                j, psi_eff_r = im, fm
            else
                i, psi_eff_l = im, fm
            end
        end
        return r_tab[i] - psi_eff_l * (r_tab[j]-r_tab[i]) / (psi_eff_r-psi_eff_l)
    end

    rp = bisect(1,il,E)
    ra = bisect(ir,N,E)

    # --- local cubic approximation if rp < r_tab[1] ---
    
    if rp < r_tab[1]

        # use first few points to build expansion
        drr = r_tab[2] - r_tab[1]
        psi0 = psi_tab[1]
        d1 = (psi_tab[2] - psi_tab[1]) / drr
        d2 = (psi_tab[3] - 2*psi_tab[2] + psi_tab[1]) / drr^2
        d3 = (psi_tab[4] - 3*psi_tab[3] + 3*psi_tab[2] - psi_tab[1]) / drr^3

        function psi_eff_local(r) # psi eff approximation 
            psi_ = psi0 + d1*r + 0.5*d2*r^2 + (1/6)*d3*r^3
            return psi_ + L^2/(2*r^2)
        end

        # Bracket root between [eps, r_tab[2]]
        eps = 1e-8
        rlo, rhi = eps, r_tab[2]
        flo, fhi = psi_eff_local(rlo)-E, psi_eff_local(rhi)-E

        if flo*fhi > 0
            rp = r_tab[1]   # clamp if no root found
        else
            # bisection 
            for _ in 1:30
                mid = 0.5*(rlo+rhi)
                fm = psi_eff_local(mid)-E
                if flo*fm <= 0
                    rhi, fhi = mid, fm
                else
                    rlo, flo = mid, fm
                end
            end
            rp = 0.5*(rlo+rhi)
        end
    end

    return rp, ra
end

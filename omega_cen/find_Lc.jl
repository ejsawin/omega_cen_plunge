function find_Lc(E)
    # leftmost two points
    eta_1 = 2*r_tab[1]^2 * (E - psi_tot_tab[1])
    eta_2 = 2*r_tab[2]^2 * (E - psi_tot_tab[2])

    # rightmost two points 
    eta_N1 = 2*r_tab[N-1]^2 * (E - psi_tot_tab[N-1])
    eta_N = 2*r_tab[N]^2 * (E - psi_tot_tab[N])

    # left, right indices
    il = 1
    ir = N-1

    # change in eta
    delta_eta_l = eta_2 - eta_1 
    delta_eta_r = eta_N - eta_N1

    if ((delta_eta_l > 0) && (delta_eta_r < 0)) # i.e. maximum in [r1,rN]

        while (ir - il) > 1 # until adjacent indices

            im = div(il+ir,2) # midpoint

            eta_m = 2*r_tab[im]^2 * (E - psi_tot_tab[im])
            eta_m1 = 2*r_tab[im+1]^2 * (E - psi_tot_tab[im+1])
            delta_eta_m = eta_m1 - eta_m 
            
            if (delta_eta_m > 0) # maximum to right 
                il = im
                delta_eta_l = delta_eta_m

            elseif (delta_eta_m < 0)
                ir = im
                delta_eta_r = delta_eta_m
                
            else # derivative ~ 0 -> found bracket
                il = im
                ir = im + 1
                break
            end
        end
        
        # [il, il + 1]
        rl = r_tab[il]
        rl1 = r_tab[il+1]

        psil = psi_tot_tab[il]
        psil1 = psi_tot_tab[il+1]

        slopel = (psil1-psil)/(1/rl1-1/rl)

        rsol_l = 0.5 * 1/((E-psil)/slopel + 1/rl)

        # [ir, ir + 1]
        rr = r_tab[ir]
        rr1 = r_tab[ir+1]

        psir = psi_tot_tab[ir]
        psir1 = psi_tot_tab[ir+1]

        sloper = (psir1-psir)/(1/rr1-1/rr)
        
        rsol_r = 0.5 * 1/((E-psir)/sloper + 1/rr)

        if (rl <= rsol_l <= rl1)
            rsol = rsol_l
            psi = psil + slopel * (1/rsol - 1/rl)
        else 
            rsol = rsol_r
            psi = psir + sloper * (1/rsol - 1/rr)
        end

        Lc2 = 2*rsol^2*(E-psi)

    elseif (delta_eta_l <= 0) # rsol < r_tab[1]
    
        if eta_2 < eta_1 # decreasing eta -> r < r1 -> analytic soln
            rsol = -G * M_bh / (2*(E - psi_tab[1]))
            psi = psi_tab[1] - G * M_bh / rsol

        else # maximum in [r1,r2]
            r1 = r_tab[1]
            r2 = r_tab[2]

            psi1 = psi_tot_tab[1]
            psi2 = psi_tot_tab[2]

            slope_12 = (psi2 - psi1) / (1/r2 - 1/r1)

            rsol = 0.5 * 1 / ((E-psi1)/slope_12 + 1/r1) # lin interp
            psi = psi1 + slope_12*(1/rsol - 1/r1)
        end

        Lc2 = 2*rsol^2*(E-psi)

    else # if (delta_eta_r > 0)
        
        if eta_N > eta_N1 # increasing eta -> r > rN -> analytic soln
            rsol = -(G*M_tot)/(2*E)
            psi = -(G*M_tot)/rsol

        else # maximum in [rN-1, rN]
            rN1 = r_tab[N-1]
            rN = r_tab[N]

            psi_N1=psi_tab_tot[N-1]
            psi_N=psi_tab_tot[N]

            slope_N1N = (psi_N - psi_N1) / (1/rN - 1/rN1)
            
            rsol = 0.5 * 1 / ((E - psi_N1)/slope_N1N + 1/r_N1) # lin interp
            psi = psi_N1 + slope_N1N * (1/rsol - 1/rN1)

        end

        Lc2 = 2*rsol^2*(E-psi)
    end 

    return rsol, sqrt(Lc2) # return r_c, L_c
end

function dLc_dE(E)
    return (find_Lc(E+epsilon)[2]-find_Lc(E-epsilon)[2])/(2*epsilon)
end

function d2Lc_dE2(E)
    return (find_Lc(E+epsilon)[2]+find_Lc(E-epsilon)[2]-2*find_Lc(E)[2])/(epsilon^2)
end

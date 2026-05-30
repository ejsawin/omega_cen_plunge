# --- Circular angular momentum, derivatives for given E --- #

## Calculate rc, Lc via effective potential ##
function find_Lc(E)

    N = length(dat.r)

    # psi(r) + (L^2)/(2*r^2) = E -> L^2 = eta = 2*r^2*(E-psi(r))

    # Value at leftmost two points
    eta_1 = 2*dat.r[1]^2 * (E - psi_tot_tab[1])
    eta_2 = 2*dat.r[2]^2 * (E - psi_tot_tab[2])

    # Value at rightmost two points
    eta_N1 = 2*dat.r[N-1]^2 * (E - psi_tot_tab[N-1])
    eta_N = 2*dat.r[N]^2 * (E - psi_tot_tab[N])

    # Starting left, right indices
    il = 1
    ir = N-1

    # Slope of eta at left, right
    del_eta_l = eta_2 - eta_1
    del_eta_r = eta_N - eta_N1

    if ((del_eta_l > 0) && (del_eta_r < 0)) # I.e. maximum in [r1,rN]

        while (ir - il) > 1 # Iterate until adjacent indices

            # Calculate eta, slope at midpoint
            im = div(il + ir, 2)

            eta_m = 2*dat.r[im]^2 * (E - psi_tot_tab[im])
            eta_m1 = 2*dat.r[im+1]^2 * (E - psi_tot_tab[im+1])

            del_eta_m = eta_m1 - eta_m

            if (del_eta_m > 0) # Maximum to right
                il = im
                del_eta_l = del_eta_m

            elseif (del_eta_m < 0) # Maximum to left
                ir = im
                del_eta_r = del_eta_m

            else # Derivative ~0 -> Found maximum
                il = im
                ir = im + 1
                break
            end
        end

        # Linear interpolation (left bracket)
        rl = dat.r[il]
        rl1 = dat.r[il+1]

        psil = psi_tot_tab[il]
        psil1 = psi_tot_tab[il+1]

        slopel = (psil1-psil)/(1/rl1-1/rl)

        rsol_l = 0.5 * 1/((E-psil)/slopel + 1/rl)

        # Linear interpolation (right bracket)
        rr = dat.r[ir]
        rr1 = dat.r[ir+1]

        psir = psi_tot_tab[ir]
        psir1 = psi_tot_tab[ir+1]

        sloper = (psir1-psir)/(1/rr1-1/rr)

        rsol_r = 0.5 * 1/((E-psir)/sloper + 1/rr)

        if (rl <= rsol_l <= rl1) # Soln in left bracket
            rsol = rsol_l
            psi = psil + slopel * (1/rsol - 1/rl)

        else # Soln in right bracket
            rsol = rsol_r
            psi = psir + sloper * (1/rsol - 1/rr)
        end

        Lc2 = max(0.0,2*rsol^2*(E-psi)) # Protect against nonphysical solns

    elseif (del_eta_l <= 0) # rsol < r_tab[1]

        if eta_2 < eta_1 # rc < r1 -> Analytic soln
            rsol = -G * M_bh / (2 * (E - psi_tab[1]))
            psi = psi_tab[1] - G * M_bh / rsol

        else # Maximum in [r1, r2] -> Linear interpolation
            r1 = dat.r[1]
            r2 = dat.r[2]

            psi1 = psi_tot_tab[1]
            psi2 = psi_tot_tab[2]

            slope_12 = (psi2 - psi1) / (1/r2 - 1/r1)

            rsol = 0.5 * 1 / ((E-psi1)/slope_12 + 1/r1)
            psi = psi1 + slope_12*(1/rsol - 1/r1)
        end

        Lc2 = max(0.0,2*rsol^2*(E-psi)) # Protect against nonphysical solns

    else # if del_eta_r > 0

        if eta_N > eta_N1 # rc > rN -> Analytic soln
            Mtot = sum(dat.m) + M_bh # total mass
            rsol = -(G*Mtot)/(2*E)
            psi = -(G*Mtot)/rsol

        else # Maximum in [rN-1,rN] -> Linear interpolation
            rN1 = dat.r[N-1]
            rN = dat.r[N]

            psi_N1 = psi_tot_tab[N-1]
            psi_N = psi_tot_tab[N]

            slope_N1N = (psi_N - psi_N1) / (1/rN - 1/rN1)

            rsol = 0.5 * 1 / ((E - psi_N1)/slope_N1N + 1/rN1)
            psi = psi_N1 + slope_N1N * (1/rsol - 1/rN1)
        end

         Lc2 = max(0.0,2*rsol^2*(E-psi)) # Protect against nonphysical solns
    end

    return rsol, sqrt(Lc2) # rc, Lc
end

## Finite difference derivatives for (E,L) -> (E,j) ##
function dLc_dE(E)
    return (find_Lc(E+epsilon)[2]-find_Lc(E-epsilon)[2])/(2*epsilon)
end

function d2Lc_dE2(E)
    return (find_Lc(E+epsilon)[2]+find_Lc(E-epsilon)[2]-2*find_Lc(E)[2])/(epsilon^2)
end



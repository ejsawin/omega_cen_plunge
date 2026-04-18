# --- Gravitational Potential Calculation --- #
# --- Working (Tested 1/22/2026) --- #

## Find shell index for given radius ##
function find_shell(rr,r_tab)
    N = length(r_tab)

    # Check edge cases
    if rr < r_tab[1] # First shell
        return 0

    elseif rr == r_tab[N] # Last shell
        return N-1

    elseif rr > r_tab[N] # Out of bounds
        return N
    end
    
    # Simple bisection
    left, right = 1, N
    while left < right
        mid = div(left + right, 2)
        if r_tab[mid] <= rr
            left = mid + 1
        else
            right = mid
        end
    end

    return left - 1  # Largest index r_tab[i] <= rr
end

## Generate psi, M arrays via Henon- field only! ##
function find_psi_arrays(r_tab,m_tab)
    N = length(r_tab)

    psi_arr=zeros(N+1) # Initialize arrays
    M_arr=zeros(N+1)
    psi_tot_arr=zeros(N)

    psi_arr[end] = 0.0 # Eq 1.17
    M_arr[end-1] = sum(m_tab) # Eq 1.18
    psi_arr[end-1] = -G*M_arr[N] / r_tab[N]

    @inbounds for k in (N-1):-1:1 # Iterate over all shells
        M_arr[k] = M_arr[k+1] - m_tab[k+1] # Eq 1.20
        psi_arr[k] = psi_arr[k+1] - G*M_arr[k] * (1/r_tab[k] - 1/r_tab[k+1]) # Eq 1.19
    end

    # Total potential (including IMBH)
    psi_tot_arr = psi_arr[1:end-1] - G * M_bh ./ r_tab

    # Return arrays
    return psi_arr, psi_tot_arr, M_arr
end

## Find exact psi ##
function find_psi_exact(rr,r_tab,psi_arr,M_arr)
    N = length(r_tab)

    k = find_shell(rr,r_tab) # Find radial shell
    imbh_term = -G*M_bh / rr # IMBH contribution, Eq 1.14

    if k == 0 # Inside first shell
        return psi_arr[1] + imbh_term

    elseif k == N # Last shell
        return -G * M_arr[N] / rr + imbh_term

    else # Linear interpolation
        inv_rk = 1 / r_tab[k]
        inv_rk1 = 1 / r_tab[k+1]
        inv_rr = 1 / rr

        frac = (inv_rk - inv_rr) / (inv_rk - inv_rk1)
        return psi_arr[k] + frac * (psi_arr[k+1] - psi_arr[k]) + imbh_term # Eq 1.21, 1.25
    end
end

## Wrapper function for psi ##
psi_exact(r_tab, psi_arr, M_arr) = rr -> find_psi_exact(rr, r_tab, psi_arr, M_arr)

## Effective potential ##
function psi_eff(r,L)
    N = length(dat.r)
    Mtot = sum(dat.m) + M_bh

    if r < dat.r[1] # r < r1 -> Analytic
        return -G*M_bh/r + psi_tab[1] + (L^2)/(2*r^2)

    elseif r > dat.r[N] # r > rN -> Keplerian
        return -(G*Mtot)/r + (L^2)/(2*r^2)

    else # r1 < r < rN
        return psi_calc(r) + (L^2)/(2*r^2)
    end
end
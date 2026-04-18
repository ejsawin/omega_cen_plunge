# --- Diffusion coefficient calculations --- #

## Calculate velocity diff coeffs using log sampling ##
function vel_coef(r,E,L)

    # Object total velocity
    vel = sqrt(find_vr(r,E,L) ^ 2 + find_vt(r,E,L) ^ 2)

    logr = log(r) # Convert into log space
    log_vel = log(vel)

    # Coulomb logarithm
    cou_log = log(M_bh / mean(dat.m))

    # Constants
    c1 = (-G^2 * cou_log) / (vel ^ 2 * r ^ 3)
    c2 = (2 * G^2 * cou_log) / (3 * r ^ 3)

    # Integrands (DF in log space!)
    int1(logv) = DF.F(logr,logv)
    int2(logv) = DF.N(logr,logv)
    int3(logv) = (v = exp(logv); (v^2 / vel^3) * DF.F(logr, logv))
    int4(logv) = (v = exp(logv); DF.F(logr, logv) / v)

    # Compute integrals via midpoint
    calc1 = midpoint(x->int1(x),log(dat.vmin),log_vel,n_steps)
    calc2 = midpoint(x->int2(x),log(dat.vmin),log_vel,n_steps)
    calc3 = midpoint(x->int3(x),log(dat.vmin),log_vel,n_steps)
    calc4 = midpoint(x->int4(x),log_vel,log(dat.vmax),n_steps)

    # Drift coefficient -> Split by mass dependence 
    vp1 = c1 * calc1 # Eq 1.52
    vp2 = c1 * calc2 # Mass dependant 

    # Diffusion coefficients -> No mass dependence
    vp_sq = c2 * (calc3 + calc4) # Eq 1.53
    vt_sq = c2 * ((3 / vel) * calc1 - calc3 + 2 * calc4) # Eq 1.54

    return vp1, vp2, vp_sq, vt_sq
end

## Calculate local (E,L) diff coeffs ##
function loc_coef(r,E,L)

    # Object total velocity
    vel = sqrt(find_vr(r,E,L) ^ 2 + find_vt(r,E,L) ^ 2)

    logr = log(r) # Convert into log space
    log_vel = log(vel)

    # Vel diff coefs
    vp1, vp2, vp_sq, vt_sq = vel_coef(r,E,L)

    # 1st order coefs split, 2nd order unchanged
    delE1 = 0.5 * (vp_sq + vt_sq) + vel*vp1 # Eq 1.55
    delE2 = vel * vp2 # Mass dependent 

    delL1 = (L / vel) * vp1 + (r ^ 2) / (4 * L) * vt_sq # Eq 1.57
    delL2 = (L / vel) * vp2 # Mass dependent

    delE_sq = vel ^ 2 * vp_sq # Eq 1.56
    delL_sq = (L ^ 2 / vel ^ 2) * vp_sq + 0.5*(r^2 - L^2 / vel^2) * vt_sq # Eq 1.58
    delEL = L * vp_sq # Eq 1.59

    return delE1, delE2, delL1, delL2, delE_sq, delL_sq, delEL
end

## Calculate orbit averaged (E,L) diff coeffs ##
function avg_coef_el(E,L)
    
    # Calculate orbit turning points, period
    rp, ra = find_rp_ra(E,L)
    T = find_period(E,L)
    h = pi/n_steps

    # Initialize variables 
    DE1 = DE2 = DL1 = DL2 = DEE = DLL = DEL = 0.0

    # Simultaneous midpoint integration over orbit #
    for i in 1:n_steps
        
        # Anomaly conversion
        tt = -pi/2 + (i - 0.5) * h
        r_hat = (ra - rp) / 2 # Eq 1.28
        t_hat = (ra + rp) / 2 # Eq 1.29
        rr = r_hat * sin(tt) + t_hat # Eq 1.27

        pref = (2h * r_hat * cos(tt)) / (T * find_vr(rr,E,L)) # Prefactor

        # Compute local diff coeffs 
        dE1, dE2, dL1, dL2, dEE, dLL, dEL = loc_coef(rr,E,L)

        DE1 += pref * dE1 # Eq 1.65
        DE2 += pref * dE2 # Mass dependent

        DL1 += pref * dL1 # Eq 1.67
        DL2 += pref * dL2 # Mass dependent 

        DEE += pref * dEE # Eq 1.66
        DLL += pref * dLL # Eq 1.68
        DEL += pref * dEL # Eq 1.69
    end

    return DE1, DE2, DL1, DL2, DEE, DLL, DEL
end

## Convert (E,L) -> (E,j) coeffs ##
function avg_coef_ej(E,j)

    # Find Jc, derivatives; J
    Jc = find_Lc(E)[2]
    dJc = dLc_dE(E)
    d2Jc = d2Lc_dE2(E)
    J = j * Jc

    # Compute avg diff coeff in (E,L)
    DE1_L, DE2_L, DL1_L, DL2_L, DEE_L, DLL_L, DEL_L = avg_coef_el(E,J)

    # Energy coeffs preserved under transformation
    DE1, DE2, DEE = DE1_L, DE2_L, DEE_L

    # Convert j coeffs -> Eq 1.73, 1.74
    Dj1_1 = (DL1_L/Jc) - (J/(Jc^2))*dJc*DE1_L - (1/(Jc^2))*dJc*DEL_L
    Dj1_2 = -(J/2)*((d2Jc/(Jc^2)) - (2*(dJc^2))/(Jc^3))*DEE_L

    Dj1 = Dj1_1 + Dj1_2
    Dj2 = (DL2_L/Jc) - (J/(Jc^2))*dJc*DE2_L # Mass dependent

    # Diffusion -> No mass dependence 
    Djj = ((J*dJc)/(Jc^2))^2 * DEE_L - 2 * ((J*dJc/(Jc^3))) * DEL_L + (DLL_L/(Jc^2)) # Eq 1.76
    DEj = -(J/(Jc^2))*dJc*DEE_L + (DEL_L / Jc) # Eq 1.75

    return DE1, DE2, Dj1, Dj2, DEE, Djj, DEj
end

## Compute diff coeff grid (E,j) for sampling ##
function coef_grid(E_tab,j_tab)

    # Initialize matrices
    DE1_tab = zeros(res,res)
    DE2_tab = zeros(res,res)
    Dj1_tab = zeros(res,res)
    Dj2_tab = zeros(res,res)
    DEE_tab = zeros(res,res)
    Djj_tab = zeros(res,res)
    DEj_tab = zeros(res,res)

    # Loop over energy, momentum
    Threads.@threads for i in 1:res
        E = E_tab[i]
        for k in 1:res
            j = j_tab[k]
            DE1, DE2, Dj1, Dj2, DEE, Djj, DEj = avg_coef_ej(E,j)
            @inbounds begin
                DE1_tab[k,i] = DE1
                DE2_tab[k,i] = DE2
                Dj1_tab[k,i] = Dj1
                Dj2_tab[k,i] = Dj2
                DEE_tab[k,i] = DEE
                Djj_tab[k,i] = Djj
                DEj_tab[k,i] = DEj
            end
        end
    end
    return DE1_tab, DE2_tab, Dj1_tab, Dj2_tab, DEE_tab, Djj_tab, DEj_tab
end

struct DiffusionCoeffs

    # Orbit-averaged diff coeff grids
    DE1_tab::Matrix{Float64}
    DE2_tab::Matrix{Float64} # Mass dependent 
    Dj1_tab::Matrix{Float64}
    Dj2_tab::Matrix{Float64} # Mass dependent
    DEE_tab::Matrix{Float64}
    Djj_tab::Matrix{Float64}
    DEj_tab::Matrix{Float64}

    # Calculated energy and bounds (quantile)
    E_tab::AbstractVector{Float64}
    j_tab::AbstractVector{Float64}
    emin::Float64
    emax::Float64

    # Interpolated functions 
    DE1_samp::Function
    DE2_samp::Function # Mass dependent 
    Dj1_samp::Function
    Dj2_samp::Function # Mass dependent
    DEE_samp::Function
    Djj_samp::Function
    DEj_samp::Function
end

## Generate diff coeffs, interpolated functions ##
function generate_coeffs()

    # Calculate energy, quantile for snapshot 
    EE = 0.5 .* (dat.vr.^2 .+ dat.vt.^2) .+ psi_calc.(dat.r)

    emin = quantile(EE,1e-5)
    emax = quantile(EE,0.99999)

    # Sampling grid for E,j 
    E_tab = range(emin,emax,res)
    j_tab = range(1e-5,0.99999,res)

    # Compute coeffs over grid 
    (DE1_tab, DE2_tab, Dj1_tab, Dj2_tab, 
     DEE_tab, Djj_tab, DEj_tab) = coef_grid(E_tab, j_tab)

    # Bilinear interpolation
    interp(Z) = (j,E) -> bilinear_interp(j, E, E_tab, j_tab, Z)
    
    return DiffusionCoeffs(
        DE1_tab, DE2_tab, Dj1_tab, Dj2_tab,
        DEE_tab, Djj_tab, DEj_tab,
        E_tab, j_tab, emin, emax,
        interp(DE1_tab), interp(DE2_tab),
        interp(Dj1_tab), interp(Dj2_tab),
        interp(DEE_tab), interp(Djj_tab),
        interp(DEj_tab))
end

## Save orbit-averaged diff coeff struct to HDF5 ##
function save_coeffs(name::String, dc::DiffusionCoeffs)
    h5open(name, "w") do file
        file["E_tab"] = collect(dc.E_tab)
        file["j_tab"] = collect(dc.j_tab)
        file["emin"] = dc.emin
        file["emax"] = dc.emax
        file["DE1_tab"] = dc.DE1_tab
        file["DE2_tab"] = dc.DE2_tab
        file["Dj1_tab"] = dc.Dj1_tab
        file["Dj2_tab"] = dc.Dj2_tab
        file["DEE_tab"] = dc.DEE_tab
        file["Djj_tab"] = dc.Djj_tab
        file["DEj_tab"] = dc.DEj_tab
    end
end

## Load diff coeff HDF5 ##
function load_coeffs(name::String)
    h5open(name, "r") do file
        # Read in coeff data
        E_tab = read(file["E_tab"])
        j_tab = read(file["j_tab"])
        emin = read(file["emin"])
        emax = read(file["emax"])
        DE1_tab = read(file["DE1_tab"])
        DE2_tab = read(file["DE2_tab"])
        Dj1_tab = read(file["Dj1_tab"])
        Dj2_tab = read(file["Dj2_tab"])
        DEE_tab = read(file["DEE_tab"])
        Djj_tab = read(file["Djj_tab"])
        DEj_tab = read(file["DEj_tab"])

        # Bilinear interpolation
        interp(Z) = (j,E) -> bilinear_interp(j, E, j_tab, E_tab, Z)

        return DiffusionCoeffs(
            DE1_tab, DE2_tab, Dj1_tab, Dj2_tab,
            DEE_tab, Djj_tab, DEj_tab,
            E_tab, j_tab, emin, emax,
            interp(DE1_tab), interp(DE2_tab),
            interp(Dj1_tab), interp(Dj2_tab),
            interp(DEE_tab), interp(Djj_tab),
            interp(DEj_tab))
    end
end
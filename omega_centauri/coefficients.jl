# --- Diffusion coefficient calculations --- #
using HDF5

## Compute velocity diffusion coefficients for given (r,E,L,m) ##
function find_coef_lin(r,E,L,m_obj) # Linear sampling

    # Object total velocity
    vel = sqrt(find_vr(r,E,L) ^ 2 + find_vt(r,E,L) ^ 2)

    # Coulomb logarithm
    cou_log = log(M_bh / mean(dat.m))

    # Constants
    c1 = -(G ^ 2 * cou_log) / (vel ^ 2 * r ^ 2)
    c2 = (2 * G ^ 2 * cou_log) / (3 * r ^ 2)

    # Integrands
    int1(v) = DF.F(r,v)
    int2(v) = DF.N(r,v)
    int3(v) = (v ^ 2) / (vel ^ 3) * DF.F(r,v)
    int4(v) = (1 / v) * DF.F(r,v)

    # Compute integrals via midpoint
    calc1 = midpoint(x->int1(x), dat.vmin, vel, n_steps)
    calc2 = midpoint(x->int2(x), dat.vmin, vel, n_steps)
    calc3 = midpoint(x->int3(x), dat.vmin, vel, n_steps)
    calc4 = midpoint(x->int4(x), vel, dat.vmax, n_steps)

    # Calculate diff coeffs
    d_par = c1 * (calc1 + m_obj * calc2) # Eq 1.32
    d_par_sq = c2 * (calc3 + calc4) # Eq 1.33
    d_tan_sq = c2 * ((3/vel)*calc1 - calc3 + 2*calc4) # Eq 1.34

    return d_par, d_par_sq, d_tan_sq
end

function find_coef_log(r,E,L,m_obj) # Logarithmic sampling

    # Object total velocity; Convert (r,v) -> (logr,logv)
    vel = sqrt(find_vr(r,E,L) ^ 2 + find_vt(r,E,L) ^ 2)

    logr = log(r)
    log_vel = log(vel)

    # Coulomb logarithm
    cou_log = log(M_bh / mean(dat.m))

    # Constants
    c1 =(-G^2 * cou_log) / (vel ^ 2 * r ^ 3)
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

    # Calculate diffusion coefficients
    d_par = c1 * (calc1 + m_obj * calc2) # Eq 1.42
    d_par_sq = c2 * (calc3 + calc4) # Eq 1.43
    d_tan_sq = c2 * ((3 / vel) * calc1 - calc3 + 2 * calc4) # Eq 1.44

    return d_par, d_par_sq, d_tan_sq
end

## Compute local diffusion coefficients ##
function find_loc_coef(r,E,L,m_obj)

    # Object total velocity 
    vel = sqrt(find_vr(r,E,L)^2 + find_vt(r,E,L)^2)

    # Calculate velocity diffusion coefficients
    vp, vp_sq, vt_sq = find_coef_log(r,E,L,m_obj)

    # Calculate local diffusion coefficients
    delE = 0.5 * (vp_sq + vt_sq) + vel * vp # Eq 1.45
    delEsq = vel ^ 2 * vp_sq # Eq 1.46

    delL = (L / vel) * vp + (r^2)/(4*L)*vt_sq # Eq 1.47
    delLsq = (L ^ 2 / vel ^ 2) * vp_sq + 0.5*(r^2 - L^2 / vel^2) * vt_sq # Eq 1.48

    delEL = L * vp_sq # Eq 1.49

    return delE, delEsq, delL, delLsq, delEL
end

## Compute orbit averaged diffusion coefficients ##
function find_avg_coef_el(E,L,m_obj)

    # Find orbit bounds, period
    rp, ra = find_rp_ra(E,L)
    T = find_period(E,L)

    # Storage array
    result = zeros(5)

    # Simultaneous midpoint
    h = pi / n_steps
    for i in 1:n_steps
        tt = -pi/2 + (i - 0.5) * h

        r_hat = (ra - rp) / 2 # Eq 1.18
        t_hat = (ra + rp) / 2 # Eq 1.19
        rr = r_hat * sin(tt) + t_hat # Eq 1.17
        vr_tilde = find_vr(rr,E,L) / (r_hat * cos(tt))

        result .+= find_loc_coef(rr,E,L,m_obj) ./ vr_tilde
    end

    integral = h .* result

    return (2/T) .* integral # Return all coeffs
end

function find_avg_coef_ej(E,j,m_obj)

    # Compute circular angular momentum, derivatives
    Jc = find_Lc(E)[2]
    dJc, dJc2 = dLc_dE(E), d2Lc_dE2(E)
    J = j * Jc

    # Compute orbit averaged coeffs in (E,L)
    DE, DEE, DJ, DJJ, DEJ = find_avg_coef_el(E,J,m_obj)

    # Compute orbit averaged coeffs in (E,j)
    DE_j = DE
    DEE_j = DEE

    # Eq 1.63
    j1 = (DJ / Jc) - (J/(Jc^2)) * dJc * DE - (1/(Jc^2)) * dJc * DEJ
    j2 = -0.5 * ((1/(Jc^2)) * dJc2 - (2/(Jc^3)) * (dJc^2)) * DEE

    DJ_j = j1 + j2

    # Eq 1.66
    jj1 = ((J * dJc / (Jc^2))^2)*DEE
    jj2 = -2 * ((J * dJc / (Jc^3))) * DEJ + DJJ / (Jc^2)

    DJJ_j = jj1 + jj2

    # Eq 1.65
    DEJ_j = -(J * dJc /(Jc^2)) * DEE + DEJ / Jc

    return DE_j, DEE_j, DJ_j, DJJ_j, DEJ_j
end

## Compute coefficients over grid for sampling ##
function find_coef_grid(E_tab,j_tab,m_obj)

    # Initialize matrices
    DE_tab = zeros(res,res)
    DEE_tab = zeros(res,res)
    DJ_tab = zeros(res,res)
    DJJ_tab = zeros(res,res)
    DEJ_tab = zeros(res,res)

    # Loop over energy, momentum
    Threads.@threads for i in 1:res
        E = E_tab[i]
        for k in 1:res
            j = j_tab[k]
            DE, DEE, DJ, DJJ, DEJ = find_avg_coef_ej(E, j, m_obj)
            @inbounds begin
                DE_tab[k,i]  = DE
                DEE_tab[k,i] = DEE
                DJ_tab[k,i]  = DJ
                DJJ_tab[k,i] = DJJ
                DEJ_tab[k,i] = DEJ
            end
        end
    end
    return DE_tab, DEE_tab, DJ_tab, DJJ_tab, DEJ_tab
end

## Struct to store sampled coefficient grid, interpolated functions ##
struct DiffusionCoeffs

    # Orbit-averaged diff coeff grid
    DE_tab::Matrix{Float64}
    DEE_tab::Matrix{Float64}
    Dj_tab::Matrix{Float64}
    Djj_tab::Matrix{Float64}
    DEj_tab::Matrix{Float64}

    # Calculated energy and bounds (quantile)
    E_tab::AbstractVector{Float64}
    j_tab::AbstractVector{Float64}

    emin::Float64
    emax::Float64

    # Test object mass 
    m_obj::Float64

    # Interpolated functions 
    DE_samp::Function
    DEE_samp::Function
    Dj_samp::Function
    Djj_samp::Function
    DEj_samp::Function
end

## Automatically generate diffusion coeffs for snapshot, given m_obj ##
function generate_coeffs(m_obj)

    # Calculate energy for entire distribution, bounds 
    EE = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
    emin = quantile(EE,1e-5)
    emax = quantile(EE,0.99999)

    # E, j sampling range 
    E_tab = range(emin,emax,res)
    j_tab = range(1e-5,0.99999,res)

    # Calculate coefficients over grid 
    DE_tab, DEE_tab, Dj_tab, Djj_tab, DEj_tab = find_coef_grid(E_tab,j_tab,m_obj)

    # Interpolate w/ bilinear interpolation
    DE_samp  = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DE_tab)
    DEE_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DEE_tab)
    Dj_samp  = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, Dj_tab)
    Djj_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, Djj_tab)
    DEj_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DEj_tab)

    # Return struct 
    return DiffusionCoeffs(
        DE_tab, DEE_tab, Dj_tab, Djj_tab, DEj_tab,
        E_tab, j_tab, emin, emax, m_obj,
        DE_samp, DEE_samp, Dj_samp, Djj_samp, DEj_samp
    )
end

## Save orbit-averaged diff coeff struct to HDF5 file ##
function save_coeffs(name::String,dc::DiffusionCoeffs)
    h5open(name, "w") do file
        file["E_tab"] = collect(dc.E_tab)
        file["j_tab"] = collect(dc.j_tab)
        file["emin"] = dc.emin
        file["emax"] = dc.emax
        file["DE_tab"] = dc.DE_tab
        file["DEE_tab"] = dc.DEE_tab
        file["Dj_tab"] = dc.Dj_tab
        file["Djj_tab"] = dc.Djj_tab
        file["DEj_tab"] = dc.DEj_tab
    end
end

## Open diff coeff HDF5 file ##
function load_coeffs(name::String)
    h5open(name, "r") do file

        # Read in coeff dat
        E_tab = read(file["E_tab"])
        j_tab = read(file["j_tab"])
        emin = read(file["emin"])
        emax = read(file["emax"])
        DE_tab = read(file["DE_tab"])
        DEE_tab = read(file["DEE_tab"])
        Dj_tab = read(file["Dj_tab"])
        Djj_tab = read(file["Djj_tab"])
        DEj_tab = read(file["DEj_tab"])

        # Interpolate coeffs for MC
        DE_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DE_tab)
        DEE_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DEE_tab)
        Dj_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, Dj_tab)
        Djj_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, Djj_tab)
        DEj_samp = (j, E) -> bilinear_interp(j, E, j_tab, E_tab, DEj_tab)
        return DiffusionCoeffs(
            DE_tab, DEE_tab, Dj_tab, Djj_tab, DEj_tab,
            E_tab, j_tab, emin, emax, m_obj,
            DE_samp, DEE_samp, Dj_samp, Djj_samp, DEj_samp
        )
    end
end
# --- Circular angular momentum, derivatives for given E --- #

## Calculate rc, Lc via effective potential ##
function find_Lc(E)

    N = length(psi_rtab)

    # psi(r) + (L^2)/(2*r^2) = E -> L^2 = eta = 2*r^2*(E-psi(r))

    # Value at leftmost two points
    eta_1 = 2*psi_rtab[1]^2 * (E - psi_tot_tab[1])
    eta_2 = 2*psi_rtab[2]^2 * (E - psi_tot_tab[2])

    # Value at rightmost two points
    eta_N1 = 2*psi_rtab[N-1]^2 * (E - psi_tot_tab[N-1])
    eta_N = 2*psi_rtab[N]^2 * (E - psi_tot_tab[N])

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

            eta_m = 2*psi_rtab[im]^2 * (E - psi_tot_tab[im])
            eta_m1 = 2*psi_rtab[im+1]^2 * (E - psi_tot_tab[im+1])

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
        rl = psi_rtab[il]
        rl1 = psi_rtab[il+1]

        psil = psi_tot_tab[il]
        psil1 = psi_tot_tab[il+1]

        slopel = (psil1-psil)/(1/rl1-1/rl)

        rsol_l = 0.5 * 1/((E-psil)/slopel + 1/rl)

        # Linear interpolation (right bracket)
        rr = psi_rtab[ir]
        rr1 = psi_rtab[ir+1]

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
            r1 = psi_rtab[1]
            r2 = psi_rtab[2]

            psi1 = psi_tot_tab[1]
            psi2 = psi_tot_tab[2]

            slope_12 = (psi2 - psi1) / (1/r2 - 1/r1)

            rsol = 0.5 * 1 / ((E-psi1)/slope_12 + 1/r1)
            psi = psi1 + slope_12*(1/rsol - 1/r1)
        end

        Lc2 = max(0.0,2*rsol^2*(E-psi)) # Protect against nonphysical solns

    else # if del_eta_r > 0

        if eta_N > eta_N1 # rc > rN -> Analytic soln
            Mtot = psi_Mtot + M_bh # total mass
            rsol = -(G*Mtot)/(2*E)
            psi = -(G*Mtot)/rsol

        else # Maximum in [rN-1,rN] -> Linear interpolation
            rN1 = psi_rtab[N-1]
            rN = psi_rtab[N]

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

## Lc(E) derivatives for the (E,L) -> (E,j) coefficient transform (Dj1 etc.) ##
#
# Two backends, selected by the LC_DERIV_MODE flag (default :fd = the original behavior):
#   :fd     -- 3-point finite differences of find_Lc (ORIGINAL; production default).
#   :spline -- analytic derivatives of a natural cubic spline fit to Lc(E) (see the spline
#              block at the bottom of this file). Motivation: d2Lc_dE2 by finite difference
#              amplifies the interpolation kinks in find_Lc by ~1/eps^2, which is a prime
#              suspect for the jagged sign-change / "hook" feature in the Dj1 map. The
#              spline smooths Lc(E) first, so its 2nd derivative is clean. Diagnostic /
#              opt-in only -- flip via set_Lc_deriv_mode!(:spline) after build_Lc_spline!.
#
# Finite-difference step is relative: an absolute one is noise-dominated for deeply bound E
# and can push E across 0 (unbound) for weakly bound E.
function dLc_dE_fd(E)
    eps_E = 1e-4 * abs(E)
    return (find_Lc(E+eps_E)[2]-find_Lc(E-eps_E)[2])/(2*eps_E)
end

function d2Lc_dE2_fd(E)
    eps_E = 1e-4 * abs(E)
    return (find_Lc(E+eps_E)[2]+find_Lc(E-eps_E)[2]-2*find_Lc(E)[2])/(eps_E^2)
end

# Spline backend: analytic 1st/2nd derivatives of the cached Lc(E) spline (LC_SPLINE).
function dLc_dE_spline(E)
    sp = LC_SPLINE[]
    sp === nothing && error("dLc_dE_spline: LC_SPLINE not built; call build_Lc_spline! first.")
    return spline_d1(sp, E)
end
function d2Lc_dE2_spline(E)
    sp = LC_SPLINE[]
    sp === nothing && error("d2Lc_dE2_spline: LC_SPLINE not built; call build_Lc_spline! first.")
    return spline_d2(sp, E)
end

# Flag-dispatched public API (call sites: avg_coef_ej, mc_drift). Default :fd reproduces
# the original finite-difference behavior exactly.
dLc_dE(E)   = LC_DERIV_MODE === :spline ? dLc_dE_spline(E)   : dLc_dE_fd(E)
d2Lc_dE2(E) = LC_DERIV_MODE === :spline ? d2Lc_dE2_spline(E) : d2Lc_dE2_fd(E)


# ======================================================================================
# Natural cubic spline of Lc(E) (diagnostic / opt-in backend for the derivatives above).
# Pure-Julia, no external deps. Numerical-Recipes tridiagonal solve; value/1st/2nd deriv.
# ======================================================================================
struct CubicSpline
    x::Vector{Float64}    # knots, strictly ascending
    y::Vector{Float64}    # values
    y2::Vector{Float64}   # second derivatives at knots (natural BC: y2[1]=y2[end]=0)
end

function CubicSpline(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    (n == length(y) && n >= 3) || error("CubicSpline: need matching lengths and >= 3 knots.")
    xx = collect(float.(x)); yy = collect(float.(y))
    u  = zeros(n); y2 = zeros(n)          # natural BC -> y2[1] = u[1] = 0
    @inbounds for i in 2:n-1
        sig = (xx[i]-xx[i-1]) / (xx[i+1]-xx[i-1])
        p   = sig*y2[i-1] + 2.0
        y2[i] = (sig-1.0)/p
        ui  = (yy[i+1]-yy[i])/(xx[i+1]-xx[i]) - (yy[i]-yy[i-1])/(xx[i]-xx[i-1])
        u[i] = (6.0*ui/(xx[i+1]-xx[i-1]) - sig*u[i-1]) / p
    end
    y2[n] = 0.0
    @inbounds for k in n-1:-1:1
        y2[k] = y2[k]*y2[k+1] + u[k]
    end
    return CubicSpline(xx, yy, y2)
end

@inline function _spline_klo(sp::CubicSpline, t::Real)
    t <= sp.x[1]   && return 1
    t >= sp.x[end] && return length(sp.x) - 1     # clamp; extrapolates on the end cubic
    return searchsortedlast(sp.x, t)              # x[klo] <= t <= x[klo+1]
end

function spline_val(sp::CubicSpline, t::Real)     # value
    klo = _spline_klo(sp, t); khi = klo+1
    h = sp.x[khi]-sp.x[klo]; a = (sp.x[khi]-t)/h; b = (t-sp.x[klo])/h
    return a*sp.y[klo] + b*sp.y[khi] + ((a^3-a)*sp.y2[klo] + (b^3-b)*sp.y2[khi])*h^2/6
end

function spline_d1(sp::CubicSpline, t::Real)      # first derivative
    klo = _spline_klo(sp, t); khi = klo+1
    h = sp.x[khi]-sp.x[klo]; a = (sp.x[khi]-t)/h; b = (t-sp.x[klo])/h
    return (sp.y[khi]-sp.y[klo])/h + (-(3a^2-1)*sp.y2[klo] + (3b^2-1)*sp.y2[khi])*h/6
end

function spline_d2(sp::CubicSpline, t::Real)      # second derivative (linear in each panel)
    klo = _spline_klo(sp, t); khi = klo+1
    h = sp.x[khi]-sp.x[klo]; a = (sp.x[khi]-t)/h; b = (t-sp.x[klo])/h
    return a*sp.y2[klo] + b*sp.y2[khi]
end

# Cached Lc(E) spline + the derivative-mode flag. LC_DERIV_MODE is a typed global (read in
# dLc_dE, i.e. avg_coef_ej / mc_drift -- not the MC hot loop, but typed for consistency).
const LC_SPLINE = Ref{Union{Nothing,CubicSpline}}(nothing)
LC_DERIV_MODE::Symbol = :fd     # :fd (finite diff, PRODUCTION DEFAULT) or :spline

# Sample Lc(E) on a dense log-|E| grid over [E_lo, E_hi] (both < 0) and spline it, so the
# spline backend has smooth analytic dLc/dE, d2Lc/dE2. Call before set_Lc_deriv_mode!(:spline).
function build_Lc_spline!(E_lo::Real, E_hi::Real; n::Integer = 4000)
    (E_lo < 0 && E_hi < 0) || error("build_Lc_spline!: E_lo, E_hi must be < 0.")
    aE = 10 .^ range(log10(min(abs(E_lo), abs(E_hi))), log10(max(abs(E_lo), abs(E_hi))), length = n)
    Egrid  = sort(-aE)                                  # ascending E
    Lcgrid = [find_Lc(E)[2] for E in Egrid]
    LC_SPLINE[] = CubicSpline(Egrid, Lcgrid)
    println(">>> build_Lc_spline!: $n knots, E in [$(Egrid[1]), $(Egrid[end])], " *
            "Lc in [$(minimum(Lcgrid)), $(maximum(Lcgrid))]")
    flush(stdout)
    return nothing
end

function set_Lc_deriv_mode!(mode::Symbol)
    mode in (:fd, :spline) || error("set_Lc_deriv_mode!: mode must be :fd or :spline.")
    global LC_DERIV_MODE = mode
    println(">>> LC_DERIV_MODE = $mode   (Lc(E) derivatives used in the (E,L)->(E,j) transform)")
    flush(stdout)
    return nothing
end



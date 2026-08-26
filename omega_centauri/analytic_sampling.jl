# ============================================================================
# Analytic single-mass, isotropic power-law cusp around an IMBH -- synthetic
# "CMC-like" snapshot generator for the EPIC analytic / convergence tests.
#
# Model (Phase-0 "case B", self-consistent Keplerian cusp):
#   * density   rho(r) = A r^(-gamma)  over [r_min, r_infl]   (gamma is an argument)
#   * potential pure Keplerian, psi(r) = G M_bh / r  (stellar self-gravity is small
#     deep in the cusp; ~factor-2 near r_infl -- the closed-form DF is exact for
#     r << r_infl and approximate at the outer edge)
#   * DF        isotropic Keplerian closed form  f(eps) ~ eps^(gamma - 3/2)
#   * outer edge r_max == r_infl, so M_star(<r_infl) = M_bh and N = M_bh/m_star
#   * code units: CMC/Henon style (G = M_0 = 1), pinned to r_infl; the physical
#     size is chosen so c(code) matches the realistic COLOSSUS runs (~2632).
#
# Built up in pieces. PIECE 1 (this file so far): the model definition, the
# derived quantities + CMC-style code units, and parameter logging.
# NEXT: radius sampler -> velocity sampler -> snapshot (.h5) + .conv.sh writers.
# ============================================================================

using Printf
using Random
using HDF5

# --- Physical constants (astro units). NB: deliberately NOT named `G` -- the
#     pipeline defines `const G = 1`, so we avoid a clash if this file is ever
#     loaded alongside main.jl. ---
const G_ASTRO   = 4.30091e-3             # G [pc (km/s)^2 / Msun]
const C_KMS     = 2.99792458e5           # speed of light [km/s]
const PC_TO_CM  = 3.0856775814913673e18  # [cm / pc]
const PC_TO_KM  = 3.0856775814913673e13  # [km / pc]
const MSUN_TO_G = 1.988409870698051e33   # [g / Msun]
const MYR_TO_S  = 3.1556952e13           # [s / Myr]

struct AnalyticCusp
    # --- physical inputs ---
    gamma      :: Float64    # density slope, rho ~ r^-gamma
    M_bh       :: Float64    # IMBH mass [Msun]
    m_star     :: Float64    # single particle mass [Msun]
    r_min      :: Float64    # inner sampling radius [pc]
    r_infl     :: Float64    # outer radius = influence radius [pc]  (== r_max)
    gamma_coul :: Float64    # Coulomb-log arg in ln(gamma_coul*N) (Eq 52; NOT the slope)
    # --- derived model ---
    N          :: Int        # particle count = round(M_bh / m_star)
    M_star     :: Float64    # total stellar mass = N*m_star [Msun] (case B ~ M_bh)
    rho_coeff  :: Float64    # A in rho(r) = A r^-gamma  [Msun/pc^3]
    # --- CMC-style code units (Henon, pinned to r_infl) ---
    massunitmsun   :: Float64
    lengthunitpc   :: Float64
    v_unit_kms     :: Float64
    c_code         :: Float64
    nbtimeunitsmyr :: Float64   # dynamical (crossing) time unit [Myr]
    timeunitsmyr   :: Float64   # relaxation time unit [Myr] (Eq 52)
    massunitcgs    :: Float64
    lengthunitcgs  :: Float64
    nbtimeunitcgs  :: Float64
    # --- diagnostics ---
    r_g            :: Float64    # IMBH gravitational radius G M_bh/c^2 [pc]
    sigma_infl_kms :: Float64    # ~ velocity dispersion near r_infl
    t_relax_myr    :: Float64    # ~ local relaxation time near r_infl (order of magnitude)
end

"""
    AnalyticCusp(; gamma=1.75, M_bh=1e5, m_star=30.0, r_min=1e-5, r_infl=0.033,
                   gamma_coul=0.11)

Build the case-B cusp model: derive N, the density amplitude, the CMC-style code
units (pinned to r_infl), and a few diagnostics. Defaults are the locked Phase-0
model (c(code) ~ 2632). `gamma`, `M_bh`, `m_star`, `r_min`, `r_infl` are the only
independent inputs; everything else follows.
"""
function AnalyticCusp(; gamma::Real = 1.75, M_bh::Real = 1.0e5, m_star::Real = 30.0,
                        r_min::Real = 1.0e-5, r_infl::Real = 0.033,
                        gamma_coul::Real = 0.11)

    gamma < 3.0 || error("AnalyticCusp: need gamma < 3 for finite enclosed mass (got $gamma).")
    0 < r_min < r_infl ||
        error("AnalyticCusp: need 0 < r_min < r_infl (got r_min=$r_min, r_infl=$r_infl).")

    # --- case-B model: total stellar mass = M_bh -> N particles of m_star ---
    N      = round(Int, M_bh / m_star)
    M_star = N * m_star

    # density amplitude A from  M_star = ∫_{r_min}^{r_infl} 4π A r^(2-γ) dr
    p         = 3.0 - gamma
    rho_coeff = M_star * p / (4π * (r_infl^p - r_min^p))

    # --- CMC-style code units (Henon), pinned to r_infl ---
    massunitmsun   = M_star
    lengthunitpc   = r_infl
    v_unit_kms     = sqrt(G_ASTRO * M_star / r_infl)                  # = U_l/U_t, code velocity unit
    c_code         = C_KMS / v_unit_kms
    nbtimeunitsmyr = (r_infl * PC_TO_KM) / (v_unit_kms * MYR_TO_S)    # dynamical time unit
    timeunitsmyr   = (N / log(gamma_coul * N)) * nbtimeunitsmyr       # relaxation unit (Eq 52)
    massunitcgs    = M_star * MSUN_TO_G
    lengthunitcgs  = r_infl * PC_TO_CM
    nbtimeunitcgs  = nbtimeunitsmyr * MYR_TO_S

    # --- diagnostics ---
    r_g      = G_ASTRO * M_bh / C_KMS^2                               # gravitational radius [pc]
    sigma    = sqrt(G_ASTRO * M_bh / ((1.0 + gamma) * r_infl))        # cusp sigma near r_infl [km/s]
    lnLambda = log(M_bh / m_star)                                     # matches the coefficient Coulomb log
    rho_infl = rho_coeff * r_infl^(-gamma)
    # Spitzer local relaxation time near r_infl (order-of-magnitude), pc/(km/s) -> Myr
    t_relax  = 0.34 * sigma^3 / (G_ASTRO^2 * rho_infl * m_star * lnLambda) *
               (PC_TO_KM / MYR_TO_S)

    return AnalyticCusp(gamma, M_bh, m_star, r_min, r_infl, gamma_coul,
                        N, M_star, rho_coeff,
                        massunitmsun, lengthunitpc, v_unit_kms, c_code,
                        nbtimeunitsmyr, timeunitsmyr,
                        massunitcgs, lengthunitcgs, nbtimeunitcgs,
                        r_g, sigma, t_relax)
end

"""
    log_params(m::AnalyticCusp)

Print every input, derived quantity, and code unit to stdout (for the job log).
"""
function log_params(m::AnalyticCusp)
    println("="^74)
    println(" Analytic power-law cusp -- model parameters")
    println("="^74)
    @printf("  gamma (density slope)         = %.4f\n",        m.gamma)
    @printf("  M_bh                          = %.6g Msun\n",   m.M_bh)
    @printf("  m_star (particle mass)        = %.6g Msun\n",   m.m_star)
    @printf("  N (= round(M_bh/m_star))      = %d\n",          m.N)
    @printf("  M_star (= N*m_star)           = %.6g Msun\n",   m.M_star)
    @printf("  r_min                         = %.6g pc\n",     m.r_min)
    @printf("  r_infl (= r_max)              = %.6g pc\n",     m.r_infl)
    @printf("  rho amplitude A (rho=A r^-g)  = %.6g Msun/pc^3\n", m.rho_coeff)
    println("-"^74)
    println(" CMC-style code units (Henon, pinned to r_infl):")
    @printf("  massunitmsun                  = %.6g Msun\n",   m.massunitmsun)
    @printf("  lengthunitparsec              = %.6g pc\n",     m.lengthunitpc)
    @printf("  v_unit                        = %.6g km/s\n",   m.v_unit_kms)
    @printf("  c (code)                      = %.6g\n",        m.c_code)
    @printf("  nbtimeunitsmyr (dynamical)    = %.6g Myr\n",    m.nbtimeunitsmyr)
    @printf("  timeunitsmyr (relaxation)     = %.6g Myr   [gamma_coul=%.3g]\n",
            m.timeunitsmyr, m.gamma_coul)
    @printf("  massunitcgs                   = %.6g g\n",      m.massunitcgs)
    @printf("  lengthunitcgs                 = %.6g cm\n",     m.lengthunitcgs)
    @printf("  nbtimeunitcgs                 = %.6g s\n",      m.nbtimeunitcgs)
    @printf("  M_bh (code units)             = %.6g\n",        m.M_bh / m.massunitmsun)
    println("-"^74)
    println(" Diagnostics:")
    @printf("  r_g = G M_bh/c^2              = %.6g pc   (r_infl/r_g = %.6g)\n",
            m.r_g, m.r_infl / m.r_g)
    @printf("  sigma near r_infl            ~ %.6g km/s\n",    m.sigma_infl_kms)
    @printf("  t_relax near r_infl          ~ %.6g Myr   (order-of-magnitude)\n", m.t_relax_myr)
    println("="^74)
    flush(stdout)
    return nothing
end


# ============================================================================
# PIECE 2: radius sampler
# ============================================================================

"""
    density(m::AnalyticCusp, r) -> Msun/pc^3

Analytic mass density rho(r) = A r^(-gamma); 0 outside [r_min, r_infl].
"""
function density(m::AnalyticCusp, r::Real)
    (m.r_min <= r <= m.r_infl) || return 0.0
    return m.rho_coeff * r^(-m.gamma)
end

"""
    radius_cdf(m::AnalyticCusp, r) -> in [0,1]

Analytic enclosed-mass fraction P(<r) = M_star(<r)/M_star over [r_min, r_infl]:
`P(<r) = (r^p - r_min^p)/(r_infl^p - r_min^p)`, p = 3 - gamma.
"""
function radius_cdf(m::AnalyticCusp, r::Real)
    p = 3.0 - m.gamma
    r <= m.r_min  && return 0.0
    r >= m.r_infl && return 1.0
    return (r^p - m.r_min^p) / (m.r_infl^p - m.r_min^p)
end

"""
    sample_radii(m::AnalyticCusp; n=m.N, rng=default_rng()) -> Vector{Float64} [pc]

Draw `n` radii from the cusp by exact inverse-CDF (no rejection). Inverting
`P(<r)` gives `r = ( r_min^p + u (r_infl^p - r_min^p) )^(1/p)`, `u ~ U(0,1)`.
Because the mass is weighted toward large r (`dN/dr ∝ r^(2-gamma)`), most radii
land near r_infl and the innermost sampled radius is ~ r_infl * N^(-1/(3-gamma)),
well above r_min for these parameters.
"""
function sample_radii(m::AnalyticCusp; n::Integer = m.N,
                      rng::AbstractRNG = Random.default_rng())
    p    = 3.0 - m.gamma
    a    = m.r_min^p
    span = m.r_infl^p - a
    invp = 1.0 / p
    radii = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        radii[i] = (a + rand(rng) * span)^invp
    end
    return radii
end

"""
    validate_radii(m::AnalyticCusp, radii; nbins=30)

Quick check that the sampled radii reproduce `rho ∝ r^(-gamma)`: fit the density
slope in log-log over log-spaced shells (populated bins only) and report the
recovered gamma, the populated radial range, and the KS max deviation between the
empirical and analytic CDFs. Returns `(gamma_fit, ks, r_lo, r_hi)`.
"""
function validate_radii(m::AnalyticCusp, radii::AbstractVector{<:Real}; nbins::Integer = 30)
    n = length(radii)
    edges = 10.0 .^ range(log10(m.r_min), log10(m.r_infl), length = nbins + 1)

    logr = Float64[]; logrho = Float64[]; wts = Float64[]
    for b in 1:nbins
        r_a, r_b = edges[b], edges[b + 1]
        cnt = count(x -> r_a <= x < r_b, radii)
        cnt == 0 && continue
        vol = (4π / 3) * (r_b^3 - r_a^3)           # shell volume [pc^3]
        push!(logr,   log10(sqrt(r_a * r_b)))       # geometric bin center
        push!(logrho, log10(m.m_star * cnt / vol))  # rho = m_star * n
        push!(wts,    float(cnt))                    # Poisson weight ~ count (down-weight sparse inner bins)
    end

    # count-weighted least-squares slope of log(rho) vs log(r)  (expect -gamma)
    nb   = length(logr)
    W    = sum(wts)
    xbar = sum(wts .* logr) / W;  ybar = sum(wts .* logrho) / W
    sxx  = sum(wts .* (logr .- xbar) .^ 2)
    sxy  = sum(wts .* (logr .- xbar) .* (logrho .- ybar))
    gamma_fit = -(sxy / sxx)

    # KS-style max |empirical CDF - analytic CDF|
    sorted = sort(radii)
    ksmax = 0.0
    @inbounds for (i, r) in enumerate(sorted)
        ksmax = max(ksmax, abs(i / n - radius_cdf(m, r)))
    end

    @printf("  radius check: N=%d, fitted density slope = %.4f (input gamma = %.4f), KS max|dCDF| = %.4f\n",
            n, gamma_fit, m.gamma, ksmax)
    @printf("               populated bins span r = %.3g .. %.3g pc (%d/%d bins), min sampled r = %.3g pc\n",
            10.0^minimum(logr), 10.0^maximum(logr), nb, nbins, minimum(radii))
    flush(stdout)
    return (gamma_fit = gamma_fit, ks = ksmax, r_lo = minimum(radii), r_hi = maximum(radii))
end


# ============================================================================
# PIECE 3: velocity sampler
# ============================================================================

## Keplerian relative potential psi(r) = G M_bh / r  [(km/s)^2], r in pc. ##
psi_kepler(m::AnalyticCusp, r::Real) = G_ASTRO * m.M_bh / r

## Isotropic 1D dispersion of a Keplerian cusp: sigma_1D^2 = G M_bh/((1+gamma) r). ##
sigma_1d(m::AnalyticCusp, r::Real) = sqrt(G_ASTRO * m.M_bh / ((1.0 + m.gamma) * r))

"""
    sample_speed(m, r, rng) -> v [km/s]

Draw a speed at radius r from the local distribution p(v|r) ∝ v^2 f(eps), with the
Keplerian closed-form DF f(eps) ∝ eps^(gamma-3/2), eps = psi(r) - v^2/2, and
v in [0, v_esc = sqrt(2 psi)]. Rejection sampling against the analytic mode
(v_peak^2 = 2 psi/(1+q), q = gamma-3/2). Assumes q >= 0 (gamma >= 3/2); for a
shallower cusp p(v) diverges at v_esc and would need a different envelope.
"""
function sample_speed(m::AnalyticCusp, r::Real, rng::AbstractRNG)
    q    = m.gamma - 1.5
    psi  = psi_kepler(m, r)
    vesc = sqrt(2 * psi)
    vp2  = 2 * psi / (1 + q)                    # mode: v_peak^2
    pmax = vp2 * max(psi - 0.5 * vp2, 0.0)^q    # envelope height p(v_peak)
    while true
        v   = vesc * rand(rng)
        eps = max(psi - 0.5 * v^2, 0.0)
        p   = v^2 * eps^q
        rand(rng) * pmax <= p && return v
    end
end

"""
    sample_velocities(m, radii; rng) -> (vr, vt)   [km/s each]

For every radius draw a speed from `sample_speed`, then an isotropic direction:
`mu = cos(angle) ~ U(-1,1)`, `vr = v*mu` (signed), `vt = v*sqrt(1-mu^2)` (>= 0).
Matches the SnapDat convention (vr signed, vt = tangential magnitude, L = r*vt).
"""
function sample_velocities(m::AnalyticCusp, radii::AbstractVector{<:Real};
                           rng::AbstractRNG = Random.default_rng())
    n  = length(radii)
    vr = Vector{Float64}(undef, n)
    vt = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        v     = sample_speed(m, radii[i], rng)
        mu    = 2 * rand(rng) - 1
        vr[i] = v * mu
        vt[i] = v * sqrt(max(1 - mu^2, 0.0))
    end
    return vr, vt
end

"""
    validate_velocities(m, radii, vr, vt; nbins=20, mincount=20)

Check the sampled velocities against the Keplerian cusp: in log-r bins compare the
measured radial dispersion to sigma_1D(r) = sqrt(G M_bh/((1+gamma) r)) and the
tangential/radial ratio (isotropy => ~1). Reports mean +/- spread of
sigma_r(measured)/sigma_1D(predicted) and of sigma_t/sigma_r over populated bins.
"""
function validate_velocities(m::AnalyticCusp, radii::AbstractVector{<:Real},
                             vr::AbstractVector{<:Real}, vt::AbstractVector{<:Real};
                             nbins::Integer = 20, mincount::Integer = 20)
    edges = 10.0 .^ range(log10(m.r_min), log10(m.r_infl), length = nbins + 1)
    ratio = Float64[]; iso = Float64[]
    for b in 1:nbins
        r_a, r_b = edges[b], edges[b + 1]
        idx = findall(x -> r_a <= x < r_b, radii)
        length(idx) >= mincount || continue
        r_c = sqrt(r_a * r_b)
        sr  = sqrt(sum(vr[i]^2 for i in idx) / length(idx))          # <vr^2>^(1/2) = sigma_1D
        st1 = sqrt(sum(vt[i]^2 for i in idx) / length(idx) / 2)      # <vt^2> = 2 sigma_1D^2
        push!(ratio, sr / sigma_1d(m, r_c))
        push!(iso,   st1 / sr)
    end
    nb     = length(ratio)
    mean_r = sum(ratio) / nb
    sd_r   = sqrt(sum((ratio .- mean_r) .^ 2) / nb)
    mean_i = sum(iso) / nb
    @printf("  velocity check: sigma_r/sigma_pred = %.4f +/- %.4f over %d bins (expect 1.0), sigma_t/sigma_r = %.4f (expect 1.0)\n",
            mean_r, sd_r, nb, mean_i)
    flush(stdout)
    return (sigma_ratio = mean_r, sigma_ratio_sd = sd_r, iso_ratio = mean_i, nbins = nb)
end


# ============================================================================
# PIECE 4: CMC-style snapshot (.h5) + .conv.sh writers
# ============================================================================

# One compound-table row, matching the fields read_file_h5 requires. r/vr/vt are
# in CODE units; m_MSUN is in solar masses (read_file_h5 divides it by massunitmsun).
struct SnapRow
    m_MSUN   :: Float64
    r        :: Float64
    vr       :: Float64
    vt       :: Float64
    startype :: Float64
    binflag  :: Float64
end

"""
    write_snapshot(m, radii, vr, vt, path; time_gyr=1.0, snap_index=1,
                   startype=14, binflag=0) -> key

Write a CMC-style compound snapshot to `path`. Physical `(radii [pc], vr, vt [km/s])`
are converted to code units via the model's unit fields; `m_MSUN = m_star` (solar);
all particles are `startype` (14 = BH), `binflag` (0 = single). Particles are sorted
by ascending radius (CMC convention; find_psi_arrays integrates shells that way). The
dataset key is `"<snap_index>(t=<time_gyr>Gyr)"`, which snapshot_keys_by_time parses.
The IMBH mass is deliberately NOT stored -- encode it in the filename and pass it to
EPIC at run time. Returns the dataset key.
"""
function write_snapshot(m::AnalyticCusp, radii::AbstractVector{<:Real},
                        vr::AbstractVector{<:Real}, vt::AbstractVector{<:Real},
                        path::AbstractString; time_gyr::Real = 1.0, snap_index::Integer = 1,
                        startype::Real = 14.0, binflag::Real = 0.0)
    n = length(radii)
    (length(vr) == n && length(vt) == n) || error("write_snapshot: radii/vr/vt length mismatch.")

    perm = sortperm(radii)             # ascending r (CMC / find_psi_arrays convention)
    rc   = m.lengthunitpc              # pc per code length
    vc   = m.v_unit_kms                # km/s per code velocity

    rows = Vector{SnapRow}(undef, n)
    @inbounds for (k, i) in enumerate(perm)
        rows[k] = SnapRow(m.m_star, radii[i] / rc, vr[i] / vc, vt[i] / vc,
                          float(startype), float(binflag))
    end

    key = "$(snap_index)(t=$(time_gyr)Gyr)"
    h5open(path, "w") do f
        f[key] = rows
    end
    @printf("wrote snapshot %s  [%d particles, key = \"%s\", startype = %d]\n",
            path, n, key, Int(startype))
    flush(stdout)
    return key
end

"""
    write_conv_sh(m, path)

Emit a CMC-style `.conv.sh` (bash `key=value` lines) with exactly the fields
set_model_constants! parses (massunitmsun, lengthunitparsec, lengthunitcgs,
nbtimeunitcgs, nbtimeunitsmyr, timeunitsmyr), plus massunitcgs. Reading it back
reproduces c(code), r_conv, etc.
"""
function write_conv_sh(m::AnalyticCusp, path::AbstractString)
    open(path, "w") do io
        println(io, "# Synthetic .conv.sh for the analytic power-law cusp (EPIC test).")
        println(io, "# CMC-style code units (G = M_0 = 1), pinned to r_infl.")
        @printf(io, "# gamma=%.4f M_bh=%.6g Msun m_star=%.6g Msun N=%d r_infl=%.6g pc c_code=%.6g\n",
                m.gamma, m.M_bh, m.m_star, m.N, m.r_infl, m.c_code)
        @printf(io, "massunitcgs=%.10g\n",      m.massunitcgs)
        @printf(io, "massunitmsun=%.10g\n",     m.massunitmsun)
        @printf(io, "lengthunitcgs=%.10g\n",    m.lengthunitcgs)
        @printf(io, "lengthunitparsec=%.10g\n", m.lengthunitpc)
        @printf(io, "nbtimeunitcgs=%.10g\n",    m.nbtimeunitcgs)
        @printf(io, "nbtimeunitsmyr=%.10g\n",   m.nbtimeunitsmyr)
        @printf(io, "timeunitsmyr=%.10g\n",     m.timeunitsmyr)
    end
    @printf("wrote conv.sh  %s  [massunitmsun=%.6g, lengthunitparsec=%.6g pc, c(code)=%.6g]\n",
            path, m.massunitmsun, m.lengthunitpc, m.c_code)
    flush(stdout)
    return nothing
end

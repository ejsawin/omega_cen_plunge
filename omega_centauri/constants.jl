# --- Constants --- #

## Physical constants ##
const G = 1
const M = 1



# --- Type-only declarations (NO value) -- KEEP THESE, do not delete --- #
# The value lines above are commented out (the real values come from the conv.sh /
# centmass at runtime), but the ::Float64 TYPE must still be declared or these become bare
# `Any` globals -- and an untyped global read in a hot loop measured ~40x slower than a
# typed one. set_model_constants!(conv.sh) fills c/t_conv/r_conv/timeunitsmyr/massunitmsun;
# get_MBH_mass fills M_bh per snapshot. `global X::T` (no value) declares the type while
# keeping them model-agnostic. Reading one before it is assigned throws UndefVarError (the
# run scripts always assign first), which is the intended fail-fast.
global M_bh::Float64
global c::Float64
global t_conv::Float64
global r_conv::Float64
global timeunitsmyr::Float64
global massunitmsun::Float64

# Active snapshot's potential tables, set by activate_orbit_psi!() (and the DC pipeline).
# find_Lc / find_rp_ra / find_period read these every MC step, so they MUST be typed too.
global psi_rtab::Vector{Float64}
global psi_tab::Vector{Float64}
global psi_tot_tab::Vector{Float64}
global M_tab::Vector{Float64}
global psi_Mtot::Float64

const b = 3*pi/16 # Normalized to Virial radius

# Model-dependent constants. These are TYPED globals (mutable but type-stable): the
# values below are the Elena defaults, but set_model_constants!(conv_path) in data_io.jl
# overrides c/t_conv/r_conv/timeunitsmyr from a CMC .conv.sh file (once per model), and
# automate_dc/mc update M_bh per snapshot from centmass.dat. The ::Float64 annotation is
# essential: a bare global here is read as `Any` in the hot loops (gw_rates, loss_cone,
# psi) and causes a ~7x slowdown from boxing / dynamic dispatch.
#M_bh::Float64         = 1.6935299e-03 # Elena 1Gyr; per-snapshot from centmass.dat
#c::Float64            = 4349.0        # Elena; code speed of light = c_phys / v_unit
#t_conv::Float64       = 0.0709232     # Elena; N-body time unit [Myr] (the MC/orbit clock)
#r_conv::Float64       = 5.0           # Elena; code length unit [pc]
#timeunitsmyr::Float64 = 67207.0       # Elena; code time unit [Myr] (converts centmass.dat times)
#massunitmsun::Float64 = 5.52429e6     # Elena; code unit of mass [M_sun] (converts particle m_MSUN -> code)


# IMBH0.1 (for reference; use set_model_constants! + centmass instead of uncommenting --
# note set_model_constants! does NOT set M_bh, so an uncommented M_bh here silently
# overrides the per-model value and corrupts loss_cone / gw_rates / find_Lc).
#c            = 2632.7836
#M_bh         = 1.062641901e-01 # IMBH0.1, 1Gyr
#t_conv       = 0.017174
#r_conv       = 2.0



## Numerical constants ##
const n_steps = 2872 # Integration sampling points, does not impact new integrals
const tol = 1e-12 # Root finding
const lc_tol = 1e-8 # Lc buffer for rp, ra finding 
const maxiter = 100 # Max iterations (bisection)
const epsilon = 1e-2 # Discrete derivative step

const dr=0.001 # Grid size (DF Sampling- Lin (r,v))
const dv=0.001

const dlogr=0.005 # Grid size (DF Sampling- Log (r,v))
const dlogv=0.005

const dE=0.005 # Grid size (DF Sampling - Lin (E,j))
const dj=0.005

const dlogE=0.005 # Grid size (DF Sampling - Log (E,j))
const dlogj=0.005

const psi_nbins = 1000000 # Log-r bins for orbit-sampled potential
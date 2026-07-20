# --- Constants --- #

## Physical constants ##
const G = 1
const M = 1

#const M_bh = 1.6935299e-03 #Elena 1Gyr
#const c = 4349#Elena
#const t_conv = 0.0709232 # Elena
#const r_conv = 5 # Elena

const b = 3*pi/16 # Normalized to Virial radius


const c = 2632.7836 # IMBH0.1
const M_bh = 1.062641901e-01 #IMBH0.1, 1Gyr 
const t_conv = 0.017174 # IMBH0.1
const r_conv = 2 # IMBH0.1



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
# --- Constants --- #

## Physical constants ##

const G = 1
const M = 1
const M_bh = 8.500250825e-03
const c = 4542.89
const b = 3*pi/16 # Normalized to Virial radius 

## Numerical constants ##
const n_steps = 200 # Integration sampling points
const tol = 1e-20 # Root finding
const lc_tol = 1e-8 # Lc buffer for rp, ra finding 
const maxiter = 200 # Max iterations (bisection)
const epsilon = 1e-3 # Discrete derivative step

const dr=0.001 # Grid size (DF Sampling- Lin)
const dv=0.001

const dlogr=0.005 # Grid size (DF Sampling- Log)
const dlogv=0.005

const res = 50 # Resolution for coeff sampling
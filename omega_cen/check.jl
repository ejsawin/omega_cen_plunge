# Import needed packages, files
using Plots
using BenchmarkTools
using LaTeXStrings
using Measures

include("constants.jl")
include("data_io.jl")
include("interpolation.jl")
include("plummer_theory.jl")
include("psi_calc.jl")
include("sampling.jl")
include("lc_calc.jl")
include("coefficients.jl")
include("gw.jl")
include("monte_carlo.jl")
include("orbits.jl")
include("gw_experiment.jl")

# Read in file
dat = read_file("omegaCenEddieNew.csv")

# Compute potential 
psi_tab, psi_tot_tab, M_tab = find_psi_arrays(dat.r,dat.m)
psi_calc = psi_exact(dat.r,psi_tab,M_tab)

# Compute distribution function 
DF = log_rv_df(dat)
Etab = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
Ltab = dat.r .* dat.vt

# Target object mass
m_obj = 5e-6 #2.8e-7 

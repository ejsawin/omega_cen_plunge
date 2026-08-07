# Import needed packages, files
using Plots
using BenchmarkTools
using LaTeXStrings
using Measures
using HDF5
using CSV
using DataFrames
using Statistics
using QuadGK
using Plots
using Random

include("constants.jl")
include("data_io.jl")
include("interpolation.jl")
include("plummer_theory.jl")
include("psi_calc.jl")
include("rv_sampling.jl") # Defines canonical SampDF struct (needed by sampling.jl)
include("sampling.jl")
include("lc_calc.jl")
include("coefficients.jl")
include("gw.jl")
include("monte_carlo.jl")
include("orbits.jl")
include("gw_experiment.jl")
include("automation.jl") # automate_dc: snapshot -> diffusion-coefficient pipeline

# Read in file
#dat = read_file("omegaCenEddieNew.csv")
#dat = read_file("~/ToProject/Omega_Cen_Plunge/1Gyr.csv")
#dat = read_file_h5("/proj/rodriguezlab/projects/omega_cen/final_models/fb_0.02/N11e6_MBH5e2_RV5_gamma_3_alpha3_2.5_fb_0.02/output.window.snapshots.h5", "10(t=1.0002531Gyr)")
#dat = read_file_h5("/proj/rodriguezlab/projects/omega_cen/final_models/fb_0.02/N11e6_MBH5e2_RV5_gamma_3_alpha3_2.5_fb_0.02/output.window.snapshots.h5", "25(t=12.005126Gyr)")
#dat = read_file_h5("/proj/rodriguezlab/projects/Big_Cluster/10_1000_timesteplim/output_10_1000_timesteplim_2.window.snapshots.h5", "9(t=1.0000013Gyr)")

#=
psi_tab, psi_tot_tab, M_tab = find_psi_arrays(dat.r,dat.m)
psi_rtab = dat.r            # active potential grid (switch via activate_orbit_psi!)
psi_Mtot = sum(dat.m)       # total stellar mass on the potential grid
psi_calc = psi_exact(psi_rtab,psi_tab,M_tab)

# Sorted per-particle energies
ETab = sort(0.5 .* (dat.vr.^2 .+ dat.vt.^2) .+ psi_calc.(dat.r))

# Influence radius: where the enclosed stellar mass equals the IMBH mass
r_infl = dat.r[searchsortedfirst(view(M_tab, 1:length(dat.r)), M_bh)]

# Compute distribution function
DF = kde_rv_df(dat)

# Target object mass
m_obj = 5e-6
rho = rho_calc(DF,dat)

# Read in diffusion coefficients
coef = load_coeffs("DC_1dot00Gyr_powerlaw.hdf5")
=#
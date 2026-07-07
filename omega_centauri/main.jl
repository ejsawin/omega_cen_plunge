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

# Read in file
#dat = read_file("omegaCenEddieNew.csv")
dat = read_file("~/ToProject/Omega_Cen_Plunge/1Gyr.csv")

# Compute potential 
psi_tab, psi_tot_tab, M_tab = find_psi_arrays(dat.r,dat.m)
psi_calc = psi_exact(dat.r,psi_tab,M_tab)

# Compute distribution function 
DF = kde_rv_df(dat)

# Target object mass
m_obj = 5e-6
rho = rho_calc(DF,dat)

# Read in diffusion coefficients
coef = load_coeffs("Jun_10_expanded_adap_sp_forplot.hdf5")
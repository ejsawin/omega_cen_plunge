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

# Target object mass
m_obj = 5e-6
rho = rho_calc(DF,dat)

# Read in diffusion coefficients
coef = load_coeffs("coef_50.hdf5")

# Calculate energy, momentum 
Etab1 = 0.5 .* (dat.vr .^ 2 .+ dat.vt .^ 2) .+ psi_calc.(dat.r)
Ltab1 = dat.r .* dat.vt

Etab = Etab1[(Etab1 .< 0) .& (dat.startype .== 14)] # Filter unbound orbits
Ltab = Ltab1[(Etab1 .< 0) .& (dat.startype .== 14)]

jtab = zeros(length(Etab))

Threads.@threads for i in eachindex(Etab)
    jtab[i] = Ltab[i] / last(find_Lc(Etab[i]))
end

function avg_coef_ej_scan(E, j_arr)
    n = length(j_arr)
    DE   = zeros(n)
    DEE  = zeros(n)
    Dj   = zeros(n)
    Djj  = zeros(n)
    DEj  = zeros(n)
    Threads.@threads for i in 1:n
        j = j_arr[i]
        DE1, DE2, Dj1, Dj2, DEE_i, Djj_i, DEj_i = avg_coef_ej(E, j)
        DE[i]  = DE1 + m_obj * DE2
        DEE[i] = DEE_i
        Dj[i]  = Dj1 + m_obj * Dj2
        Djj[i] = Djj_i
        DEj[i] = DEj_i
    end
    return DE, DEE, Dj, Djj, DEj
end
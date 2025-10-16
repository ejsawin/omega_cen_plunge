# Necessary modules
include("make_globals.jl")
include("vel_coef.jl")
include("vel_calc.jl")
include("midpoint.jl")
include("loc_coef.jl")
include("orbit_bounds.jl")
include("period_find.jl")
include("avg_coef.jl")
include("find_Lc.jl")

using .MakeGlobals
using Plots
using KernelDensity
using HCubature
using Base.Threads
using LaTeXStrings
using BenchmarkTools 

# Read in data
read_file("omegaCenEddieNew.csv")
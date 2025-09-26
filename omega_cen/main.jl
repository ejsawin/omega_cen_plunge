# Necessary modules
include("make_globals.jl")
include("vel_coef.jl")
include("vel_calc.jl")
include("midpoint.jl")
include("loc_coef.jl")
include("orbits.jl")
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
read_file("Plummer.dat")
#read_file("omegaCenEddie.csv")


e_range=range(-0.99,-0.01,length=1000)

p=plot(e_range,find_Lc.(e_range),title=L"L_{c}")
p1=plot(e_range,dLc_dE.(e_range),title=L"\frac{dL{c}}{dE}")
p2=plot(e_range,d2Lc_dE2.(e_range),title=L"\frac{d2L{c}}{dE2}")

p_comb=plot(p,p1,p2,size=(3000,1000),layout=(1,3),xlabel="E")
display(p_comb)
readline()
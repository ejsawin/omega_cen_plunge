### Graphing / experiment functions - not needed for code ###

"""
emperical_total = sum(tabN_rv) * dr * dv
theoretical_total = sum(tabDFth_rv) * dr * dv
scale=theoretical_total/emperical_total


#graph1 = heatmap(tabr_samp, tabv_samp, scale*tabN_rv',xlim=(0,4),c=:thermal,title="N_1D (norm)")
#graph2 = heatmap(tabr_samp, tabv_samp, tabDFth_rv',xlim=(0,4),c=:thermal,title="1D DF (Theory)")
#graph3 = plot(graph1,graph2,layout=(1,2))

graph=surface(tabr_samp, tabv_samp, (scale*tabN_rv .- tabDFth_rv)',title="Difference",camera=(90,-5),xlabel="r",ylabel="v",zlabel="Difference",c=:viridis,lw=0.2)
display(graph)
readline()
"""

# Plummer density comparison
"""
function plum_rho(r)
    rho=3/(4*pi*(b^3)) * (1+(r^2)/(b^2))^(-5/2)
    return rho
end

test_r=exp.(range(log(rmin),log(2.5),length=500))
tab_rho=zeros(Float64,500)
tab_rho_th=zeros(Float64,500)

Threads.@threads for i in 1:500
    r=test_r[i]

    a=midpoint(x->N_1D(r,x),vmin,vmax,int_steps)
    tab_rho[i]=a

    tab_rho_th[i]=4*pi*(r^2)*plum_rho(r)
end


graph=scatter(test_r,tab_rho,label="Interp",ylabel=L" \rho(r) ",xlabel="r")
plot!(graph,test_r,tab_rho_th,label="Theory",linestyle=:dash)

graph2=plot(test_r,abs.(tab_rho_th .- tab_rho) ./ tab_rho_th,label="Rel Error",linecolor=:green)
graph3=plot(graph,graph2,layout=(1,2),size=(1000,500))
gui()
display(graph3)
#savefig("plum_rho_comp.png")
readline()
"""

# Plot coefficients in (vr,vt)
function graph_coef_2D()

    r_=0.1
    psi_val=psi_calc(r_)
    vmax = sqrt(abs(2*psi_val))

    vr_vals=range(0.5,1.3,length=50)
    vt_vals=range(0.5,1.3,length=50)

    nv_r=length(vr_vals)
    nv_t=length(vt_vals)

    Dpar_vals=Array{Float64}(undef,nv_r,nv_t)
    Dpar2_vals=similar(Dpar_vals)
    Dtan2_vals=similar(Dpar_vals)

    Dpar_vals_th  = similar(Dpar_vals)
    Dpar2_vals_th = similar(Dpar_vals)
    Dtan2_vals_th = similar(Dpar_vals)

    Dpar_comp  = similar(Dpar_vals)
    Dpar2_comp = similar(Dpar_vals)
    Dtan2_comp = similar(Dpar_vals)

    Threads.@threads for i in 1:nv_r
        vr = vr_vals[i]
        for j in 1:nv_t
            vt = vt_vals[j]

            # Check if physical
            if sqrt(vr^2 + vt^2) > vmax
                Dpar_vals[i,j]     = NaN
                Dpar2_vals[i,j]    = NaN
                Dtan2_vals[i,j]    = NaN

                Dpar_vals_th[i,j]  = NaN
                Dpar2_vals_th[i,j] = NaN
                Dtan2_vals_th[i,j] = NaN

                Dpar_comp[i,j]     = NaN
                Dpar2_comp[i,j]    = NaN
                Dtan2_comp[i,j]    = NaN
                continue
            end

            E_test = 0.5 * (vr^2 + vt^2) + psi_calc(r_)
            L_test = vt * r_

            d1, d2, d3    = v_coef_lin(r_, E_test, L_test)
            d11, d22, d33 = v_coef_lin_th(r_, E_test, L_test)

            Dpar_vals[i,j]     = d1
            Dpar2_vals[i,j]    = d2
            Dtan2_vals[i,j]    = d3

            Dpar_vals_th[i,j]  = d11
            Dpar2_vals_th[i,j] = d22
            Dtan2_vals_th[i,j] = d33

            Dpar_comp[i,j]  = (d1 - d11) / d11
            Dpar2_comp[i,j] = (d2 - d22) / d22
            Dtan2_comp[i,j] = (d3 - d33) / d33
        end
    end

    # plot heatmaps
    graph1 = heatmap(vt_vals, vr_vals, Dpar_vals,
                     xlabel="v_t", ylabel="v_r", title="⟨Δv∥⟩ (interp)",showemptybins=true,color=:viridis)
    graph2 = heatmap(vt_vals, vr_vals, Dpar2_vals,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv∥)²⟩ (interp)",showemptybins=true,color=:viridis)
    graph3 = heatmap(vt_vals, vr_vals, Dtan2_vals,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv⊥)²⟩ (interp)",showemptybins=true,color=:viridis)

    graph4 = heatmap(vt_vals, vr_vals, Dpar_vals_th,
                     xlabel="v_t", ylabel="v_r", title="⟨Δv∥⟩ (theory)",showemptybins=true,color=:viridis)
    graph5 = heatmap(vt_vals, vr_vals, Dpar2_vals_th,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv∥)²⟩ (theory)",showemptybins=true,color=:viridis)
    graph6 = heatmap(vt_vals, vr_vals, Dtan2_vals_th,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv⊥)²⟩ (theory)",showemptybins=true,color=:viridis)

    graph7 = heatmap(vt_vals, vr_vals, Dpar_comp,
                     xlabel="v_t", ylabel="v_r", title="⟨Δv∥⟩ Rel Error",showemptybins=true,color=:viridis)
    graph8 = heatmap(vt_vals, vr_vals, Dpar2_comp,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv∥)²⟩ Rel Error",showemptybins=true,color=:viridis)
    graph9 = heatmap(vt_vals, vr_vals, Dtan2_comp,
                     xlabel="v_t", ylabel="v_r", title="⟨(Δv⊥)²⟩ Rel Error",showemptybins=true,color=:viridis)
    graph10=plot(graph1, graph2, graph3, graph4, graph5, graph6, graph7, graph8, graph9, layout=(3,3),size=(1500,1500))
    gui()
    display(graph10)
    #savefig("diff_coef_comp.png")
    readline()
end

#savefig("plummer_locv_comp_lin.png")

# Plot coefficients in 1D 
function graph_coef_1D()
    vt_vals = range(0.5, 7.0, length=50)

    vr = 1.0
    r_ = 1.0

    n = length(vt_vals)

    # preallocate result arrays
    Dpar_vals      = Vector{Float64}(undef, n)
    Dpar2_vals     = Vector{Float64}(undef, n)
    Dtan2_vals     = Vector{Float64}(undef, n)

    Dpar_vals_th   = Vector{Float64}(undef, n)
    Dpar2_vals_th  = Vector{Float64}(undef, n)
    Dtan2_vals_th  = Vector{Float64}(undef, n)

    Dpar_comp      = Vector{Float64}(undef, n)
    Dpar2_comp     = Vector{Float64}(undef, n)
    Dtan2_comp     = Vector{Float64}(undef, n)


    Threads.@threads for i in 1:n
        vt = vt_vals[i]

        E_test = 0.5 * (vr^2 + vt^2) + psi_calc(r_)
        L_test = vt * r_

        d1, d2, d3   = vel_coef_lin(r_, E_test, L_test)
        d11, d22, d33 = vel_coef_lin_th(r_, E_test, L_test)

        Dpar_vals[i]     = d1
        Dpar2_vals[i]    = d2
        Dtan2_vals[i]    = d3

        Dpar_vals_th[i]  = d11
        Dpar2_vals_th[i] = d22
        Dtan2_vals_th[i] = d33

        Dpar_comp[i]     = (d1 ./ d11) #/ d11
        Dpar2_comp[i]    = (d2 ./ d22) #/ d22
        Dtan2_comp[i]    = (d3 ./ d33) #/ d33
    end


    graph1=plot(vt_vals, Dpar_vals, label="⟨Δv∥⟩ Itp", title="Comparison",xlabel="v_t",linestyle=:dash)
    plot!(vt_vals, Dpar2_vals, label="⟨(Δv∥)²⟩ Itp",linestyle=:dash)
    plot!(vt_vals, Dtan2_vals, label="⟨(Δv⊥)²⟩ Itp",linestyle=:dash)

    plot!(vt_vals, Dpar_vals_th, label="⟨Δv∥⟩ Th")
    plot!(vt_vals, Dpar2_vals_th, label="⟨(Δv∥)²⟩ Th")
    plot!(vt_vals, Dtan2_vals_th, label="⟨(Δv⊥)²⟩ Th")

    graph2=plot(vt_vals, Dpar_comp, label="⟨Δv∥⟩", title = "Rel Error", xlabel="v_t")
    plot!(vt_vals, Dpar2_comp, label="⟨(Δv∥)²⟩")
    plot!(vt_vals, Dtan2_comp, label="⟨(Δv⊥)²⟩")

    graph5=plot(graph1,graph2,layout=(1,2),size=(1000,500))

    gui()
    display(graph5)
    readline()
end

"""
function safe_DE(E, j)
    try
        val = find_DEj(E,j)
        return isnan(val) ? NaN : val
    catch
        return NaN
    end
end

function dee_heatmap(E_range, j_range; steps_E=40, steps_j=40)
    ee = range(E_range[1], E_range[2], length=steps_E)
    jj = range(j_range[1], j_range[2], length=steps_j)

    z = Array{Float64}(undef, length(ee), length(jj))

    Threads.@threads for i in eachindex(ee)
        for j in eachindex(jj)
            z[i,j] = abs(safe_DE(ee[i], jj[j]))
        end
    end

    p = heatmap(jj, ee, z; xlabel="j", ylabel="E", title="DEj(E, j)")

    p2 = contour(jj, ee, z; xlabel="j", ylabel="E", title="DEj(E, j)")

    p1=plot(p,p2,layout=(1,2),size=(2000,1100))
    display(p1)
    savefig("DEj_heatmap.png")
    readline()
end

# run it
dee_heatmap((-0.99, -0.01), (0.01,0.99), steps_E=50, steps_j=50)
"""
using SpecialFunctions

## Estimate velocity dispersion from surrounding objects ##
function sigma(r, rad)
    idx = find_shell(r, dat.r) # Raw indices

    i_low = max(1, idx - rad) # Ensure > rad from edge #
    i_high = min(length(dat.r), idx + rad)

    indices = [i_low:idx-1; idx+1:i_high] # Safe indices

    tabm = dat.m[indices]
    tabv = dat.v[indices]

    sigma = (1/sqrt(3)) * sqrt(sum(tabm .* (tabv .^ 2)) / sum(tabm))
    return sigma
end

## Find circular velocity at given radius ##
function vc(r)
    idx = max(1,find_shell(r,dat.r))
    Menc = M_tab[idx] # Enclosed mass 

    return sqrt((G*Menc)/r) # vc^2 / r = GM(r) / r^2
end

## Bracket function ##
function bracket(r,rad)
    X = vc(r) / (sqrt(2) * sigma(r,rad))
    return erf(X) - (2X/sqrt(pi))*exp(-X^2)
end

## Frictional force ##
function fric(r,Mbh,rad)

    coulog=log(0.11*length(dat.r)) # Coulomb logarithm

    num = 4*pi*(G^2)*(Mbh^2)*rho(r)*coulog*bracket(r,rad)
    den = vc(r)^2
    return num / den
end

## Num derivative of r*vc(r) w/ adaptive step ##
function drvc_dr(r; dr=max(r*0.05, 0.01))

    r1 = max(r - dr, 1e-6) # Protect against negative radii 
    r2 = r + dr

    return (r2 * vc(r2) - r1 * vc(r1)) / (r2 - r1)
end

## dr/dt ##
function rdot(r,Mbh,rad)
    return -(fric(r,Mbh,rad) * r) / (Mbh * drvc_dr(r))
end

## Calculate orbit decay time ##
function t_fric(ri, rf, Mbh, rad; n=1000)
    return midpoint(r -> 1/abs(rdot(r, Mbh, rad)), rf, ri, n)
end

## Calculate plunge (t,r) with Euler's method ##
function plunge(ri, rf, Mbh, rad; dt=0.1)

    r = Float64(ri)
    t = 0.0

    r_tab = [r]
    t_tab = [t]

    while r > rf

        dr = rdot(r, Mbh, rad) * dt # Evolve by timestep
        r += dr
        t += dt

        println(t, "  ", r)
        push!(r_tab, r)
        push!(t_tab, t)

    end

    # Convert to Myr, pc
    return t_tab .* 0.07453, r_tab .* 5
end

## Plot plunge in 2D ##
function plot_plunge(ri, rf, Mbh, rad; dt=0.1, filename="plunge.png")
    
    t_myr, r_pc = plunge(ri, rf, Mbh, rad; dt=dt)

    n = length(r_pc)
    θ = zeros(n)

    for i in 2:n
        
        r_henon = r_pc[i] / 5.0 # Convert into physical units 
        dt_henon = (t_myr[i] - t_myr[i-1]) / 0.07453

        θ[i] = θ[i-1] + (vc(r_henon) / r_henon) * dt_henon # Assuming circular velocity, same as code
    end

    x = r_pc .* cos.(θ)
    y = r_pc .* sin.(θ)

    plot(x, y, color=:cyan, linewidth=0.8, aspect_ratio=:equal,
         background=:black, foreground=:white,
         xlabel="x (pc)", ylabel="y (pc)", legend=false,dpi=300)
    scatter!([0], [0], color=:gold, markersize=2)
    savefig(filename)
end
# --- Simple midpoint integrator ---

function midpoint(func,x_start,x_end,steps)

    delx=(x_end-x_start)/steps
    integral=0.0

    for i in 0:steps-1
        #Define subinterval
        x0=x_start+i*delx 
        x1=x0+delx

        #Find midpoint
        mid=(x0+x1)/2

        integral+=func(mid)*delx 
    end

    return integral 
end
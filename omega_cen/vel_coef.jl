# --- Calculate velocity diffusion coefficients at (r, E, L) ---
function vel_coef_lin(rr,EE,LL) # Eq 1.32 - 1.36

    # determine velocity
    v_step=sqrt((v_r(rr,EE,LL))^2 + (v_t(rr,EE,LL))^2)

    # constants
    c1=-(G^2*coulomb_log)/(v_step^2 * rr^2)
    c2=(2*G^2*coulomb_log)/(3*rr^2)

    # integrands
    vp_int(vv)=F_1D(rr,vv)
    vp_int1(vv)=N_1D(rr,vv)

    vpsq_int(vv)=(vv^2)/(v_step^3)*F_1D(rr,vv)
    vpsq_int1(vv)=(1 / vv) * F_1D(rr,vv)

    # integrals
    vp_calc=midpoint(x->vp_int(x),vmin,v_step,int_steps)
    vp_calc1=midpoint(x->vp_int1(x),vmin,v_step,int_steps)

    vpsq_calc=midpoint(x->vpsq_int(x),vmin,v_step,int_steps)
    vpsq_calc1=midpoint(x->vpsq_int1(x),v_step,vmax,int_steps)

    # calculate diffusion coefficients 
    D_par=c1*(vp_calc+m_test*vp_calc1)
    D_par2=c2*(vpsq_calc+vpsq_calc1)
    D_tan2=c2*((3/v_step)*vp_calc - vpsq_calc + 2 * vpsq_calc1)

    return D_par,D_par2,D_tan2
end

function vel_coef_log(rr,EE,LL) # Eq 1.44 - 1.48

    # convert r, v into log 
    logr=log(rr)

    v_step=sqrt((v_r(rr,EE,LL))^2 + (v_t(rr,EE,LL))^2)
    logv=log(v_step)

    # constants
    c1 = -(G^2 * coulomb_log) / (v_step^2 * rr^3)
    c2 = (2 * G^2 * coulomb_log) / (3 * rr^3)

    # integrands
    vp_int(log_va) = F_1D(logr,log_va) # 1D log DF
    vp_int1(log_va) = N_1D(logr,log_va)

    vpsq_int(log_va) = (((exp(log_va))^2) / (v_step^3)) * F_1D(logr,log_va)
    vpsq_int1(log_va) = (1 / exp(log_va)) * F_1D(logr,log_va)

    # integrals 
    vp_calc=midpoint(x->vp_int(x),log(vmin),logv,int_steps)
    vp_calc1=midpoint(x->vp_int1(x),log(vmin),logv,int_steps)

    vpsq_calc=midpoint(x->vpsq_int(x),log(vmin),logv,int_steps)
    vpsq_calc1=midpoint(x->vpsq_int1(x),logv,log(vmax),int_steps)

    # calculate diffusion coefficients 
    D_par=c1*(vp_calc+m_test*vp_calc1)
    D_par2=c2*(vpsq_calc+vpsq_calc1)
    D_tan2=c2*((3/v_step)*vp_calc - vpsq_calc + 2 * vpsq_calc1)

    return D_par,D_par2,D_tan2
end

# --- Theoretical velocity coefficients ---
function vel_coef_th(rr,EE,LL) # Eq 1.22 - 1.25

    # determine velocity
    v_step=sqrt((v_r(rr,EE,LL))^2 + (v_t(rr,EE,LL))^2)

    # constants
    c1=-(G^2*(m_test+m_field)*coulomb_log)/(v_step^2*rr^2)
    c2=(2*G^2*m_field*coulomb_log)/(3*rr^2)

    # integrands
    vp_int(vv) = F_th(rr,vv)
    vpar2_int1(vv) = (vv^2/v_step^3)*F_th(rr,vv)
    vpar2_int2(vv) = (1/vv)*F_th(rr,vv)
    vtan2_int1(vv) = ((3/v_step) - (vv^2/v_step^3))*F_th(rr,vv)
    vtan2_int2(vv) = (1/vv)*F_th(rr,vv)

    # calculate integrals
    vp_calc = midpoint(x->vp_int(x),vmin,v_step,int_steps)
    vpar2_calc1 = midpoint(x->vpar2_int1(x),vmin,v_step,int_steps)
    vpar2_calc2 = midpoint(x->vpar2_int2(x),v_step,vmax,int_steps)
    vtan2_calc1 = midpoint(x->vtan2_int1(x),vmin,v_step,int_steps)
    vtan2_calc2 = midpoint(x->vtan2_int2(x),v_step,vmax,int_steps)

    # calculate diffusion coefficients
    D_par = c1*vp_calc
    D_par2 = c2*(vpar2_calc1+vpar2_calc2)
    D_tan2 = c2*(vtan2_calc1+2*vtan2_calc2)

    return D_par,D_par2,D_tan2
end
# --- Calculate local diffusion coefficients ---

function loc_coef(rr,EE,LL)
    v_step=sqrt((v_r(rr,EE,LL)) ^ 2 + (v_t(rr,EE,LL)) ^ 2)

    vp,vp2,vt2=vel_coef_lin(rr,EE,LL)

    del_E=0.5*vp2+0.5*vt2+v_step*vp # Eq 2.45
    del_E2=(v_step^2)*vp2 # Eq 2.46

    del_L=(LL/v_step)*vp+((rr^2)/(4*LL))*vt2 # Eq 2.47
    del_L2=((LL^2)/(v_step^2))*vp2+0.5*(rr^2-(LL/v_step)^2)*vt2 # Eq 2.48

    del_EL=LL*vp2 # Eq 2.49

    return del_E, del_E2, del_L, del_L2, del_EL
end 
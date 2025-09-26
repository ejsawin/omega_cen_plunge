# --- Orbit averaged diffusion coefficients ---
function avg_coef_int(tt,EE,LL,ii,rp,ra)
    r_hat=(ra-rp)/2 # Eq 1.18
    t_hat=(ra+rp)/2 # Eq 1.19

    rrr=r_hat*sin(tt)+t_hat # Eq 1.17
    vr_tilde=v_r(rrr,EE,LL)/(r_hat*cos(tt)) # Eq 1.20

    return (loc_coef(rrr,EE,LL)[ii])/vr_tilde
end

function avg_coef(EE_,LL_,indx) # Eq 1.59 - 1.63
    rp,ra=compute_rp_ra(EE_,LL_)
    T=period(EE_,LL_)
    return (2/T)*midpoint(x -> avg_coef_int(x,EE_,LL_,indx,rp,ra),-pi/2,pi/2,int_steps)
end

### D(E,L) -> D(E,j)
function find_DE(E,j)
    Jc=find_Lc(E)
    J=j*Jc
    return avg_coef(E,J,1)
end

function find_DEE(E,j)
    Jc=find_Lc(E)
    J=j*Jc
    return avg_coef(E,J,2)
end

function find_Dj(E,j) # Eq 1.67
    Jc=find_Lc(E)
    J=j*Jc

    DE=avg_coef(E,J,1)
    DEE=avg_coef(E,J,2)
    DJ=avg_coef(E,J,3)
    DEJ=avg_coef(E,J,5)

    dJc=dLc_dE(E)
    dJc2=d2Lc_dE2(E)

    dj_1 = (DJ/Jc) - (J/(Jc^2)) * dJc * DE - (1/(Jc^2)) * dJc * DEJ
    dj_2 = -0.5 * ((1/(Jc^2)) * dJc2 - (2/(Jc^3)) * (dJc^2)) * DEE

    return dj_1+dj_2
end

function find_DEj(E,j) # Eq 1.69
    Jc = find_Lc(E)
    J  = j*Jc

    DEE = avg_coef(E,J,2)
    DEJ = avg_coef(E,J,5)

    dJc  = dLc_dE(E)

    dej = -(J/(Jc^2)) * dJc * DEE + (DEJ/Jc)

    return dej
end

function find_Djj(E,j) # Eq 1.70
    Jc = find_Lc(E)
    J  = j*Jc

    DEE = avg_coef(E,J,2)
    DJJ = avg_coef(E,J,4)
    DEJ = avg_coef(E,J,5)

    dJc = dLc_dE(E)

    djj_1 = ((J/(Jc^2))*dJc)^2 * DEE
    djj_2 = -2 * (J/(Jc^3)) * dJc * DEJ + DJJ/(Jc^2)

    return djj_1 + djj_2
end

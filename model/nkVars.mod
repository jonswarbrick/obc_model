parameters Z_bar Disp_bar MC_bar C_bar H_bar Y_bar y_bar GSS;

Z_bar = 1 / betta - ( 1 -  deltaSS );
Disp_bar = ( (1-xi)*( ( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)) )^(-zzeta) )/(1-xi*Pi_bar^zzeta);
MC_bar =  ( ( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)) )*((zzeta-1)/zzeta)*(1-xi*betta*Pi_bar^zzeta)/(1- xi*betta*Pi_bar^(zzeta-1));

@#if utility_type == 1
    H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH))* ( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ) / (MC_bar*(1-alp)) );
@#endif
@#if utility_type == 2   
    H_bar  = call_csolve_nk;
@#endif
@#if utility_type == 3   
    H_bar  =  ( 1 / ( (1 - epsilonH)^(psi_h) * (1 - epsilonC)*( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ) / (MC_bar*(1-alp)) ) )^(1/(psi_h+1));
@#endif
@#if utility_type == 4   
    H_bar = 1 / ( 1 + varrho/(1-varrho) * ( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ) / (MC_bar*(1-alp)) );
@#endif
@#if utility_type == 5
    H_bar = (((1-alp)*((MC_bar*alp/(1 / betta - ( 1 -  deltaSS )))^(alp/(1-alp)))^(1-gam_jr) ) / ( varrho*( (1-alp)*gam_jr*(( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ))^(gam_jr-1) + (theta_jr+1-gam_jr)*(( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ))^gam_jr) ) )^(1/theta_jr);
@#endif


C_bar = ( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar )*( Ass * H_bar ) * ( MC_bar*alp/Z_bar ) ^ ( alp / ( 1 - alp ) );
Y_bar = ( Ass * H_bar ) * ( MC_bar*alp/Z_bar ) ^ ( alp / ( 1 - alp ) ) / Disp_bar;
y_bar = log(Y_bar);
GSS = gySS*Y_bar;

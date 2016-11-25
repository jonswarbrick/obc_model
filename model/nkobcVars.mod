//required
var b D mv nu kappa;
@#for lag in [1:7]
var prodR@{lag} lagD@{lag} SZ@{lag};
@#endfor
//analysis
var E B;

parameters y_bar Disp_bar MC_bar C_bar C_by_YW YW_bar GSS nubar kappa_new H_bar;
   
nubar = parameter_nubar;
kappa_new = parameter_kappa_new;

parameters lambdaB_bar SZ7_bar SZ6_bar SZ5_bar SZ4_bar SZ3_bar SZ2_bar 
SZ1_bar MD1_bar MD2_bar MD3_bar MD4_bar MD5_bar MD6_bar MD7_bar MV_bar;

lambdaB_bar = gam*(1-(1-gam)*(1-Theta));
SZ7_bar = lambdaB_bar;
SZ6_bar = lambdaB_bar + (1-gam)*SZ7_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ5_bar = lambdaB_bar + (1-gam)*SZ6_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ4_bar = lambdaB_bar + (1-gam)*SZ5_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ3_bar = lambdaB_bar + (1-gam)*SZ4_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ2_bar = lambdaB_bar + (1-gam)*SZ3_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ1_bar = lambdaB_bar + (1-gam)*SZ2_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
@#for lag in [1:7]
MD@{lag}_bar = SZ@{lag}_bar*(Pi_bar/ betta)^@{lag}*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_bar );
@#endfor
MV_bar = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam); 
MC_bar =  (( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)))*((zzeta-1)/zzeta)*(1-xi*((1-gam)*betta*MV_bar)*Pi_bar^zzeta)/(1- xi*((1-gam)*betta*MV_bar)*Pi_bar^(zzeta-1));
Disp_bar = ( (1-xi)*(( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)))^(-zzeta) )/(1-xi*Pi_bar^zzeta);
C_by_YW = 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/((1/((1-gam)*(betta/Pi_bar)*MV_bar))/Pi_bar-1+deltaSS);

@#if utility_type == 1
    H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH))* C_by_YW / (MC_bar*(1-alp)) );
@#endif
@#if utility_type == 2   
    H_bar  = call_csolve_nk_obc;
@#endif
@#if utility_type == 3   
    H_bar  =  ( 1 / ( (1 - epsilonH)^(psi_h) * (1 - epsilonC)*C_by_YW / (MC_bar*(1-alp)) ) )^(1/(psi_h+1));
@#endif
@#if utility_type == 4   
    H_bar = 1 / ( 1 + varrho/(1-varrho) * C_by_YW / (MC_bar*(1-alp)) );
@#endif
@#if utility_type == 5
    H_bar = (((1-alp)*((MC_bar*alp/((1/((1-gam)*(betta/Pi_bar)*MV_bar))/Pi_bar-1+deltaSS))^(alp/(1-alp)))^(1-gam_jr) ) / ( varrho*( (1-alp)*gam_jr*(C_by_YW)^(gam_jr-1)+  (theta_jr+1-gam_jr)*(C_by_YW)^gam_jr) ) )^(1/theta_jr);
@#endif


YW_bar = ( Ass * H_bar ) * ( MC_bar*alp/((1/((1-gam)*(betta/Pi_bar)*MV_bar))/Pi_bar-1+deltaSS) ) ^ ( alp / ( 1 - alp ) );
GSS = gySS*YW_bar/Disp_bar;
C_bar = C_by_YW*YW_bar;
y_bar = log( YW_bar/Disp_bar );

//required
var r_PLUS_b E mK jK jB;
//analysis
var D;

parameters C_bar H_bar GSS;

@#if utility_type == 1
    H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH)) * (1 - gySS - deltaSS * alp/( 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) ) )) / (1-alp) );
@#endif
@#if utility_type == 2 
    H_bar  = call_csolve_obc;
@#endif
@#if utility_type == 3 
    H_bar  =  ( 1 / ( (1 - epsilonH)^(psi_h) * (1 - epsilonC)*(1 - gySS - deltaSS * alp/( 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) ) )) / (1-alp) ) )^(1/(psi_h+1));
@#endif
@#if utility_type == 4
    H_bar = 1 / ( 1 + varrho/(1-varrho) * (1 - gySS - deltaSS * alp/( 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) ) )) / (1-alp) );
@#endif

GSS = gySS*( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) );    

C_bar = (1 - gySS - deltaSS * alp/( 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) ) )) * ( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) );

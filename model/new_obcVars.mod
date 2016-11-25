//required
var b D mv  nu kappa;
@#for lag in [1:7]
var prodR@{lag} lagD@{lag} SZ@{lag};
@#endfor
//analysis
var E_rate D_rate;


parameters C_bar H_bar GSS nubar kappa_new;
   
nubar = parameter_nubar;
kappa_new = parameter_kappa_new;

parameters lambdaB_bar SZ7_bar SZ6_bar SZ5_bar SZ4_bar SZ3_bar SZ2_bar 
SZ1_bar MD1_bar MD2_bar MD3_bar MD4_bar MD5_bar MD6_bar MD7_bar;

lambdaB_bar = gam*(1-(1-gam)*(1-Theta));
SZ7_bar = lambdaB_bar;
SZ6_bar = lambdaB_bar + (1-gam)*SZ7_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ5_bar = lambdaB_bar + (1-gam)*SZ6_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ4_bar = lambdaB_bar + (1-gam)*SZ5_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ3_bar = lambdaB_bar + (1-gam)*SZ4_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ2_bar = lambdaB_bar + (1-gam)*SZ3_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ1_bar = lambdaB_bar + (1-gam)*SZ2_bar*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
@#for lag in [1:7]
MD@{lag}_bar = SZ@{lag}_bar*(1 / betta)^@{lag}*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_bar );
@#endfor

@#if utility_type == 1
    H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH))* ( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) / (1-alp) );
@#endif
@#if utility_type == 2   
    H_bar  = call_csolve_new_obc;
@#endif
@#if utility_type == 3   
    H_bar  =  ( 1 / ( (1 - epsilonH)^(psi_h) * (1 - epsilonC)*( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) / (1-alp) ) )^(1/(psi_h+1));
@#endif
@#if utility_type == 4   
    H_bar = 1 / ( 1 + varrho/(1-varrho) * ( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) / (1-alp) );
@#endif
@#if utility_type == 5
    H_bar = (((1-alp)*((alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ))^(alp/(1-alp)))^(1-gam_jr) ) / ( varrho*( (1-alp)*gam_jr*(( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ))^(gam_jr-1) + (theta_jr+1-gam_jr)*(( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ))^gam_jr) ) )^(1/theta_jr);
@#endif

GSS = gySS*( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) ); 

C_bar = ( 1 - gySS - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) * ( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) );



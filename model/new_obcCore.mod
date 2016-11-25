
# lag_R = exp( r(-1) );
# lead_R = exp( r(+1) );
# SG = 0;
# lead_SG = 0;
# MV = exp( mv );
# lead_MV = exp( mv(+1) );
# lead_kappa = kappa(+1);
# lead_Xi = (1-gam)*lead_Lambda*lead_MV*(1-kappa)/(1-lead_kappa); 
# lambdaB = (1-(1-gam)*(1-Theta))*(MV-1)/(MV-(1-kappa)*(1-(1-gam)*(1-Theta)));
# lead_lambdaB = (1-(1-gam)*(1-Theta))*(lead_MV-1)/(lead_MV-(1-lead_kappa)*(1-(1-gam)*(1-Theta)));

prodR1 = lag_R;
prodR2 = lag_R*prodR1(-1);
prodR3 = lag_R*prodR2(-1);
prodR4 = lag_R*prodR3(-1);
prodR5 = lag_R*prodR4(-1);
prodR6 = lag_R*prodR5(-1);
prodR7 = lag_R*prodR6(-1);
lagD1 = D(-1)/Pi;
lagD2 = lagD1(-1)/Pi;
lagD3 = lagD2(-1)/Pi;
lagD4 = lagD3(-1)/Pi;
lagD5 = lagD4(-1)/Pi;
lagD6 = lagD5(-1)/Pi;
lagD7 = lagD6(-1)/Pi;
SZ1 = lambdaB + (1-gam)*SZ2(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ2 = lambdaB + (1-gam)*SZ3(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ3 = lambdaB + (1-gam)*SZ4(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ4 = lambdaB + (1-gam)*SZ5(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ5 = lambdaB + (1-gam)*SZ6(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ6 = lambdaB + (1-gam)*SZ7(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ7 = lambdaB; 
@#for lag in [1:7]
# MD@{lag} = SZ@{lag}*prodR@{lag}*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB );
# mD@{lag} = ( MD@{lag} + (1-gam)*(1-Theta)*prodR@{lag} )/( 1 - (1-gam)*(1-Theta) );
@#endfor
# mV = MV / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
# summ_times_D = mD1*lagD1 + mD2*lagD2 + mD3*lagD3 + mD4*lagD4 + mD5*lagD5 + mD6*lagD6 + mD7*lagD7;
# sumM_times_D = MD1*lagD1 + MD2*lagD2 + MD3*lagD3 + MD4*lagD4 + MD5*lagD5 + MD6*lagD6 + MD7*lagD7;
# sumprodR_times_D = prodR1*lagD1 + prodR2*lagD2 + prodR3*lagD3 + prodR4*lagD4 + prodR5*lagD5 + prodR6*lagD6 + prodR7*lagD7;
# lead_MD1 = SZ1(+1)*prodR1(+1)*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lead_lambdaB );
# B = exp( b );
# lag_B = exp( b(-1) );
# Vhat = ( RK*lag_S/Pi - (lag_R-SG)*lag_B/Pi )/(1-kappa);
# V = MV*Vhat + sumM_times_D;
# E = ( D + S - B )/(1-kappa) - Vhat;
kappa = kappa_new*( 1 - exp(-nu*E/V ) );
# lambdaD = kappa - (1-gam)*(1-kappa)*(lead_Lambda/lead_Pi)*lead_MD1;
# N = S - B;
D_rate = D/N;
E_rate = E/N;

(1-gam)*(lead_Lambda/lead_Pi)*( lead_MV/MV )*(R - lead_SG)*(1-kappa)/(1-kappa(+1)) + (MV-1) / ( MV - (1-kappa)*(1-(1-gam)*(1-Theta)) ) = 1;
1 = lead_Xi*lead_RK/lead_Pi;

0 = min( D , lambdaD );
0 = min( mV*Vhat + summ_times_D - B , lambdaB);
%%0 = min( SG , V + (1-Theta)*sumprodR_times_D );

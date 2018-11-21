var k a logit_delta q c r psi spread inv h b D mv nu kappa Xjr lambdaX;

@#for lag in [1:7]
var prodR@{lag} lagD@{lag} SZ@{lag};
@#endfor

varexo epsA;

parameters varrho alp zzeta betta deltaSS sigma_c gy rhodelta rhoA
rho_psi Ass Phi sigmaB xiB Theta kappaSS logit_deltaSS epsilon kappa_GK gam
sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr C_bar H_bar 
nubar kappa_new lambdaB_bar SZ7_bar SZ6_bar SZ5_bar SZ4_bar SZ3_bar 
SZ2_bar SZ1_bar MD1_bar MD2_bar MD3_bar MD4_bar MD5_bar MD6_bar MD7_bar;

load('../opts.mat');

varrho=2.6;
alp=0.3;
zzeta=7.0;
betta=0.995;
deltabar=0.025;
logit_deltabar = log(deltabar/(1-deltabar));
sigma_c=2;
gam_jr = 0.001;
theta_jr = 1.4;
sigma_h=2.37;
psi_h=1;
chi = .7;
epsilonC = 0;
epsilonH = 0;
Ass=1;
rhoA=parameter_rhoA;
rhodelta=.85;
rho_psi = 0;
sigmaB=0.975;
xiB=0.003;
epsilon = -2;
kappa_GK = 13;
gam=1e-8;
kappaSS=.05;
gy = 0;

Theta=parameter_Theta;
Phi=2;
sigma_a = parameter_sigma_a;
sigma_psi = 0;
sigma_delta = .9;

logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));

nubar = 400;
kappa_new = .1;

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

H_bar = csolve_H_obc;
C_bar = ( 1 - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) * ( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) );

model;

# delta = 1/(1+exp(-logit_delta));
# lead_delta = 1/(1+exp(-logit_delta(+1)));
# lag_delta = 1/(1+exp(-logit_delta(-1)));
# K = exp( k );
# lead_K = exp( k(+1) );
# lag_K = exp( k(-1) );
# A = exp( a );
# lag_C = exp( c(-1) );
# C = exp( c );
# lead_C = exp( c(+1) );
# I = exp( inv );
# lag_I = exp( inv(-1) );
# lead_I = exp( inv(+1) );
# Q = exp( q );
# lead_Q = exp( q(+1) );
# lag_Q = exp( q(-1) );
# Y = C+I;
# lead_Y = lead_C+lead_I;
# y = log(Y);
# Z = alp*Y/(psi*lag_K);
# lead_Z = alp*lead_Y/(psi(+1)*K);
# H = exp( h );
# lead_H = exp( h(+1) );
# UH = - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * theta_jr*varrho*Xjr*H^(theta_jr-1);
# UX =  - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
# UC = (C - varrho*H^theta_jr*Xjr)^(-sigma_c);
# lead_UC = (lead_C - varrho*lead_H^theta_jr*Xjr(+1))^(-sigma_c);
# lambdaC = UC + lambdaX*gam_jr*Xjr/C;
# lead_lambdaC = lead_UC + lambdaX(+1)*gam_jr*Xjr(+1)/lead_C;
# lead_Lambda = betta*lead_lambdaC/lambdaC;
# w = log(1-alp) + y - h;
# R = exp( r );
# S = Q*K;
# lag_S = lag_Q*lag_K;
# RK = psi*(Z + (1-delta)*Q)/lag_Q;
# lead_RK = psi(+1)*(lead_Z + (1-lead_delta)*lead_Q)/Q;

# lag_R = exp( r(-1) );
# lead_R = exp( r(+1) );
# SG = 0;
# lead_SG = 0;
# MV = exp( mv );
# lead_MV = exp( mv(+1) );
# lead_kappa = kappa(+1);
# lead_Xi = (1-gam)*lead_Lambda*lead_MV*(1-kappa)/(1-lead_kappa);
# lambdaB = (1-(1-gam)*(1-Theta))*(1-kappa)*(MV-1)/(MV-(1-kappa)*(1-(1-gam)*(1-Theta)));
# lead_lambdaB = (1-(1-gam)*(1-Theta))*(1-lead_kappa)*(lead_MV-1)/(lead_MV-(1-lead_kappa)*(1-(1-gam)*(1-Theta)));

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
# Vhat = ( RK*lag_S - (lag_R-SG)*lag_B )/(1-kappa);
# V = MV*Vhat + sumM_times_D;
# E = ( D + S - B )/(1-kappa) - Vhat;
# lambdaD = kappa - (1-gam)*(1-kappa)*(lead_Lambda)*lead_MD1;
# N = S - B;
# D_rate = D/N;
# E_rate = E/N;
# lev = k - log(Vhat);

prodR1 = lag_R;
prodR2 = lag_R*prodR1(-1);
prodR3 = lag_R*prodR2(-1);
prodR4 = lag_R*prodR3(-1);
prodR5 = lag_R*prodR4(-1);
prodR6 = lag_R*prodR5(-1);
prodR7 = lag_R*prodR6(-1);
lagD1 = D(-1);
lagD2 = lagD1(-1);
lagD3 = lagD2(-1);
lagD4 = lagD3(-1);
lagD5 = lagD4(-1);
lagD6 = lagD5(-1);
lagD7 = lagD6(-1);
SZ1 = lambdaB + (1-gam)*SZ2(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ2 = lambdaB + (1-gam)*SZ3(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ3 = lambdaB + (1-gam)*SZ4(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ4 = lambdaB + (1-gam)*SZ5(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ5 = lambdaB + (1-gam)*SZ6(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ6 = lambdaB + (1-gam)*SZ7(+1)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lead_lambdaB );
SZ7 = lambdaB;

kappa = kappa_new*( 1 - exp(-nu*E/V ) );

(1-gam)*(lead_Lambda)*( lead_MV/MV )*(R - lead_SG)*(1-kappa)/(1-kappa(+1)) + (MV-1) / ( MV - (1-kappa)*(1-(1-gam)*(1-Theta)) ) = 1;
1 = lead_Xi*lead_RK;
Xjr = C^gam_jr * Xjr(-1)^(1-gam_jr);
lambdaX = UX + betta * (1-gam_jr) * lambdaX(+1) * Xjr(+1) / Xjr;

0 = min( D , lambdaD );
0 = min( mV*Vhat + summ_times_D - B , lambdaB);

nu = nubar;

spread = lead_RK - R;

I*(1-Phi*(1-I/lag_I)^2) = K - (1-delta)*psi*lag_K;
1=Q*(1-Phi*(I/lag_I-1)^2 - (I/lag_I)*2*Phi*(I/lag_I-1))+lead_Lambda*2*Phi*(lead_I/I-1)*(lead_I/I)^2*lead_Q;

lead_Lambda*R=1;
y = (1-alp)*((a) + h ) + alp*( k(-1) + log(psi));

log( -UH/lambdaC) = w;

a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
logit_delta = logit_deltabar;
psi = 1;

end;

steady_state_model;
    psi = 1;
    logit_delta = logit_deltaSS;
    nu = nubar;
    q = 0;
    R_ = 1 / betta;
    r = log( R_ );
    a = log( Ass );
    lambdaD_ = 0;
    lambdaB_ = gam*(1-(1-gam)*(1-Theta));
    SZ7 = lambdaB_;
    SZ6 = lambdaB_ + (1-gam)*SZ7*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ5 = lambdaB_ + (1-gam)*SZ6*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ4 = lambdaB_ + (1-gam)*SZ5*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ3 = lambdaB_ + (1-gam)*SZ4*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ2 = lambdaB_ + (1-gam)*SZ3*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ1 = lambdaB_ + (1-gam)*SZ2*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    @#for lag in [1:7]
    prodR@{lag} = R_^@{lag};
    MD@{lag}_ = SZ@{lag}*R_^@{lag}*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_ );
    mD@{lag}_ = ( MD@{lag}_ + (1-gam)*(1-Theta)*R_^@{lag} )/( 1 - (1-gam)*(1-Theta) );
    @#endfor
    kappa = (1-gam)*betta*MD1_/(1+(1-gam)*betta*MD1_);
    MV_ = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam);
    mv = log( MV_ );
    mV_ = MV_ / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
    RK_ = 1/((1-gam)*betta*MV_);
    Z_ = RK_-1+deltaSS;
    K_over_Y = alp/Z_;
    I_over_Y = deltaSS * K_over_Y;
    C_over_Y = 1 - I_over_Y;
    h = log( H_bar );
    H_ = H_bar;
    Y_ = ( Ass * H_ ) * K_over_Y ^ ( alp / ( 1 - alp ) );
    K_ = K_over_Y * Y_;
    k = log( K_ );
    inv = log(deltaSS) + k;
    C_ = C_over_Y * Y_;
    c = log( C_ );
    Xjr = C_;
    UX_ =  - (C_ - varrho*H_^theta_jr*Xjr)^(-sigma_c) * varrho * H_^(theta_jr);
    lambdaX = UX_ / ( 1 - betta*( (1-gam_jr) ) );
    sum_mD = mD1_ + mD2_ + mD3_ + mD4_ + mD5_ + mD6_ + mD7_;
    sum_MD = MD1_ + MD2_ + MD3_ + MD4_ + MD5_ + MD6_ + MD7_;
    AUXBD = sum_mD/(1+mV_*R_/(1-kappa));
    AUXBS = (mV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa));
    AUXED = (1/nubar)*( log( 1 - kappa/kappa_new ) )*( ((MV_*R_/(1-kappa))/(1+mV_*R_/(1-kappa)))*sum_mD - sum_MD );
    AUXES = (1/nubar)*( log( 1 - kappa/kappa_new ) )*((MV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa)));
    D = K_*(RK_ - 1 - (R_-1)*AUXBS - (1-kappa)*AUXES)/(1-(1-kappa)*AUXED + (R_-1)*AUXBD);
    E_ = AUXED*D - AUXES*K_;
    B_ = AUXBD*D + AUXBS*K_;
    b = log( B_ );
    @#for lag in [1:7]
    lagD@{lag} = D;
    @#endfor
    spread = RK_ - R_;

end;



steady;
check;

shocks;
    var epsA = 1; 
end;

stoch_simul( order = 2, irf = 0, periods = 10000 ) y inv h spread D_rate E_rate lev;

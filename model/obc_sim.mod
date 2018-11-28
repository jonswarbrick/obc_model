var k a logit_delta q c r psi spread inv h b D_star md mv kappa Xjr lambdaX;

varexo epsA;

parameters varrho alp zzeta betta deltaSS sigma_c gy rhodelta rhoA
rho_psi Ass Phi sigmaB xiB Theta kappaSS logit_deltaSS epsilon kappa_GK gam
sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr C_bar H_bar 
nubar kappa_new rho;


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
rhoA=.95;
rhodelta=.85;
rho_psi = 0;
sigmaB=0.975;
xiB=0.003;
epsilon = -2;
kappa_GK = 13;
gam=1e-8;
kappaSS=.05;
gy = 0;
rho = betta * ( 1 - 1 / ( 1 + betta^(-1) + betta^(-2) + betta^(-3) + betta^(-4) + betta^(-5) + betta^(-6) + betta^(-7) ) );  

Theta=0.669734151867773;
Phi=2;
sigma_a = 0.0060729167797153;
sigma_psi = 0;
sigma_delta = .9;

logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));

nubar = 400;
kappa_new = .1;

H_bar = csolve_H_obc;
C_bar = ( 1 - deltaSS * alp/(((1/(betta*( 1 - gam*(1-(1-gam)*(gam*(1-Theta)*rho/(1-rho))/(1+(1-gam)*(gam*(1-Theta)*rho/(1-rho))))*(1-(1-gam)*(1-Theta)) )))-1+deltaSS )) ) * ( Ass * H_bar )*(alp/((1/(betta*( 1 - gam*(1-(1-gam)*(gam*(1-Theta)*rho/(1-rho))/(1+(1-gam)*(gam*(1-Theta)*rho/(1-rho))))*(1-(1-gam)*(1-Theta)) )))-1+deltaSS ))^( alp / ( 1 - alp ) );

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
# MV = exp( mv );
# lead_MV = exp( mv(+1) );
# MD = exp( md );
# lead_MD = exp( md(+1) );
# lead_kappa = kappa(+1);
# lead_Xi = (1-gam)*lead_Lambda*lead_MV*(1-kappa)/(1-lead_kappa);
# lambdaB = lead_Xi*(lead_RK-R)/(1-kappa);
# mV = MV/(1-(1-gam)*(1-Theta)) - (1-kappa);
# mD = ( MD + (1-gam)*(1-Theta)*rho)/(1-(1-gam)*(1-Theta));

# B = exp( b );
# lag_B = exp( b(-1) );
# Vhat = ( RK*lag_S - lag_R*lag_B )/(1-kappa);
# V = MV*Vhat + MD*D_star(-1);
# D = D_star/R - rho*D_star(-1);
# E = ( D + S - B )/(1-kappa) - Vhat;
# lambdaE = 1 - lead_Xi * lead_RK;
# lambdaD = kappa - lambdaE - (1-gam)*(1-kappa)*lead_Lambda*lead_MD*R;
# N = S - B;
# D_rate = D/Vhat;
# E_rate = E/Vhat;
# lev = q + k - log(Vhat);

kappa = kappa_new*( 1 - exp(-nubar*E/V ) );

MD = rho*(1-gam)*(  (1-Theta)*lambdaB +  lead_Lambda * lead_MD * R *( 1-(1-gam)*(1-Theta) ) ) / ( 1-(1-gam)*(1-Theta) - lambdaB );
(1-gam)*(lead_Lambda)*( lead_MV/MV )*R*(1-kappa)/(1-kappa(+1)) + lambdaB / (1-(1-gam)*(1-Theta)) = 1;
Xjr = C^gam_jr * Xjr(-1)^(1-gam_jr);
lambdaX = UX + betta * (1-gam_jr) * lambdaX(+1) * Xjr(+1) / Xjr;

0 = min( D_rate , lambdaD );
%0 = min( E_rate , lambdaE );
lambdaE = 0;
0 = min( mV*Vhat + mD*D_star(-1) - B , lambdaB);

spread = lead_RK - R;

I*(1-Phi*(1-I/lag_I)^2) = K - (1-delta)*psi*lag_K;
1=Q*(1-Phi*(I/lag_I-1)^2 - (I/lag_I)*2*Phi*(I/lag_I-1))+lead_Lambda*2*Phi*(lead_I/I-1)*(lead_I/I)^2*lead_Q;

lead_Lambda*R=1;
y = (1-alp)*((a) + h ) + alp*( k(-1) + log(psi));

log( -UH/lambdaC ) = w;

a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
logit_delta = logit_deltabar;
psi = 1;

end;

steady_state_model;
    psi = 1;
    logit_delta = logit_deltaSS;
    q = 0;
    R_ = 1 / betta;
    r = log( R_ );
    a = log( Ass );
    lambdaD_ = 0;
    lambdaB_ = gam*(1-(1-gam)*(1-Theta));
    MD_ = gam*(1-Theta)*rho/(1-rho);
    kappa = (1-gam)*MD_/(1+(1-gam)*MD_);
    MV_ = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam);
    mv = log( MV_ );
    mV_ = MV_ / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
    mD_ = ( MD_ + (1-gam)*(1-Theta)*rho)/(1-(1-gam)*(1-Theta));
    md = log( MD_ );
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
    D_star = K_ * ( ( mV_ + 1-kappa - (1-kappa)*(1/nubar)*(log(1-kappa/kappa_new))*MV_ )*RK_/(1-kappa+R_*mV_) - 1) / ( betta-rho - mD_ + (1-kappa)*(1/nubar)*(log( 1 - kappa/kappa_new ))*MD_ + ( mV_+1-kappa-(1-kappa)*(1/nubar)*log( 1 - kappa/kappa_new )*MV_ )*R_*mD_/(1-kappa+R_*mV_)); 
    
    %D_star = K_ * ( 1 + (1 / (1+R_*mV_/(1-kappa)))*( (1/nubar)*log( 1 - kappa/kappa_new )*MV_ - mV_/(1-kappa) - 1 ) * RK_ ) / ( mD_ - betta*(1-rho) - (1-kappa)*(1/nubar)*log( 1 - kappa/kappa_new )*MD_ + ( (1/nubar)*log( 1 - kappa/kappa_new )*MV_ - mV_/(1-kappa) - 1 )*mD_/( 1 + R_*mV_/(1-kappa) ) );
    Vhat_ = ( RK_*K_ - R_ * mD_ * D_star )/ (1-kappa+R_*mV_);
    b = log( mV_ * Vhat_ + mD_ * D_star );
    spread = RK_ - R_;
end;


steady;
check(qz_zero_threshold=1e-9);

shocks;
    var epsA = 1; 
end;

stoch_simul( order = 2, irf = 0, periods = 10000 ) y inv h spread D_rate E_rate lev;


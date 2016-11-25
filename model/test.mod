var k a g logit_delta q c r pi mc disp psi X1 X2 Y inv I C K R Q spread H;

varexo epsA eps_psi epsM;

parameters gySS varrho alp zzeta betta deltaSS sigma_c rhodelta rhoG rhoA 
Ass Phi sigmaB xiB Theta kappaSS logit_deltaSS epsilon kappa_GK gam
sigma_g sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH xi gamR gamPi gamY Pi_bar pi_bar 
r_bar sigmaM Z_bar Disp_bar MC_bar C_bar H_bar Y_bar y_bar GSS;

load('../opts.mat');

gySS=0.2;
varrho=0.684;
alp=0.3;
zzeta=7.0;
betta=0.995;
deltabar=0.025;
logit_deltabar = log(deltabar/(1-deltabar));
sigma_c=2.0;
sigma_h=parameter_sigma_h;
psi_h=parameter_psi_h;
chi = .7;
epsilonC = parameter_habits_C; 
epsilonH = parameter_habits_H;
Ass=1;
rhoA=parameter_rhoA;
rhodelta=parameter_rhodelta;
rhoG=parameter_rhoG;
sigmaB=0.975;
xiB=0.00017;
epsilon = -2;
kappa_GK = 13;
gam=1e-8;
kappaSS=parameter_kappa;
Theta=parameter_Theta;
Phi=parameter_Phi;
sigma_g = parameter_sigma_g;
sigma_a = parameter_sigma_a;
sigma_psi = parameter_sigma_psi;
sigma_delta = parameter_sigma_delta;
logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));
Pi_bar = 1.005;
pi_bar = log(Pi_bar);
r_bar = log(Pi_bar/betta);
xi = 0.7;
sigmaM = .01;
gamR = parameter_gamR;
gamPi = parameter_gamPi;
gamY = parameter_gamY;
Z_bar = 1 / betta - ( 1 -  deltaSS );
Disp_bar = ( (1-xi)*( ( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)) )^(-zzeta) )/(1-xi*Pi_bar^zzeta);
MC_bar =  ( ( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta)) )*((zzeta-1)/zzeta)*(1-xi*betta*Pi_bar^zzeta)/(1- xi*betta*Pi_bar^(zzeta-1));
H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH))* ( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar ) / (MC_bar*(1-alp)) );
C_bar = ( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar )*( Ass * H_bar ) * ( MC_bar*alp/Z_bar ) ^ ( alp / ( 1 - alp ) );
Y_bar = ( Ass * H_bar ) * ( MC_bar*alp/Z_bar ) ^ ( alp / ( 1 - alp ) ) / Disp_bar;
y_bar = log(Y_bar);
GSS = gySS*Y_bar;

model;

# delta = 1/(1+exp(-logit_delta));
# lead_delta = 1/(1+exp(-logit_delta(+1)));
# lag_delta = 1/(1+exp(-logit_delta(-1)));
K = exp( k );
# lead_K = exp( k(+1) );
# lag_K = exp( k(-1) );
# A = exp( a );
# lag_C = exp( c(-1) );
C = exp( c );
# lead_C = exp( c(+1) );
I = exp( inv );
# lag_I = exp( inv(-1) );
# lead_I = exp( inv(+1) );
Q = exp( q );
# lead_Q = exp( q(+1) );
# lag_Q = exp( q(-1) );
# G = exp( g );
# lead_G = exp( g(+1) );
# Pi = exp( pi );
# lead_Pi = exp( pi(+1) );
# Disp = exp( disp );
# MC = exp( mc );
# lead_MC = exp( mc(+1) );
Y = C+I*(1-Phi*(1-I/lag_I)^2)+G;
# lead_Y = lead_C+lead_I*(1-Phi*(1-lead_I/I)^2) + lead_G;
# YW = Y*exp(disp);
# lead_YW = lead_Y*exp(disp(+1));
# Z = MC*alp*YW/(psi*lag_K);
# lead_Z = lead_MC*alp*lead_YW/(psi(+1)*K);
# lag_H = H(-1);
# lead_H = H(+1);
# UH = -varrho*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)-1));
# UC = (1-varrho)*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)-1))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)));
# lead_UC = (1-varrho)*((lead_C - epsilonC*C)^((1-varrho)*(1-sigma_c)-1))*((1-lead_H - epsilonH*(1-H))^(varrho*(1-sigma_c)));
# lead_Lambda = betta*lead_UC/UC; 
# W = MC*(1-alp)*YW/H;
R = exp( r );
# S = Q*K;
# lag_S = lag_Q*lag_K;
# RK = Pi*psi*(Z + (1-delta)*Q)/lag_Q;
# lead_RK = lead_Pi*psi(+1)*(lead_Z + (1-lead_delta)*lead_Q)/Q;
# Y_alt = (A*H)^(1-alp)*(psi*lag_K)^alp / Disp;
# lag_Disp = exp( disp(-1) );
# y = log(Y);
# E = 0;
# B = K*Q;
# N = 0;
# D = 0;
# D_rate = 0;
# E_rate = 0;
# lag_R = exp( r(-1) );
# kappa = 0;

1 = xi*Pi^(zzeta-1) + (1-xi)*(X1/X2)^(1-zzeta);
Disp =  lag_Disp*xi*Pi^zzeta + (1-xi)*(X1/X2)^(-zzeta);
X1 = (zzeta/(zzeta-1))*Y*MC + xi*X1(+1)*lead_Lambda*lead_Pi^zzeta;
X2 = Y + xi*X2(+1)*lead_Lambda*lead_Pi^(zzeta-1);
r = max( 0 , gamR*r(-1) + (1-gamR)*( r_bar + gamPi*(pi - pi_bar) + gamY*( y - y_bar) ) + sigmaM*epsM );
lead_Lambda*lead_RK/lead_Pi = 1;
spread = lead_RK - R;
I*(1-Phi*(1-I/lag_I)^2) = K - (1-delta)*psi*lag_K; 
1=Q*(1-Phi*(I/lag_I-1)^2 - (I/lag_I)*2*Phi*(I/lag_I-1))+lead_Lambda*2*Phi*(lead_I/I-1)*(lead_I/I)^2*lead_Q;
lead_Lambda*R/lead_Pi=1; 
log(Y) = (1-alp)*((a) + log(H)) + alp*( k(-1) + log(psi)) - disp; 
UH/UC = - W;

a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
g = log(GSS);
psi = exp(sigma_psi*eps_psi);
logit_delta = logit_deltabar;

end;

steady_state_model;
    psi = 1;
    pi = log( Pi_bar );
    q = 0;
    R_ = Pi_bar / betta;
    r = log( R_ );
    Z_ = 1 / betta - ( 1 -  deltaSS );
    X1_by_X2 = ( (1 - xi*Pi_bar^(zzeta-1))/(1-xi) )^(1/(1-zzeta));
    Disp_ = ( (1-xi)*X1_by_X2^(-zzeta) )/(1-xi*Pi_bar^zzeta);
    disp = log( Disp_ );
    MC_ =  X1_by_X2*((zzeta-1)/zzeta)*(1-xi*betta*Pi_bar^zzeta)/(1- xi*betta*Pi_bar^(zzeta-1));
    mc = log( MC_ );
    a = log( Ass );
    g = log( GSS );
    K_over_YW = MC_*alp/Z_;
    I_over_YW = deltaSS * K_over_YW;
    C_over_YW = 1/Disp_ - gySS/Disp_ - I_over_YW;
    H_ = H_bar;
    H = H_bar;
    YW_ = ( Ass * H_ ) * K_over_YW ^ ( alp / ( 1 - alp ) ) ;
    Y_ = YW_/Disp_;
    K_ = K_over_YW * YW_;
    k = log( K_ );
    inv = log(deltaSS) + k;
    c = log( C_over_YW * YW_ );
    logit_delta = logit_deltaSS;
    B_ = K_;
    X1 = (zzeta/(zzeta-1))*Y_*MC_ / (1 - xi*betta*Pi_bar^zzeta);
    X2 = Y_ / ( 1 - xi*betta*Pi_bar^(zzeta-1));
    r_PLUS_b = log( R_ * B_ );
    spread = (Z_+1-deltaSS)*Pi_bar-R_;

    Y = Y_;
    I = exp(inv);
    C = exp(c);
    R = R_;
    K = K_;
    Q = 1;
end;

steady;
check;

shocks;
    var epsA = 1; 
    var eps_psi = 1; 
    var epsM = 1;
end;


stoch_simul( pruning, replic = 500, order = 3, irf = 60 ) ;%Y_alt Y H I C K R B D E N D_rate E_rate spread Q S A delta G pi kappa; 


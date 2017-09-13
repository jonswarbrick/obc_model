@#include "sim_type.mod"
@#include "utility_type.mod"
@#include "shock_choice.mod"
@#include "adj_type.mod"

// Required
var k a logit_delta q c r pi mc disp psi;
// For analysis
var spread inv Y;
// with habits:
var H;

@#if utility_type == 5
var Xjr lambdaX;
@#endif

@#include "shocks_1.mod"

parameters varrho alp zzeta betta deltaSS sigma_c rhodelta rhoA
rho_psi Ass Phi sigmaB xiB Theta kappaSS logit_deltaSS epsilon kappa_GK gam
sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr gy;

load('../opts.mat');

varrho=2.6;
alp=0.3;
zzeta=7.0;
betta=0.995;
deltabar=0.025;
logit_deltabar = log(deltabar/(1-deltabar));
sigma_c=parameter_sigma_c;
gam_jr = parameter_gam_jr;
theta_jr = parameter_theta_jr;
sigma_h=parameter_sigma_h;
psi_h=parameter_psi_h;
chi = .7;
epsilonC = parameter_habits_C;
epsilonH = parameter_habits_H;
Ass=1;
rhoA=parameter_rhoA;
rhodelta=parameter_rhodelta;
rho_psi = parameter_rho_psi;
gy = 0;

sigmaB=0.975;
xiB=0.003;
epsilon = -2;
kappa_GK = 13;

gam=1e-8;
kappaSS=parameter_kappa;

Theta=parameter_Theta;
Phi=parameter_Phi;
sigma_a = parameter_sigma_a;
sigma_psi = parameter_sigma_psi;
sigma_delta = parameter_sigma_delta;

%logit_deltaSS = logit_deltabar -  sigma_delta/(1-rhodelta);
logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));

parameters C_bar H_bar;

@#if utility_type == 1
    H_bar = 1 / ( 1 + varrho/(1-varrho) *((1-epsilonC)/(1-epsilonH)) * (1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))) / (1-alp) );
@#endif
@#if utility_type == 2
    H_bar  = call_csolve_rbc;
@#endif
@#if utility_type == 3
    H_bar  =  ( 1 / ( (1 - epsilonH)^(psi_h) * (1 - epsilonC)*(1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))) / (1-alp) ) )^(1/(psi_h+1));
@#endif
@#if utility_type == 4
    H_bar = 1 / ( 1 + varrho/(1-varrho) * (1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))) / (1-alp) );
@#endif
@#if utility_type == 5
    H_bar = call_csolve_jr_rbc;
    %H_bar = (((1-alp)*((alp/(1 / betta - ( 1 -  deltaSS )))^(alp/(1-alp)))^(1-gam_jr) ) / ( varrho*( (1-alp)*gam_jr*((1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))))^(gam_jr-1) + (theta_jr+1-gam_jr)*((1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))))^gam_jr) ) )^(1/theta_jr);
@#endif

C_bar = (  (1 - deltaSS * alp/(1 / betta - ( 1 -  deltaSS ))) ) * ( ( Ass * H_bar ) * (alp/(1 / betta - ( 1 -  deltaSS ))) ^ ( alp / ( 1 - alp ) ) );

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
# Pi = exp( pi );
# lead_Pi = exp( pi(+1) );
# Disp = exp( disp );
# MC = exp( mc );
# lead_MC = exp( mc(+1) );

@#if adj_type == 1
% CEE
%Y = (C+I*(1-Phi*(1-I/lag_I)^2))/(1-gy);
%# lead_Y = (lead_C+lead_I*(1-Phi*(1-lead_I/I)^2))/(1-gy);
Y = (C+I)/(1-gy);
# lead_Y = (lead_C+lead_I)/(1-gy);
@#endif
@#if adj_type == 2
% Ireland (2003) costs
Y = (C+I-Phi*psi*lag_K*(K/(psi*lag_K)-1)^2)/(1-gy);
# lead_Y = (lead_C + lead_I - Phi*psi(+1)*K*(lead_K/(psi(+1)*K)-1)^2)/(1-gy);
@#endif

# YW = Y*exp(disp);
# lead_YW = lead_Y*exp(disp(+1));
# Z = MC*alp*YW/(psi*lag_K);
# lead_Z = lead_MC*alp*lead_YW/(psi(+1)*K);

@#if utility_type == 1
    # lag_H = H(-1);
    # lead_H = H(+1);
    # UH = -varrho*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)-1));
    # lambdaC = (1-varrho)*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)-1))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)));
    # lead_lambdaC = (1-varrho)*((lead_C - epsilonC*C)^((1-varrho)*(1-sigma_c)-1))*((1-lead_H - epsilonH*(1-H))^(varrho*(1-sigma_c)));
@#endif
@#if utility_type == 2
    # lag_H = H(-1);
    # UH = -varrho*(1-H - epsilonH*(1-lag_H))^(-sigma_h);
    # lambdaC = (C - epsilonC*lag_C)^(-1);
    # lead_lambdaC = (lead_C - epsilonC*C)^(-1);
@#endif
@#if utility_type == 3
    # lag_H = H(-1);
    # UH = -(H - epsilonH*lag_H)^(psi_h);
    # lambdaC = (C - epsilonC*lag_C)^(-1);
    # lead_lambdaC = (lead_C - epsilonC*C)^(-1);
@#endif
@#if utility_type == 4
    # lag_H = H(-1);
    # lead_H = H(+1);
    # UH = -varrho*((1-H)/C)^(varrho-1)*( C^(1-varrho)*(1-H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
    # lambdaC = (1-varrho)*((1-H)/C)^varrho*( C^(1-varrho)*(1-H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
    # lead_lambdaC = (1-varrho)*((1-lead_H)/lead_C)^varrho*( lead_C^(1-varrho)*(1-lead_H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
@#endif
@#if utility_type == 5
    # lead_H = H(+1);
    # UH = - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * theta_jr*varrho*Xjr*H^(theta_jr-1);
    # UX =  - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
    # UC = (C - varrho*H^theta_jr*Xjr)^(-sigma_c);
    # lead_UC = (lead_C - varrho*lead_H^theta_jr*Xjr(+1))^(-sigma_c);
    %   # lambdaC = UC + lambdaX*gam_jr*C^(gam_jr-1);
    %   # lead_lambdaC = lead_UC + lambdaX(+1)*gam_jr*lead_C^(gam_jr-1);
    # lambdaC = UC + lambdaX*gam_jr*Xjr/C;
    # lead_lambdaC = lead_UC + lambdaX(+1)*gam_jr*Xjr(+1)/lead_C;
    Xjr = C^gam_jr * Xjr(-1)^(1-gam_jr);
    lambdaX = UX + betta * (1-gam_jr) * lambdaX(+1) * Xjr(+1) / Xjr;
@#endif

# lead_Lambda = betta*lead_lambdaC/lambdaC;
# W = MC*(1-alp)*YW/H;
# R = exp( r );
# S = Q*K;
# lag_S = lag_Q*lag_K;
# RK = Pi*psi*(Z + (1-delta)*Q)/lag_Q;
# lead_RK = lead_Pi*psi(+1)*(lead_Z + (1-lead_delta)*lead_Q)/Q;
# Y_alt = (A*H)^(1-alp)*(psi*lag_K)^alp / Disp;


# E = 0;
# B = K*Q;
# N = 0;
# D = 0;
# D_rate = 0;
# E_rate = 0;
# kappa = 0;

lead_Lambda*R = lead_Lambda*lead_RK;

pi = 0;
disp = 0;
mc = 0;



spread = lead_RK - R;

@#if adj_type == 1
% CEE
I*(1-Phi*(1-I/lag_I)^2) = K - (1-delta)*psi*lag_K;
1=Q*(1-Phi*(I/lag_I-1)^2 - (I/lag_I)*2*Phi*(I/lag_I-1))+lead_Lambda*2*Phi*(lead_I/I-1)*(lead_I/I)^2*lead_Q;
@#endif
@#if adj_type == 2
%% Ireland (2003) adjustment costs
I = K - (1-delta)*psi*lag_K;
1=Q-2*Phi*(K/(psi*lag_K)-1) - lead_Lambda*( (lead_Q-1)*(1-lead_delta) + Phi*(lead_K/(psi(+1)*K)-1)^2 - 2*Phi*(lead_K/(psi(+1)*K)-1)*(lead_K/(psi(+1)*K))  ); //q
@#endif


lead_Lambda*R/lead_Pi=1;
log(Y) = (1-alp)*((a) + log(H)) + alp*( k(-1) + log(psi)) - disp;

UH/lambdaC = - W;


@#include "shock_processes.mod"

end;

steady_state_model;
    mc = 0;
    disp = 0;
    psi = 1;
    q = 0;
    R_ = 1 / betta;
    RK_ = R_;
    r = log( R_ );
    Z_ = 1 / betta - ( 1 -  deltaSS ) ;
    a = log( Ass );
    K_over_Y = alp/Z_;
    I_over_Y = deltaSS * K_over_Y;
    C_over_Y = 1 - gy - I_over_Y;
    Y_over_H = Ass * K_over_Y ^ ( alp / ( 1 - alp ) );
    W_ = (1-alp)*Y_over_H;
    H = H_bar;
    H_ = H;
    Y_ = ( Ass * H_ ) * K_over_Y ^ ( alp / ( 1 - alp ) );
    K_ = K_over_Y * Y_;
    k = log( K_ );
    inv = log(deltaSS) + k;
    C_ = C_over_Y * Y_;
    c = log( C_ );
    Xjr = C_;
    UX_ =  - (C_ - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
    lambdaX = UX_ / ( 1 - betta*( (1-gam_jr) ) );
    logit_delta = logit_deltaSS;
    B_ = K_;
    r_PLUS_b = log( R_ * B_ );

    spread = 0;
    Y = Y_;
    //I = exp(inv);
    //C = exp(c);
    //G = GSS;
    //R = R_;
    //K = K_;
    //G = GSS;
    //B = K_;
    //Q = 1;
    //D = 0;
    //E = 0;
end;


steady;
check;

shocks;
@#include "shocks_2.mod"
end;


@#if sim_type == 1
    stoch_simul( order = 2, irf = 0, periods = 1000 );
@#endif
@#if sim_type == 2
    stoch_simul( nograph , replic = 256, order = 3, irf = 60, periods = 0 , irf_shocks = ( @#include "shocks_3.mod" ) );
@#endif

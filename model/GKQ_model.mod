@#include "which_model.mod"
@#include "sim_type.mod"
@#include "utility_type.mod"
@#include "shock_choice.mod"
@#include "adj_type.mod"

// Required
var k a g logit_delta q c r pi mc disp psi;
// For analysis
var spread inv Y;
// with habits:
var H;

@#if shock_choice == 1
    varexo epsA eps_psi;
@#endif
@#if shock_choice == 2
    varexo epsA epsdelta;
@#endif
@#if shock_choice == 4
    varexo epsA;
@#endif

parameters gySS varrho alp zzeta betta deltaSS sigma_c rhodelta rhoG rhoA 
Ass Phi sigmaB xiB Theta kappaSS logit_deltaSS epsilon kappa_GK gam
sigma_g sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr ;

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
gam_jr = 0.001;
theta_jr = 1.5;
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

%logit_deltaSS = logit_deltabar -  sigma_delta/(1-rhodelta);
logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));

//required
var mue mus nub E Thetax QE;  
//for analysis
var D_rate E_rate;

parameters C_bar K_by_Y H_bar GSS;

H_bar  = call_Hbar_gkq;
K_by_Y  = call_KbyY_gkq;

C_bar = ( 1 - gySS - deltaSS * K_by_Y )*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );
GSS = gySS*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );

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
# G = exp( g );
# lead_G = exp( g(+1) );
# Pi = exp( pi );
# lead_Pi = exp( pi(+1) );
# Disp = exp( disp );
# MC = exp( mc );
# lead_MC = exp( mc(+1) );

@#if adj_type == 1
% CEE
Y = C+I*(1-Phi*(1-I/lag_I)^2)+G;
# lead_Y = lead_C+lead_I*(1-Phi*(1-lead_I/I)^2) + lead_G;
@#endif
@#if adj_type == 2
% Ireland (2003) costs
Y = C+I+G-Phi*psi*lag_K*(K/(psi*lag_K)-1)^2;
# lead_Y = lead_C + lead_I + lead_G - Phi*psi(+1)*K*(lead_K/(psi(+1)*K)-1)^2;
@#endif

# YW = Y*exp(disp);
# lead_YW = lead_Y*exp(disp(+1));
# Z = MC*alp*YW/(psi*lag_K);
# lead_Z = lead_MC*alp*lead_YW/(psi(+1)*K);

@#if utility_type == 1
    # lag_H = H(-1);
    # lead_H = H(+1);
    # UH = -varrho*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)-1));
    # UC = (1-varrho)*( (C - epsilonC*lag_C)^((1-varrho)*(1-sigma_c)-1))*((1-H - epsilonH*(1-lag_H))^(varrho*(1-sigma_c)));
    # lead_UC = (1-varrho)*((lead_C - epsilonC*C)^((1-varrho)*(1-sigma_c)-1))*((1-lead_H - epsilonH*(1-H))^(varrho*(1-sigma_c)));
@#endif
@#if utility_type == 2
    # lag_H = H(-1);
    # UH = -varrho*(1-H - epsilonH*(1-lag_H))^(-sigma_h);
    # UC = (C - epsilonC*lag_C)^(-1);
    # lead_UC = (lead_C - epsilonC*C)^(-1);
@#endif
@#if utility_type == 3
    # lag_H = H(-1);
    # UH = -(H - epsilonH*lag_H)^(psi_h);
    # UC = (C - epsilonC*lag_C)^(-1);
    # lead_UC = (lead_C - epsilonC*C)^(-1);
@#endif
@#if utility_type == 4
    # lag_H = H(-1);
    # lead_H = H(+1);
    # UH = -varrho*((1-H)/C)^(varrho-1)*( C^(1-varrho)*(1-H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
    # UC = (1-varrho)*((1-H)/C)^varrho*( C^(1-varrho)*(1-H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
    # lead_UC = (1-varrho)*((1-lead_H)/lead_C)^varrho*( lead_C^(1-varrho)*(1-lead_H)^varrho - epsilonC*C_bar^(1-varrho)*(1-H_bar)^varrho )^(-sigma_c);
@#endif
@#if utility_type == 5
    # lag_H = H(-1);
    # lead_H = H(+1);   
    # habits = epsilonC*( lag_C - varrho*lag_H^theta_jr*(lag_C^gam_jr*lag_H^(1-gam_jr)));
    # lead_habits = epsilonC*(C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)));
    # UH = - (theta_jr+1-gam_jr)*varrho*C^gam_jr*H^(theta_jr-gam_jr)*(C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)) - habits )^(-sigma_c);
    # UC = (1 - gam_jr*varrho*H^theta_jr*(H/C)^(1-gam_jr))*(C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)) - habits )^(-sigma_c);
    # lead_UC = (1 - gam_jr*varrho*lead_H^theta_jr*(lead_H/lead_C)^(1-gam_jr))*(lead_C - varrho*lead_H^theta_jr*(lead_C^gam_jr*lead_H^(1-gam_jr)) - lead_habits )^(-sigma_c);
@#endif

# lead_Lambda = betta*lead_UC/UC; 
# W = MC*(1-alp)*YW/H;
# R = exp( r );
# S = Q*K;
# lag_S = lag_Q*lag_K;
# RK = Pi*psi*(Z + (1-delta)*Q)/lag_Q;
# lead_RK = lead_Pi*psi(+1)*(lead_Z + (1-lead_delta)*lead_Q)/Q;
# Y_alt = (A*H)^(1-alp)*(psi*lag_K)^alp / Disp;

# lag_R = exp( r(-1) );
# lead_R = exp( r(+1) );
# xE = E*QE/S;		
# lag_xE = E(-1)*QE(-1)/S;	
# DThetax = Theta*(epsilon + kappa_GK*xE);
# phi_GK = nub/(Thetax-(mus+mue*xE));
# lag_phi_GK = nub(-1)/(Thetax(-1) - (mus(-1) + mue(-1)*lag_xE));
# N = S/phi_GK;
# lag_N = lag_S/lag_phi_GK;	
# B = S - N - QE*E;
# lag_B = lag_S - lag_N - QE(-1)*E(-1);	
# Omega = 1 - sigmaB + sigmaB*Thetax*phi_GK;	
# lead_Omega = 1 - sigmaB + sigmaB*Thetax(+1)*phi_GK;	
# RE =(Z + (1 - delta)*QE)/QE(-1);	
# lead_RE =(lead_Z + (1 - lead_delta)*QE(+1))/QE;
# D = (1-sigmaB)*(RK*lag_S-lag_R*lag_B - RE*QE(-1)*E(-1));
D_rate = D/N;
E_rate = E/N; 

N = (sigmaB + xiB)*RK*lag_S - lag_R*sigmaB*lag_B - RE*sigmaB*QE(-1)*E(-1);
mus = lead_Lambda*lead_Omega*spread;	
nub = lead_Lambda*lead_Omega*R;
mue = lead_Lambda*lead_Omega*(R-lead_RE);
Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);	
lead_Lambda*lead_RE = 1;
(mus+mue*xE)*DThetax = mue * Thetax;


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

UH/UC = - W;


a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
%g-log(GSS) = rhoG*(g(-1)-log(GSS))-sigma_g*epsG;
g = log(GSS);

@#if shock_choice == 1
    psi = exp(sigma_psi*eps_psi);
    logit_delta = logit_deltabar;
@#endif
@#if shock_choice == 2
    logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*epsdelta;
    %logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*epsdelta^3;
    %logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*(epsdelta^2-1);
    psi = 1;
@#endif
@#if shock_choice == 4
    logit_delta = logit_deltabar;
    psi = 1;
@#endif



end;


steady;
check;

steady_state_model;
disp = 0;
pi = 0;
mc = 0;
q = 0;
R_ = 1 / betta;
r = log( R_ );
a = log( Ass );
logit_delta = logit_deltaSS;
psi = 1;

Z_ = alp/K_by_Y;
I_by_Y = deltaSS * K_by_Y;
C_by_Y = 1 - gySS - I_by_Y;
H = H_bar;
Y = ( Ass * H ) * K_by_Y ^ ( alp / ( 1 - alp ) );
g = log ( GSS );
K_ = K_by_Y * Y;
k = log( K_ );
c = log( C_by_Y * Y );
inv = log( I_by_Y * Y );
I_ = I_by_Y * Y;
C_ = C_by_Y * Y;

RK_ = (Z_ + (1 - deltaSS));
RE_ = R_;
QE =Z_/(RE_  - (1 - deltaSS));
S_ = K_;
mue = 0;
spread = RK_ - R_;
xE_ = -epsilon/kappa_GK;
Thetax = Theta*(1 + epsilon*xE_ + kappa_GK*xE_^2/2);
E = xE_*S_/QE;
N_ = ((sigmaB + xiB)*RK_ - R_*sigmaB)*S_/(1-R_*sigmaB);
phi_GK_ = S_/N_;
Omega_ = 1 - sigmaB + sigmaB*Thetax*phi_GK_;	
mus = betta*Omega_*spread;	
nub = Omega_;
D_ = (1-sigmaB)*(RK_*S_-R_*(S_-N_-QE*E) - RE_*QE*E); 
D_rate = D_/N_;
E_rate = E/N_;
end;



shocks;
@#include "shocks_2.mod"
end;


@#if sim_type == 1
    stoch_simul( order = 2, irf = 0, periods = 1000 ) Y I C K spread Q D_rate E_rate delta G A kappa;
@#endif
@#if sim_type == 2
    stoch_simul( nograph , replic = 500, order = 3, irf = 60, periods = 0 , irf_shocks = ( eps_psi ) ) Y H inv c k r D_rate E_rate spread q S A delta G pi; 
@#endif



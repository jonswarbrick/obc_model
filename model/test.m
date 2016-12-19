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

sigmaB=0.975;
xiB=0.00017;
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


H_bar  = call_Hbar_gk
K_by_Y  = call_KbyY_gk;

C_bar = ( 1 - deltaSS * K_by_Y )*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );
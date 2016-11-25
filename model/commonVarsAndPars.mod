// Required
var k a g logit_delta q c r pi mc disp psi;
// For analysis
var spread inv Y;
// with habits:
var H;

@#include "shocks_1.mod"

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
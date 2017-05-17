clear;
% Deterministic steady state of obc model


cd('model');
parameter_sigma_a = 0.003500638344089;
parameter_rhoA = 0.828157065533977;
parameter_Theta = 0.1;% 0.569692624538134;
parameter_sigma_psi = 0;
parameter_rho_psi = 0;
parameter_Phi = 2;
parameter_kappa = 0.05;
parameter_kappa_new = .1;
parameter_nubar = 400;
parameter_sigma_c =2;
parameter_gam_jr = 0.001;
parameter_theta_jr = 1.4;
parameter_habits_C = 0;
parameter_habits_H = 0;
parameter_sigma_g = 0;
parameter_sigma_delta = .9;
parameter_rhoG = 0.95;
parameter_rhodelta = 0.85;
parameter_gamR = 0.9;
parameter_gamPi = 2;
parameter_gamY = 0.4;
parameter_sigma_h = 2.37;
parameter_psi_h = 1;

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
logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));
nubar = parameter_nubar;
kappa_new = parameter_kappa_new;
lambdaB_bar = gam*(1-(1-gam)*(1-Theta));

SZ_bar = zeros(1,7);
SZ = SZ_bar;
MD_bar = zeros(1,7);
prodR = zeros(1,7);
MD = zeros(1,7);
mD = zeros(1,7);
lagD = zeros(1,7);

SZ_bar(7) = lambdaB_bar;
SZ_bar(6) = lambdaB_bar + (1-gam)*SZ_bar(7)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ_bar(5) = lambdaB_bar + (1-gam)*SZ_bar(6)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ_bar(4) = lambdaB_bar + (1-gam)*SZ_bar(5)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ_bar(3) = lambdaB_bar + (1-gam)*SZ_bar(4)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ_bar(2) = lambdaB_bar + (1-gam)*SZ_bar(3)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
SZ_bar(1) = lambdaB_bar + (1-gam)*SZ_bar(2)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_bar );
for lag=1:7
    MD_bar(lag) = SZ_bar(lag)*(1 / betta)^(lag)*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_bar );
end

H_bar = call_csolve_jr_obc;
C_bar = ( 1 - deltaSS * alp/( ( 1/((1-gam)*betta*( ( 1 - gam*(1-( (1-gam)*betta*MD_bar(1)/(1+(1-gam)*betta*MD_bar(1)) ))*(1-(1-gam)*(1-Theta)) )/(1-gam) )) )-1+deltaSS ) ) * ( Ass * H_bar )*(alp/(1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) )))^( alp / ( 1 - alp ) );

mc = 0;
disp = 0;
pi = 0;
psi = 1;
logit_delta = logit_deltaSS;
nu = nubar;
q = 0;
R_ = 1 / betta;
r = log( R_ );
a = log( Ass );
lambdaD_ = 0;
lambdaB_ = gam*(1-(1-gam)*(1-Theta));
SZ(7) = lambdaB_;
SZ(6) = lambdaB_ + (1-gam)*SZ(7)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
SZ(5) = lambdaB_ + (1-gam)*SZ(6)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
SZ(4) = lambdaB_ + (1-gam)*SZ(5)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
SZ(3) = lambdaB_ + (1-gam)*SZ(4)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
SZ(2) = lambdaB_ + (1-gam)*SZ(3)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
SZ(1) = lambdaB_ + (1-gam)*SZ(2)*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
for lag=1:7
    prodR(lag) = R_^(lag);
    MD(lag) = SZ(lag)*R_^(lag)*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_ );
    mD(lag) = ( MD(lag) + (1-gam)*(1-Theta)*R_^(lag) )/( 1 - (1-gam)*(1-Theta) );
end
kappa = (1-gam)*betta*MD(1)/(1+(1-gam)*betta*MD(1));
MV_ = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam);
mv = log( MV_ );
mV_ = MV_ / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
RK_ = 1/((1-gam)*betta*MV_);
Z_ = RK_-1+deltaSS;
K_over_Y = alp/Z_;
I_over_Y = deltaSS * K_over_Y;
C_over_Y = 1 - I_over_Y;
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
sum_mD = sum(mD);
sum_MD = sum(MD);
AUXBD = sum_mD/(1+mV_*R_/(1-kappa));
AUXBS = (mV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa));
AUXED = (1/nubar)*( log( 1 - kappa/kappa_new ) )*( ((MV_*R_/(1-kappa))/(1+mV_*R_/(1-kappa)))*sum_mD - sum_MD );
AUXES = (1/nubar)*( log( 1 - kappa/kappa_new ) )*((MV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa)));
D = K_*(RK_ - 1 - (R_-1)*AUXBS - (1-kappa)*AUXES)/(1-(1-kappa)*AUXED + (R_-1)*AUXBD);
E_ = AUXED*D - AUXES*K_;
B_ = AUXBD*D + AUXBS*K_;
b = log( B_ );
for lag=1:7
    lagD(lag) = D;
end
spread = RK_ - R_;

Y = Y_;
E_rate = E_/(K_-B_);
D_rate = D/(K_-B_);
N = K_-B_;
cap_ass = N/K_
liab_ass = B_/K_
cd('../');

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


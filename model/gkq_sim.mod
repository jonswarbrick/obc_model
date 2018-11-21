var k a nub mus mue Thetax QE E q logit_delta c r inv spread psi h Xjr lambdaX;

varexo epsA;

parameters varrho alp zzeta betta deltaSS sigma_c rhodelta rhoA Ass Phi
rho_psi sigmaB xiB Theta kappaSS logit_deltaSS  epsilon kappa_GK gam
sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr ;
parameters C_bar K_by_Y H_bar;

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
xiB = 0.003;
sigmaB = 0.982;
xiB=0.001;
epsilon = -1.21;
kappa_GK = 13.41;

gam=1e-8;
kappaSS=.05;

Theta=0.893717341213108;
Phi=2;
sigma_a = 0.005710743233638;
sigma_psi = 0;
sigma_delta = .9;

logit_deltaSS = logit_deltabar;
deltaSS = 1/(1+exp(-logit_deltaSS));

H_bar  = solve_H_gkq;
K_by_Y  = solve_KbyY_gkq;
C_bar = ( 1 - deltaSS * K_by_Y )*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );


model;
	# delta = 1/(1+exp(-logit_delta));
	# lead_delta = 1/(1+exp(-logit_delta(+1)));
	# lag_delta = 1/(1+exp(-logit_delta(-1)));
	# K = exp( k );
	# lag_K = exp( k(-1) );
	# lead_K = exp( k(+1) );
	# A = exp( a );
	# lag_A = exp( a(-1) );
	# I = exp( inv );
	# lead_I = exp( inv(+1) );
	# lag_I = exp( inv(-1) );
	# C = exp( c );
	# lead_C = exp( c(+1) );
	# lag_C = exp( c(-1) );
    # Q = exp( q );
    # lead_Q = exp( q(+1) );
    # lag_Q = exp( q(-1) );
    # Y = C+I;
    # lead_Y = lead_C+lead_I;
    # y = log( Y );
    # H = exp( h );
    # lead_H = exp( h(+1) );
    # UH = - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * theta_jr*varrho*Xjr*H^(theta_jr-1);
    # UX =  - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
    # UC = (C - varrho*H^theta_jr*Xjr)^(-sigma_c);
    # lead_UC = (lead_C - varrho*lead_H^theta_jr*Xjr(+1))^(-sigma_c);
    # lambdaC = UC + lambdaX*gam_jr*Xjr/C;
    # lead_lambdaC = lead_UC + lambdaX(+1)*gam_jr*Xjr(+1)/lead_C;
	# lead_Lambda = betta*lead_lambdaC/lambdaC;
    # Z = alp*Y/(psi*lag_K);
    # lead_Z = alp*lead_Y/(psi(+1)*K);
	# w = log(1-alp) + y - h;
	# R = exp( r );
	# lag_R = exp( r(-1) );
    # S = Q*K;
    # lag_S = lag_Q*lag_K;
    # lead_S = lead_Q*lead_K;
    # RK = psi*(Z + (1-delta)*Q)/lag_Q;
    # lead_RK = psi(+1)*(lead_Z + (1-lead_delta)*lead_Q)/Q;
    # xE = E*QE/S;
    # lag_xE = E(-1)*QE(-1)/lag_S;
    # lead_xE = E(+1)*QE(+1)/lead_S;
    # DThetax = Theta*(epsilon + kappa_GK*xE);
    # phi_GK = nub/(Thetax-(mus+mue*xE));
    # lag_phi_GK = nub(-1)/(Thetax(-1)-(mus(-1)+mue(-1)*lag_xE));
    # lead_phi_GK = nub(+1)/(Thetax(+1)-(mus(+1)+mue(+1)*lead_xE));
    # N = S/phi_GK;
    # lag_N = lag_S/lag_phi_GK;
    # B = S - N - QE*E;
    # lag_B = lag_S - lag_N - QE(-1)*E(-1);
    # Omega = 1 - sigmaB + sigmaB*Thetax*phi_GK;
    # lead_Omega = 1 - sigmaB + sigmaB*Thetax(+1)*lead_phi_GK;
    # RE =(Z + (1 - delta)*QE)/QE(-1);
    # lead_RE =(lead_Z + (1 - lead_delta)*QE(+1))/QE;
    # D = (1-sigmaB)*(RK*lag_S-lag_R*lag_B - RE*QE(-1)*E(-1));
    # D_rate = D/N;
    # E_rate = (E-E(-1))/N;
    # lev = k - log(N);

% Model equations

    I*(1-Phi*(1-I/lag_I)^2) = K - (1-delta)*psi*lag_K;
    1=Q*(1-Phi*(I/lag_I-1)^2 - (I/lag_I)*2*Phi*(I/lag_I-1))+lead_Lambda*2*Phi*(lead_I/I-1)*(lead_I/I)^2*lead_Q;
 
    lead_Lambda*R=1;
    y = (1-alp)*((a) + h) + alp*( k(-1) + log(psi));
    Xjr = C^gam_jr * Xjr(-1)^(1-gam_jr);
    lambdaX = UX + betta * (1-gam_jr) * lambdaX(+1) * Xjr(+1) / Xjr;
    log( -UH/lambdaC ) = w;

    N = (sigmaB + xiB)*RK*lag_S - lag_R*sigmaB*lag_B - RE*sigmaB*QE(-1)*E(-1);
    mus = lead_Lambda*lead_Omega*spread;
    nub = lead_Lambda*lead_Omega*R;
    mue = lead_Lambda*lead_Omega*(R-lead_RE);
    Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);
    (mus+mue*xE)*DThetax = mue * Thetax;
    lead_Lambda*lead_RE = 1;
    spread = lead_RK - R;

    a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
    logit_delta = logit_deltabar;
    psi = 1;

end;

steady_state_model;
    q = 0;
    R_ = 1 / betta;
    r = log( R_ );
    a = log( Ass );
    logit_delta = logit_deltaSS;
    psi = 1;

    Z_ = alp/K_by_Y;
    I_by_Y = deltaSS * K_by_Y;
    C_by_Y = 1 - I_by_Y;
    h = log( H_bar );
    H_ = H_bar;
    Y_ = ( Ass * H_ ) * K_by_Y ^ ( alp / ( 1 - alp ) );
    K_ = K_by_Y * Y_;
    k = log( K_ );
    C_ = C_by_Y * Y_;
    c = log( C_ );
    Xjr = C_;
    UX_ =  - (C_ - varrho*H_^theta_jr*Xjr)^(-sigma_c) * varrho * H_^(theta_jr);
    lambdaX = UX_ / ( 1 - betta*( (1-gam_jr) ) );
    inv = log( I_by_Y * Y_ );
    I_ = I_by_Y * Y_;

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

end;

steady;
check;

shocks;
    var epsA = 1; 
end;

stoch_simul( order = 2, irf = 0, periods = 10000 ) y inv h spread D_rate E_rate lev;


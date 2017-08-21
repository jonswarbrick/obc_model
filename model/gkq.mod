@#include "which_model.mod"
@#include "sim_type.mod"
@#include "utility_type.mod"
@#include "shock_choice.mod"
@#include "adj_type.mod"

var k a nub mus mue Thetax QE E q logit_delta c r inv spread psi Y H D_rate E_rate;

@#if shock_choice == 1
    varexo epsA eps_psi;
@#endif
@#if shock_choice == 2
    varexo epsA epsdelta;
@#endif
@#if shock_choice == 4
    varexo epsA epsG;
@#endif


@#if utility_type == 5
var Xjr lambdaX;
@#endif

parameters varrho alp zzeta betta deltaSS sigma_c rhodelta rhoA Ass Phi
rho_psi sigmaB xiB Theta kappaSS logit_deltaSS  epsilon kappa_GK gam
sigma_a sigma_psi sigma_delta deltabar logit_deltabar
chi sigma_h psi_h epsilonC epsilonH gam_jr theta_jr ;
parameters C_bar K_by_Y H_bar;

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
xiB=0.0045;
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


H_bar  = call_Hbar_gkq;
K_by_Y  = call_KbyY_gkq;

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
@#if adj_type == 1
% CEE
Y = C+I*(1-Phi*(1-I/lag_I)^2);
# lead_Y = lead_C+lead_I*(1-Phi*(1-lead_I/I)^2);
@#endif
@#if adj_type == 2
% Ireland (2003) costs
Y = C+I-Phi*psi*lag_K*(K/(psi*lag_K)-1)^2;
# lead_Y = lead_C + lead_I - Phi*psi(+1)*K*(lead_K/(psi(+1)*K)-1)^2;
@#endif


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
    # Z = alp*Y/(psi*lag_K);
    # lead_Z = alp*lead_Y/(psi(+1)*K);
	# W = (1-alp)*Y/H;
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


% Model equations

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

    lead_Lambda*R=1;
    log(Y) = (1-alp)*((a) + log(H)) + alp*( k(-1) + log(psi));
    UH/lambdaC = - W;

    N = (sigmaB + xiB)*RK*lag_S - lag_R*sigmaB*lag_B - RE*sigmaB*QE(-1)*E(-1);
    mus = lead_Lambda*lead_Omega*spread;
    nub = lead_Lambda*lead_Omega*R;
    mue = lead_Lambda*lead_Omega*(R-lead_RE);
    Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);
    (mus+mue*xE)*DThetax = mue * Thetax;
    lead_Lambda*lead_RE = 1;
    D_rate = D/N;
    E_rate = (E-E(-1))/N;
    spread = lead_RK - R;

%% Shock Process
@#include "shock_processes.mod"

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
    H = H_bar;
    Y = ( Ass * H ) * K_by_Y ^ ( alp / ( 1 - alp ) );
    K_ = K_by_Y * Y;
    k = log( K_ );
    C_ = C_by_Y * Y;
    c = log( C_ );
    Xjr = C_;
    UX_ =  - (C_ - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
    lambdaX = UX_ / ( 1 - betta*( (1-gam_jr) ) );
    inv = log( I_by_Y * Y );
    I_ = I_by_Y * Y;

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
    E_rate = 0;

end;

steady;
check;

shocks;
@#if shock_choice == 1
    var epsA = 1;
    %var epsdelta = 1;
    %var epsG = 1;
    var eps_psi = 1;
    %var epsM = 1;
@#endif
@#if shock_choice == 2
    var epsA = 1;
    var epsdelta = 1;
    %var epsG = 1;
    %var eps_psi = 1;
    %var epsM = 1;
@#endif
@#if shock_choice == 4
    var epsA = 1;
    %var epsdelta = 1;
    var epsG = 1;
    %var eps_psi = 1;
    %var epsM = 1;
@#endif
end;


@#if sim_type == 1
    stoch_simul( order = 2, irf = 0, periods = 1000 );
@#endif
@#if sim_type == 2
    stoch_simul( nograph , replic = 500, order = 3, irf = 60, periods = 0 , irf_shocks = ( eps_psi ) );
@#endif

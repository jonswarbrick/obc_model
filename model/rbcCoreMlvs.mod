
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

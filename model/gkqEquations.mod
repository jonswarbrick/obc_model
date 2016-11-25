# lag_R = exp( r(-1) );

@#if model_type == 2  
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

    N = (sigmaB + xiB)*RK*lag_S - lag_R*sigmaB*lag_B - RE*sigmaB*QE(-1)*E(-1);
        mus = lead_Lambda*lead_Omega*spread;	
        nub = lead_Lambda*lead_Omega*R;
        mue = lead_Lambda*lead_Omega*(R-lead_RE);
        Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);	
        QE lead_Lambda*lead_RE = 1;
    (mus+mue*xE)*DThetax = mue * Thetax;
    E_rate = E/N;
@#endif


@#if model_type == GK
        # phi_GK = m/Theta;
        # lag_phi_GK = m(-1)/Theta;
        # N = S/phi_GK;
        # lag_N = lag_S/lag_phi_GK;
        # B = S - N;
        # lag_B = lag_S - lag_N;
        # Omega = 1 - sigmaB + sigmaB*m;
        # lead_Omega = 1 - sigmaB + sigmaB*m(+1);
        # lag_R = exp( r(-1) );

        N = (sigmaB + xiB)*RK*lag_S-lag_R*sigmaB*lag_B;
        m*(1-lead_Lambda*lead_Omega*(lead_RK-R)/Theta)=lead_Lambda*lead_Omega*R;
@#endif

D_rate = D/N;
pi = 0;
disp = 0;
mc = 0;

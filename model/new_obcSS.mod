
steady_state_model;
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
    g = log( GSS );
    lambdaD_ = 0;
    lambdaB_ = gam*(1-(1-gam)*(1-Theta));
    SZ7 = lambdaB_;
    SZ6 = lambdaB_ + (1-gam)*SZ7*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ5 = lambdaB_ + (1-gam)*SZ6*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ4 = lambdaB_ + (1-gam)*SZ5*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ3 = lambdaB_ + (1-gam)*SZ4*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ2 = lambdaB_ + (1-gam)*SZ3*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    SZ1 = lambdaB_ + (1-gam)*SZ2*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
    @#for lag in [1:7]
    prodR@{lag} = R_^@{lag};
    MD@{lag}_ = SZ@{lag}*R_^@{lag}*( (1-gam)*(1-Theta) )/( (1-(1-gam)*(1-Theta)) - lambdaB_ );
    mD@{lag}_ = ( MD@{lag}_ + (1-gam)*(1-Theta)*R_^@{lag} )/( 1 - (1-gam)*(1-Theta) );
    @#endfor
    kappa = (1-gam)*betta*MD1_/(1+(1-gam)*betta*MD1_);
    MV_ = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam);
    mv = log( MV_ );
    mV_ = MV_ / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
    RK_ = 1/((1-gam)*betta*MV_);
    Z_ = RK_-1+deltaSS;
    K_over_Y = alp/Z_;
    I_over_Y = deltaSS * K_over_Y;
    C_over_Y = 1 - gySS - I_over_Y;
    H = H_bar;
    H_ = H;
    Y_ = ( Ass * H_ ) * K_over_Y ^ ( alp / ( 1 - alp ) );
    K_ = K_over_Y * Y_;
    k = log( K_ );
    inv = log(deltaSS) + k;
    c = log( C_over_Y * Y_ );
    sum_mD = mD1_ + mD2_ + mD3_ + mD4_ + mD5_ + mD6_ + mD7_;
    sum_MD = MD1_ + MD2_ + MD3_ + MD4_ + MD5_ + MD6_ + MD7_;
    AUXBD = sum_mD/(1+mV_*R_/(1-kappa));
    AUXBS = (mV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa));
    AUXED = (1/nubar)*( log( 1 - kappa/kappa_new ) )*( ((MV_*R_/(1-kappa))/(1+mV_*R_/(1-kappa)))*sum_mD - sum_MD );
    AUXES = (1/nubar)*( log( 1 - kappa/kappa_new ) )*((MV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa)));
    D = K_*(RK_ - 1 - (R_-1)*AUXBS - (1-kappa)*AUXES)/(1-(1-kappa)*AUXED + (R_-1)*AUXBD);
    E_ = AUXED*D - AUXES*K_;
    B_ = AUXBD*D + AUXBS*K_;
    b = log( B_ );
    @#for lag in [1:7]
    lagD@{lag} = D;
    @#endfor
    spread = RK_ - R_;

    Y = Y_;
    //I = exp(inv);
    //C = exp(c);
    //R = R_;
    //K = K_;
    //Q = 1;
    //E = E_;
    //B = B_;
    E_rate = E_/(K_-B_);
    D_rate = D/(K_-B_);
end;

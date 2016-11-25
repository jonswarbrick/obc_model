steady_state_model;
    psi = 1;
    q = 0;
    R_ = 1 / betta;
    r = log( R_ );
    %Z_ = 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * ( betta * Theta * ( 1 - deltaSS ) );
    Z_ = 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) );
    a = log( Ass );
    g = log( GSS );
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
    logit_delta = logit_deltaSS;
    %B_ = ( 1 - ( deltaSS + betta * Theta * ( 1 - gam ) * ( 1 - deltaSS ) ) ) / ( 1 - gam ) * K_;
    B_ = ( 1 - ( deltaSS + Theta * ( 1 - gam ) * ( 1 - deltaSS ) ) ) / ( 1 - gam ) * K_;
    r_PLUS_b = log( R_ * B_ );
    jB = 0;
    jK = -log( 1 - gam );
    mK = jK + r;
    spread = Z_ + 1 - deltaSS - R_;
    RK_ = spread + R_;
    
    Y = Y_;
    //I = exp(inv);
    //C = exp(c);
    //R = R_;
    //K = K_;
    //B = K_;
    //Q = 1;
    D = RK_*K_ - K_ - (R_-1)*B_;
    //E = 0;
end;
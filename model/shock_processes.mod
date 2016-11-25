
a-log(Ass) = rhoA*(a(-1)-log(Ass))+sigma_a*epsA;
g-log(GSS) = rhoG*(g(-1)-log(GSS))-sigma_g*epsG;
%g = log(GSS);

@#if shock_choice == 1
    psi = ( exp(sigma_psi*eps_psi^3) );
    logit_delta = logit_deltabar;
@#endif
@#if shock_choice == 2
    %logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*epsdelta;
    logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*epsdelta^3;
    %logit_delta-logit_deltabar = rhodelta*(logit_delta(-1)-logit_deltabar)+sigma_delta*(epsdelta^2-1);
    psi = 1;
@#endif
@#if shock_choice == 4
    logit_delta = logit_deltabar;
    psi = 1;
@#endif


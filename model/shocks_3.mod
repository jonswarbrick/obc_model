
% which shocks to compute IRFs for - comma separated: epsA , epsG , epsdelta , eps_psi
@#if shock_choice == 1
    eps_psi %, epsM
@#endif
@#if shock_choice == 2
    epsdelta %, epsM
@#endif
@#if shock_choice == 4
    epsA
@#endif

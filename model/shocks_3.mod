
% which shocks to compute IRFs for - comma separated: epsA , epsG , epsdelta , eps_psi
@#if shock_choice == 1
    eps_psi
@#endif
@#if shock_choice == 2
    epsdelta
@#endif
@#if shock_choice == 4
    epsA
@#endif

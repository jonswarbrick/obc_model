% declare shocks

@#if shock_choice == 1
    varexo epsA eps_psi;
@#endif
@#if shock_choice == 2
    varexo epsA epsdelta;
@#endif
@#if shock_choice == 4
    varexo epsA;
@#endif

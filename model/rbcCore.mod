
spread = lead_RK - R;

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


lead_Lambda*R/lead_Pi=1; 
log(Y) = (1-alp)*((a) + log(H)) + alp*( k(-1) + log(psi)) - disp; 

UH/UC = - W;

close all;
Y = (oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:))./(mean(oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:)));
y = log(oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:));
inv = (oo_.endo_simul(strmatch('inv',M_.endo_names,'exact'),:));
c = (oo_.endo_simul(strmatch('c',M_.endo_names,'exact'),:));
I = exp(inv)./(mean(exp(inv)));
C = exp(c)./(mean(exp(c)));
[~,hp.Y] = hpfilter(Y,1600);
[~,hp.I] = hpfilter(I,1600);
spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));
K = exp((oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),:)));
Q = exp((oo_.endo_simul(strmatch('q',M_.endo_names,'exact'),:)));
[~,hp.C] = hpfilter(C,1600);
[~,hp.y] = hpfilter(y,1600);
[~,hp.inv] = hpfilter(inv,1600);
[~,hp.c] = hpfilter(c,1600);
if strcmp(dynareOBC_.BaseFileName,'obc')
B = exp((oo_.endo_simul(strmatch('b',M_.endo_names,'exact'),:)));
D_rate = ((oo_.endo_simul(strmatch('D_rate',M_.endo_names,'exact'),:)));
E_rate = ((oo_.endo_simul(strmatch('E_rate',M_.endo_names,'exact'),:)));
elseif strcmp(dynareOBC_.BaseFileName,'rbc')
B = Q .* K;
D_rate = zeros(1,length(Y));
E_rate = zeros(1,length(Y));
elseif strcmp(dynareOBC_.BaseFileName,'gk')
B = Q .* K .* ( 1 -  Theta ./ ( exp((oo_.endo_simul(strmatch('m',M_.endo_names,'exact'),:))))); 
D_rate = ((oo_.endo_simul(strmatch('D_rate',M_.endo_names,'exact'),:)));
E_rate = zeros(1,length(Y));
end
D = ((oo_.endo_simul(strmatch('D',M_.endo_names,'exact'),:)));
N = Q.*K-B;
[y_ac,~,~] = autocorr(hp.Y,1);
disp(horzcat('S.D. Y = ',num2str(std(hp.Y)),'| Target = 0.010563'));
disp(horzcat('AC(1) Y = ',num2str(y_ac(2)),'| Target = 0.86255'));
disp(horzcat('mean spread = ',num2str(mean(spread)),'| Target = 0.0057369'));
disp(horzcat('S.D. spread = ',num2str(std(spread)),'| Target = 0.0017812'));
disp(horzcat('Investment skewness: ',num2str(skewness(I))));
disp(horzcat('Spread skewness: ',num2str(skewness(spread))));
disp(horzcat('Capital-asset ratio: ',num2str(mean(N./(Q.*K)))));

% Table for paper
variables = {'Y','I','C','D','E','spread'};
data = [ hp.Y.*100 , hp.I.*100 , hp.C.*100 , D_rate'.*100 , E_rate'.*100 , spread'.*100  ];
data2 = [ hp.y.*100 , hp.inv.*100 , hp.c.*100 , D_rate'.*100 , E_rate'.*100 , spread'.*100  ];
data3 = [ Y'.*100 , I'.*100 , C'.*100 , D_rate'.*100 , E_rate'.*100 , spread'.*100  ];
corr_coeff = corrcoef(data);
corr_coeff2 = corrcoef(data2);
corr_coeff3 = corrcoef(data3);
moments = [ mean(data)' std(data)' skewness(data)' kurtosis(data)' ];
moments2 = [ mean(data2)' std(data2)' skewness(data2)' kurtosis(data2)' ];
moments3 = [ mean(data3)' std(data3)' skewness(data3)' kurtosis(data3)' ];

paper_data = [ ...
    corr_coeff(1,1) , corr_coeff(1,2) , corr_coeff(1,3) ,...
    corr_coeff(1,4) , corr_coeff(1,5) , corr_coeff(1,6) ; ...
    ( moments(:,2) )' ; ( moments(:,3) )' ]; 


paper_data2 = [ ...
    corr_coeff2(1,1) , corr_coeff2(1,2) , corr_coeff2(1,3) ,...
    corr_coeff2(1,4) , corr_coeff2(1,5) , corr_coeff2(1,6) ; ...
    ( moments2(:,2) )' ; ( moments2(:,3) )' ]; 


paper_data3 = [ ...
    corr_coeff3(1,1) , corr_coeff3(1,2) , corr_coeff3(1,3) ,...
    corr_coeff3(1,4) , corr_coeff3(1,5) , corr_coeff3(1,6) ; ...
    ( moments3(:,2) )' ; ( moments3(:,3) )' ]; 


tab_paper = table(paper_data(:,1),paper_data(:,2),paper_data(:,3),...
    paper_data(:,4),paper_data(:,5),paper_data(:,6),...
    'VariableNames',variables,'RowNames',{'corr','sd','skew'});

tab_paper2 = table(paper_data2(:,1),paper_data2(:,2),paper_data2(:,3),...
    paper_data2(:,4),paper_data2(:,5),paper_data2(:,6),...
    'VariableNames',variables,'RowNames',{'corr','sd','skew'});


tab_paper3 = table(paper_data3(:,1),paper_data3(:,2),paper_data3(:,3),...
    paper_data3(:,4),paper_data3(:,5),paper_data3(:,6),...
    'VariableNames',variables,'RowNames',{'corr','sd','skew'});

disp('**-------------------------------------------------------------------------------**');
disp('Exact percentage');
disp(tab_paper);
disp('**-------------------------------------------------------------------------------**');


disp('**-------------------------------------------------------------------------------**');
disp('Logs');
disp(tab_paper2);
disp('**-------------------------------------------------------------------------------**');


disp('**-------------------------------------------------------------------------------**');
disp('unfiltered');
disp(tab_paper3);
disp('**-------------------------------------------------------------------------------**');

Y = ((oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:)));
R = exp((oo_.endo_simul(strmatch('r',M_.endo_names,'exact'),:)));
if strcmp(dynareOBC_.BaseFileName,'obc')
            mv = (oo_.endo_simul(strmatch('mv',M_.endo_names,'exact'),:));
            mv(mv<1e-8) = 0;
            mv(mv>1e-8) = 1;
            binding_periods = mean(mv);
disp(horzcat('Constraint binding in ',num2str(100*binding_periods),'% of periods'));
kappa = ((oo_.endo_simul(strmatch('kappa',M_.endo_names,'exact'),:)));
psi = ((oo_.endo_simul(strmatch('psi',M_.endo_names,'exact'),:)));
MV = exp((oo_.endo_simul(strmatch('mv',M_.endo_names,'exact'),:)));
ZLB1 = ((oo_.endo_simul(strmatch('dynareOBCZeroLowerBounded1',M_.endo_names,'exact'),:)));
ZLB2 = ((oo_.endo_simul(strmatch('dynareOBCZeroLowerBounded2',M_.endo_names,'exact'),:)));
prodR1 = ((oo_.endo_simul(strmatch('prodR1',M_.endo_names,'exact'),:)));
prodR2 = ((oo_.endo_simul(strmatch('prodR2',M_.endo_names,'exact'),:)));
prodR3 = ((oo_.endo_simul(strmatch('prodR3',M_.endo_names,'exact'),:)));
prodR4 = ((oo_.endo_simul(strmatch('prodR4',M_.endo_names,'exact'),:)));
prodR5 = ((oo_.endo_simul(strmatch('prodR5',M_.endo_names,'exact'),:)));
prodR6 = ((oo_.endo_simul(strmatch('prodR6',M_.endo_names,'exact'),:)));
prodR7 = ((oo_.endo_simul(strmatch('prodR7',M_.endo_names,'exact'),:)));
SZ1 = ((oo_.endo_simul(strmatch('SZ1',M_.endo_names,'exact'),:)));
SZ2 = ((oo_.endo_simul(strmatch('SZ2',M_.endo_names,'exact'),:)));
SZ3 = ((oo_.endo_simul(strmatch('SZ3',M_.endo_names,'exact'),:)));
SZ4 = ((oo_.endo_simul(strmatch('SZ4',M_.endo_names,'exact'),:)));
SZ5 = ((oo_.endo_simul(strmatch('SZ5',M_.endo_names,'exact'),:)));
SZ6 = ((oo_.endo_simul(strmatch('SZ6',M_.endo_names,'exact'),:)));
SZ7 = ((oo_.endo_simul(strmatch('SZ7',M_.endo_names,'exact'),:)));
lagD1 = ((oo_.endo_simul(strmatch('lagD1',M_.endo_names,'exact'),:)));
lagD2 = ((oo_.endo_simul(strmatch('lagD2',M_.endo_names,'exact'),:)));
lagD3 = ((oo_.endo_simul(strmatch('lagD3',M_.endo_names,'exact'),:)));
lagD4 = ((oo_.endo_simul(strmatch('lagD4',M_.endo_names,'exact'),:)));
lagD5 = ((oo_.endo_simul(strmatch('lagD5',M_.endo_names,'exact'),:)));
lagD6 = ((oo_.endo_simul(strmatch('lagD6',M_.endo_names,'exact'),:)));
lagD7 = ((oo_.endo_simul(strmatch('lagD7',M_.endo_names,'exact'),:)));
Vhat = ones(1,length(Y));
RK = ones(1,length(Y));
K_0 = exp((oo_.steady_state(strmatch('k',M_.endo_names,'exact'),:)));
Q_0 = exp((oo_.steady_state(strmatch('q',M_.endo_names,'exact'),:)));
B_0 = exp((oo_.steady_state(strmatch('b',M_.endo_names,'exact'),:)));
R_0 = exp((oo_.steady_state(strmatch('r',M_.endo_names,'exact'),:)));
RK(1) = psi(1) * ( alp * Y(1) / (psi(1) * K_0) + (1-deltabar)*Q(1)) / Q_0;
Vhat(1) = ( RK(1) * Q_0 * K_0 - R_0 * B_0 ) ./ (1-kappa(1));
for t=2:length(Y)
RK(t) = psi(t) * ( alp * Y(t) / (psi(t) * K(t-1)) + (1-deltabar)*Q(t)) / Q(t-1);
Vhat(t) = ( RK(t) * Q(t-1) * K(t-1) - R(t-1) * B(t-1) ) ./ (1-kappa(t));
end
E = ( D + Q .* K - B )./(1-kappa) - Vhat;
lambdaB = (1-(1-gam).*(1-Theta)).*(1-kappa).*(MV-1)./(MV-(1-kappa).*(1-(1-gam).*(1-Theta)));
MD1 = SZ1.* prodR1 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD1 = ( MD1 + (1-gam)*(1-Theta)*prodR1 )/( 1 - (1-gam)*(1-Theta) );
MD2 = SZ2.* prodR2 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD2 = ( MD2 + (1-gam)*(1-Theta)*prodR2 )/( 1 - (1-gam)*(1-Theta) );
MD3 = SZ3.* prodR3 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD3 = ( MD3 + (1-gam)*(1-Theta)*prodR3 )/( 1 - (1-gam)*(1-Theta) );
MD4 = SZ4.* prodR4 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD4 = ( MD4 + (1-gam)*(1-Theta)*prodR4 )/( 1 - (1-gam)*(1-Theta) );
MD5 = SZ5.* prodR5 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD5 = ( MD5 + (1-gam)*(1-Theta)*prodR5 )/( 1 - (1-gam)*(1-Theta) );
MD6 = SZ6.* prodR6 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD6 = ( MD6 + (1-gam)*(1-Theta)*prodR6 )/( 1 - (1-gam)*(1-Theta) );
MD7 = SZ7.* prodR7 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB );
mD7 = ( MD7 + (1-gam)*(1-Theta)*prodR7 )/( 1 - (1-gam)*(1-Theta) );
mV = MV ./ ( 1 - (1-gam).*(1-Theta) ) - (1-kappa);
summ_times_D = mD1.*lagD1 + mD2.*lagD2 + mD3.*lagD3 + mD4.*lagD4 + mD5.*lagD5 + mD6.*lagD6 + mD7.*lagD7;
ZLB_LHS2 = mV .* Vhat + summ_times_D - B;
ZLB_RHS2 = lambdaB;
plot_rows = 3;
plot_cols = 3;
time_st = 300;
time_end = 400;
h = figure;
subplot(plot_rows,plot_cols,1); plot(N(time_st:time_end)); title('N');
subplot(plot_rows,plot_cols,2); plot(Y(time_st:time_end)); title('Y');
subplot(plot_rows,plot_cols,3); plot(spread(time_st:time_end)); title('spread');
subplot(plot_rows,plot_cols,4); plot(D(time_st:time_end)./N(time_st:time_end)); title('D rate');
subplot(plot_rows,plot_cols,5); plot(E(time_st:time_end)./N(time_st:time_end)); title('E rate');
subplot(plot_rows,plot_cols,6); plot(MV(time_st:time_end)); title('MV');
subplot(plot_rows,plot_cols,7); plot(lambdaB(time_st:time_end)); title('lambdaB');
subplot(plot_rows,plot_cols,8); plot(N(time_st:time_end)./(Q(time_st:time_end).*K(time_st:time_end))); title('Capital-Asset');
subplot(plot_rows,plot_cols,9); plot(kappa(time_st:time_end)); title('kappa');
g = figure;
subplot(plot_rows,plot_cols,1); plot(ZLB1(time_st:time_end)); title('dynareOBCZLB1');
subplot(plot_rows,plot_cols,2); plot(ZLB2(time_st:time_end)); title('dynareOBCZLB2');
subplot(plot_rows,plot_cols,3); plot(D(time_st:time_end)); title('ZLB LHS1');
subplot(plot_rows,plot_cols,4); plot(ZLB_LHS2(time_st:time_end)); title('ZLB LHS2');
subplot(plot_rows,plot_cols,5); plot(ZLB_RHS2(time_st:time_end)); title('ZLB RHS2 (lambdaB)');
end



%%
% pi = 0;
% mc = 0;
% disp = 0;
% psi = 1;
% logit_delta = logit_deltaSS;
% nu = nubar;
% q = 0;
% R_ = 1 / betta;
% r = log( R_ );
% a = log( Ass );
% lambdaD_ = 0;
% lambdaB_ = gam*(1-(1-gam)*(1-Theta));
% SZ7 = lambdaB_;
% SZ6 = lambdaB_ + (1-gam)*SZ7*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
% SZ5 = lambdaB_ + (1-gam)*SZ6*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
% SZ4 = lambdaB_ + (1-gam)*SZ5*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
% SZ3 = lambdaB_ + (1-gam)*SZ4*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
% SZ2 = lambdaB_ + (1-gam)*SZ3*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
% SZ1 = lambdaB_ + (1-gam)*SZ2*(1-(1-gam)*(1-Theta))/( (1-(1-gam)*(1-Theta))-lambdaB_ );
%
% MD1_ = SZ1.* R_^1 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD2_ = SZ2.* R_^2 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD3_ = SZ3.* R_^3 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD4_ = SZ4.* R_^4 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD5_ = SZ5.* R_^5 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD6_ = SZ6.* R_^6 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
% MD7_ = SZ7.* R_^7 .* ( (1-gam) .* (1-Theta) ) ./( (1-(1-gam) .* (1-Theta)) - lambdaB_ );
%
% mD1_ = ( MD1_ + (1-gam)*(1-Theta)*R_^1 )/( 1 - (1-gam)*(1-Theta) );
% mD2_ = ( MD2_ + (1-gam)*(1-Theta)*R_^2 )/( 1 - (1-gam)*(1-Theta) );
% mD3_ = ( MD3_ + (1-gam)*(1-Theta)*R_^3 )/( 1 - (1-gam)*(1-Theta) );
% mD4_ = ( MD4_ + (1-gam)*(1-Theta)*R_^4 )/( 1 - (1-gam)*(1-Theta) );
% mD5_ = ( MD5_ + (1-gam)*(1-Theta)*R_^5 )/( 1 - (1-gam)*(1-Theta) );
% mD6_ = ( MD6_ + (1-gam)*(1-Theta)*R_^6 )/( 1 - (1-gam)*(1-Theta) );
% mD7_ = ( MD7_ + (1-gam)*(1-Theta)*R_^7 )/( 1 - (1-gam)*(1-Theta) );
%
%
% kappa = (1-gam)*betta*MD1_/(1+(1-gam)*betta*MD1_);
% MV_ = ( 1 - gam*(1-kappa)*(1-(1-gam)*(1-Theta)) )/(1-gam);
% mv = log( MV_ );
% mV_ = MV_ / ( 1 - (1-gam)*(1-Theta) ) - (1-kappa);
% RK_ = 1/((1-gam)*betta*MV_);
% Z_ = RK_-1+deltaSS;
% K_over_Y = alp/Z_;
% I_over_Y = deltaSS * K_over_Y;
% C_over_Y = 1 - I_over_Y;
% H = H_bar;
% H_ = H;
% Y_ = ( Ass * H_ ) * K_over_Y ^ ( alp / ( 1 - alp ) );
% K_ = K_over_Y * Y_;
% k = log( K_ );
% inv = log(deltaSS) + k;
% C_ = C_over_Y * Y_;
% c = log( C_ );
% Xjr = C_;
% UX_ = - (C_ - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
% lambdaX = UX_ / ( 1 - betta*( (1-gam_jr) ) );
% sum_mD = mD1_ + mD2_ + mD3_ + mD4_ + mD5_ + mD6_ + mD7_;
% sum_MD = MD1_ + MD2_ + MD3_ + MD4_ + MD5_ + MD6_ + MD7_;
% AUXBD = sum_mD/(1+mV_*R_/(1-kappa));
% AUXBS = (mV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa));
% AUXED = (1/nubar)*( log( 1 - kappa/kappa_new ) )*( ((MV_*R_/(1-kappa))/(1+mV_*R_/(1-kappa)))*sum_mD - sum_MD );
% AUXES = (1/nubar)*( log( 1 - kappa/kappa_new ) )*((MV_*RK_/(1-kappa))/(1+mV_*R_/(1-kappa)));
% D = K_*(RK_ - 1 - (R_-1)*AUXBS - (1-kappa)*AUXES)/(1-(1-kappa)*AUXED + (R_-1)*AUXBD);
% E_ = AUXED*D - AUXES*K_;
% B_ = AUXBD*D + AUXBS*K_;
% b = log( B_ );
%
% spread = RK_ - R_;
%
% Y = Y_;
% E_rate = E_/(K_-B_);
% D_rate = D/(K_-B_);
%
% Vhat_ = ( RK_*K_ - (R_)*B_ )/(1-kappa);
% mV_*Vhat_ + sum_mD*D - B_
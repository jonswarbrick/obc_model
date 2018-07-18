load('results_sim/obc.mat');

lag_K = exp(oo_.endo_simul(1,1:end-1));
lag_Q = exp(oo_.endo_simul(4,1:end-1));
lag_R = exp(oo_.endo_simul(6,1:end-1));
lag_B = exp(oo_.endo_simul(15,1:end-1));
Q = exp(oo_.endo_simul(4,2:end));
Y = exp(oo_.endo_simul(13,2:end));
logit_delta = exp(oo_.endo_simul(3,2:end));
kappa = (oo_.endo_simul(19,2:end));
psi = exp(oo_.endo_simul(10,2:end));
Z = alp.*Y./(psi.*lag_K);
delta = 1./(1+exp(-logit_delta));
RK = psi.*(Z + (1-delta).*Q)./lag_Q;

Vhat = ( RK.*lag_Q.*lag_K - lag_R.*lag_B ) ./ (1-kappa);

plot(Vhat);

%%
close all; clear; load('results_sim/obc.mat');
lag_K = exp(oo_.endo_simul(1,1:end-1));
lag_Q = exp(oo_.endo_simul(4,1:end-1));
MV = exp(oo_.endo_simul(17,2:end));
kappa = exp(oo_.endo_simul(19,2:end));
Q = exp(oo_.endo_simul(4,2:end));
K = exp(oo_.endo_simul(1,2:end));
D = (oo_.endo_simul(16,2:end));
A = MV ./ ( 1 - (1-gam).*(1-Theta) ) - (1-kappa);
weight_1 = MV./(1+A-kappa);
Y = exp(oo_.endo_simul(13,2:end));
psi = exp(oo_.endo_simul(10,2:end));
Z = alp.*Y./(psi.*lag_K);
logit_delta = exp(oo_.endo_simul(3,2:end));
delta = 1./(1+exp(-logit_delta));
RK = psi.*(Z + (1-delta).*Q)./lag_Q;
lag_R = exp(oo_.endo_simul(6,1:end-1));
B =  exp(oo_.endo_simul(15,2:end));
lag_B =  exp(oo_.endo_simul(15,1:end-1));
Vhat = ( RK.*lag_Q.*lag_K - lag_R.*lag_B )./(1-kappa);
E = ( D + Q.*K - B )./(1-kappa) - Vhat;

term_1 = weight_1.*(Q.*K + D - E.*(1-kappa));
theta1 = term_1./(Q.*K);

figure; 
subplot(4,1,1); plot(MV(5000:5500))
subplot(4,1,2); plot(A(5000:5500))
subplot(4,1,3); plot(theta1(5000:5500))
subplot(4,1,4); plot(Y(5000:5500))

MV_hat = MV./mean(MV);
theta1_hat = theta1./mean(theta1);

figure;
plot(MV_hat(5000:5500)); hold on; plot(theta1_hat(5000:5500))
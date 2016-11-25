%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'GKQ_model';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('GKQ_model.log');
M_.exo_names = 'epsA';
M_.exo_names_tex = 'epsA';
M_.exo_names_long = 'epsA';
M_.exo_names = char(M_.exo_names, 'eps_psi');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_psi');
M_.exo_names_long = char(M_.exo_names_long, 'eps_psi');
M_.endo_names = 'k';
M_.endo_names_tex = 'k';
M_.endo_names_long = 'k';
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names_long = char(M_.endo_names_long, 'a');
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'g');
M_.endo_names = char(M_.endo_names, 'logit_delta');
M_.endo_names_tex = char(M_.endo_names_tex, 'logit\_delta');
M_.endo_names_long = char(M_.endo_names_long, 'logit_delta');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'disp');
M_.endo_names_tex = char(M_.endo_names_tex, 'disp');
M_.endo_names_long = char(M_.endo_names_long, 'disp');
M_.endo_names = char(M_.endo_names, 'psi');
M_.endo_names_tex = char(M_.endo_names_tex, 'psi');
M_.endo_names_long = char(M_.endo_names_long, 'psi');
M_.endo_names = char(M_.endo_names, 'spread');
M_.endo_names_tex = char(M_.endo_names_tex, 'spread');
M_.endo_names_long = char(M_.endo_names_long, 'spread');
M_.endo_names = char(M_.endo_names, 'inv');
M_.endo_names_tex = char(M_.endo_names_tex, 'inv');
M_.endo_names_long = char(M_.endo_names_long, 'inv');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'H');
M_.endo_names_tex = char(M_.endo_names_tex, 'H');
M_.endo_names_long = char(M_.endo_names_long, 'H');
M_.endo_names = char(M_.endo_names, 'mue');
M_.endo_names_tex = char(M_.endo_names_tex, 'mue');
M_.endo_names_long = char(M_.endo_names_long, 'mue');
M_.endo_names = char(M_.endo_names, 'mus');
M_.endo_names_tex = char(M_.endo_names_tex, 'mus');
M_.endo_names_long = char(M_.endo_names_long, 'mus');
M_.endo_names = char(M_.endo_names, 'nub');
M_.endo_names_tex = char(M_.endo_names_tex, 'nub');
M_.endo_names_long = char(M_.endo_names_long, 'nub');
M_.endo_names = char(M_.endo_names, 'E');
M_.endo_names_tex = char(M_.endo_names_tex, 'E');
M_.endo_names_long = char(M_.endo_names_long, 'E');
M_.endo_names = char(M_.endo_names, 'Thetax');
M_.endo_names_tex = char(M_.endo_names_tex, 'Thetax');
M_.endo_names_long = char(M_.endo_names_long, 'Thetax');
M_.endo_names = char(M_.endo_names, 'QE');
M_.endo_names_tex = char(M_.endo_names_tex, 'QE');
M_.endo_names_long = char(M_.endo_names_long, 'QE');
M_.endo_names = char(M_.endo_names, 'D_rate');
M_.endo_names_tex = char(M_.endo_names_tex, 'D\_rate');
M_.endo_names_long = char(M_.endo_names_long, 'D_rate');
M_.endo_names = char(M_.endo_names, 'E_rate');
M_.endo_names_tex = char(M_.endo_names_tex, 'E\_rate');
M_.endo_names_long = char(M_.endo_names_long, 'E_rate');
M_.param_names = 'gySS';
M_.param_names_tex = 'gySS';
M_.param_names_long = 'gySS';
M_.param_names = char(M_.param_names, 'varrho');
M_.param_names_tex = char(M_.param_names_tex, 'varrho');
M_.param_names_long = char(M_.param_names_long, 'varrho');
M_.param_names = char(M_.param_names, 'alp');
M_.param_names_tex = char(M_.param_names_tex, 'alp');
M_.param_names_long = char(M_.param_names_long, 'alp');
M_.param_names = char(M_.param_names, 'zzeta');
M_.param_names_tex = char(M_.param_names_tex, 'zzeta');
M_.param_names_long = char(M_.param_names_long, 'zzeta');
M_.param_names = char(M_.param_names, 'betta');
M_.param_names_tex = char(M_.param_names_tex, 'betta');
M_.param_names_long = char(M_.param_names_long, 'betta');
M_.param_names = char(M_.param_names, 'deltaSS');
M_.param_names_tex = char(M_.param_names_tex, 'deltaSS');
M_.param_names_long = char(M_.param_names_long, 'deltaSS');
M_.param_names = char(M_.param_names, 'sigma_c');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_c');
M_.param_names_long = char(M_.param_names_long, 'sigma_c');
M_.param_names = char(M_.param_names, 'rhodelta');
M_.param_names_tex = char(M_.param_names_tex, 'rhodelta');
M_.param_names_long = char(M_.param_names_long, 'rhodelta');
M_.param_names = char(M_.param_names, 'rhoG');
M_.param_names_tex = char(M_.param_names_tex, 'rhoG');
M_.param_names_long = char(M_.param_names_long, 'rhoG');
M_.param_names = char(M_.param_names, 'rhoA');
M_.param_names_tex = char(M_.param_names_tex, 'rhoA');
M_.param_names_long = char(M_.param_names_long, 'rhoA');
M_.param_names = char(M_.param_names, 'Ass');
M_.param_names_tex = char(M_.param_names_tex, 'Ass');
M_.param_names_long = char(M_.param_names_long, 'Ass');
M_.param_names = char(M_.param_names, 'Phi');
M_.param_names_tex = char(M_.param_names_tex, 'Phi');
M_.param_names_long = char(M_.param_names_long, 'Phi');
M_.param_names = char(M_.param_names, 'sigmaB');
M_.param_names_tex = char(M_.param_names_tex, 'sigmaB');
M_.param_names_long = char(M_.param_names_long, 'sigmaB');
M_.param_names = char(M_.param_names, 'xiB');
M_.param_names_tex = char(M_.param_names_tex, 'xiB');
M_.param_names_long = char(M_.param_names_long, 'xiB');
M_.param_names = char(M_.param_names, 'Theta');
M_.param_names_tex = char(M_.param_names_tex, 'Theta');
M_.param_names_long = char(M_.param_names_long, 'Theta');
M_.param_names = char(M_.param_names, 'kappaSS');
M_.param_names_tex = char(M_.param_names_tex, 'kappaSS');
M_.param_names_long = char(M_.param_names_long, 'kappaSS');
M_.param_names = char(M_.param_names, 'logit_deltaSS');
M_.param_names_tex = char(M_.param_names_tex, 'logit\_deltaSS');
M_.param_names_long = char(M_.param_names_long, 'logit_deltaSS');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names_long = char(M_.param_names_long, 'epsilon');
M_.param_names = char(M_.param_names, 'kappa_GK');
M_.param_names_tex = char(M_.param_names_tex, 'kappa\_GK');
M_.param_names_long = char(M_.param_names_long, 'kappa_GK');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, 'gam');
M_.param_names_long = char(M_.param_names_long, 'gam');
M_.param_names = char(M_.param_names, 'sigma_g');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_g');
M_.param_names_long = char(M_.param_names_long, 'sigma_g');
M_.param_names = char(M_.param_names, 'sigma_a');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_a');
M_.param_names_long = char(M_.param_names_long, 'sigma_a');
M_.param_names = char(M_.param_names, 'sigma_psi');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_psi');
M_.param_names_long = char(M_.param_names_long, 'sigma_psi');
M_.param_names = char(M_.param_names, 'sigma_delta');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_delta');
M_.param_names_long = char(M_.param_names_long, 'sigma_delta');
M_.param_names = char(M_.param_names, 'deltabar');
M_.param_names_tex = char(M_.param_names_tex, 'deltabar');
M_.param_names_long = char(M_.param_names_long, 'deltabar');
M_.param_names = char(M_.param_names, 'logit_deltabar');
M_.param_names_tex = char(M_.param_names_tex, 'logit\_deltabar');
M_.param_names_long = char(M_.param_names_long, 'logit_deltabar');
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names_long = char(M_.param_names_long, 'chi');
M_.param_names = char(M_.param_names, 'sigma_h');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_h');
M_.param_names_long = char(M_.param_names_long, 'sigma_h');
M_.param_names = char(M_.param_names, 'psi_h');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_h');
M_.param_names_long = char(M_.param_names_long, 'psi_h');
M_.param_names = char(M_.param_names, 'epsilonC');
M_.param_names_tex = char(M_.param_names_tex, 'epsilonC');
M_.param_names_long = char(M_.param_names_long, 'epsilonC');
M_.param_names = char(M_.param_names, 'epsilonH');
M_.param_names_tex = char(M_.param_names_tex, 'epsilonH');
M_.param_names_long = char(M_.param_names_long, 'epsilonH');
M_.param_names = char(M_.param_names, 'gam_jr');
M_.param_names_tex = char(M_.param_names_tex, 'gam\_jr');
M_.param_names_long = char(M_.param_names_long, 'gam_jr');
M_.param_names = char(M_.param_names, 'theta_jr');
M_.param_names_tex = char(M_.param_names_tex, 'theta\_jr');
M_.param_names_long = char(M_.param_names_long, 'theta_jr');
M_.param_names = char(M_.param_names, 'C_bar');
M_.param_names_tex = char(M_.param_names_tex, 'C\_bar');
M_.param_names_long = char(M_.param_names_long, 'C_bar');
M_.param_names = char(M_.param_names, 'K_by_Y');
M_.param_names_tex = char(M_.param_names_tex, 'K\_by\_Y');
M_.param_names_long = char(M_.param_names_long, 'K_by_Y');
M_.param_names = char(M_.param_names, 'H_bar');
M_.param_names_tex = char(M_.param_names_tex, 'H\_bar');
M_.param_names_long = char(M_.param_names_long, 'H_bar');
M_.param_names = char(M_.param_names, 'GSS');
M_.param_names_tex = char(M_.param_names_tex, 'GSS');
M_.param_names_long = char(M_.param_names_long, 'GSS');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 23;
M_.param_nbr = 37;
M_.orig_endo_nbr = 23;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('GKQ_model_static');
erase_compiled_function('GKQ_model_dynamic');
M_.lead_lag_incidence = [
 1 13 36;
 2 14 0;
 0 15 37;
 0 16 38;
 3 17 39;
 4 18 40;
 5 19 0;
 0 20 41;
 0 21 42;
 0 22 43;
 0 23 44;
 0 24 0;
 0 25 45;
 0 26 0;
 6 27 46;
 7 28 0;
 8 29 0;
 9 30 0;
 10 31 0;
 11 32 47;
 12 33 48;
 0 34 0;
 0 35 0;]';
M_.nstatic = 4;
M_.nfwrd   = 7;
M_.npred   = 6;
M_.nboth   = 6;
M_.nsfwrd   = 13;
M_.nspred   = 12;
M_.ndynamic   = 19;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(23, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(37, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 211;
M_.NNZDerivatives(2) = 2436;
M_.NNZDerivatives(3) = 35979;
load('../opts.mat');
M_.params( 1 ) = 0.2;
gySS = M_.params( 1 );
M_.params( 2 ) = 0.684;
varrho = M_.params( 2 );
M_.params( 3 ) = 0.3;
alp = M_.params( 3 );
M_.params( 4 ) = 7.0;
zzeta = M_.params( 4 );
M_.params( 5 ) = 0.995;
betta = M_.params( 5 );
M_.params( 25 ) = 0.025;
deltabar = M_.params( 25 );
M_.params( 26 ) = log(M_.params(25)/(1-M_.params(25)));
logit_deltabar = M_.params( 26 );
M_.params( 7 ) = 2.0;
sigma_c = M_.params( 7 );
M_.params( 28 ) = parameter_sigma_h;
sigma_h = M_.params( 28 );
M_.params( 29 ) = parameter_psi_h;
psi_h = M_.params( 29 );
M_.params( 27 ) = .7;
chi = M_.params( 27 );
M_.params( 32 ) = 0.001;
gam_jr = M_.params( 32 );
M_.params( 33 ) = 1.5;
theta_jr = M_.params( 33 );
M_.params( 30 ) = parameter_habits_C;
epsilonC = M_.params( 30 );
M_.params( 31 ) = parameter_habits_H;
epsilonH = M_.params( 31 );
M_.params( 11 ) = 1;
Ass = M_.params( 11 );
M_.params( 10 ) = parameter_rhoA;
rhoA = M_.params( 10 );
M_.params( 8 ) = parameter_rhodelta;
rhodelta = M_.params( 8 );
M_.params( 9 ) = parameter_rhoG;
rhoG = M_.params( 9 );
M_.params( 13 ) = 0.975;
sigmaB = M_.params( 13 );
M_.params( 14 ) = 0.00017;
xiB = M_.params( 14 );
M_.params( 18 ) = (-2);
epsilon = M_.params( 18 );
M_.params( 19 ) = 13;
kappa_GK = M_.params( 19 );
M_.params( 20 ) = 1e-8;
gam = M_.params( 20 );
M_.params( 16 ) = parameter_kappa;
kappaSS = M_.params( 16 );
M_.params( 15 ) = parameter_Theta;
Theta = M_.params( 15 );
M_.params( 12 ) = parameter_Phi;
Phi = M_.params( 12 );
M_.params( 21 ) = parameter_sigma_g;
sigma_g = M_.params( 21 );
M_.params( 22 ) = parameter_sigma_a;
sigma_a = M_.params( 22 );
M_.params( 23 ) = parameter_sigma_psi;
sigma_psi = M_.params( 23 );
M_.params( 24 ) = parameter_sigma_delta;
sigma_delta = M_.params( 24 );
M_.params( 17 ) = M_.params(26);
logit_deltaSS = M_.params( 17 );
M_.params( 6 ) = 1/(1+exp((-M_.params(17))));
deltaSS = M_.params( 6 );
M_.params( 36 ) = call_Hbar_gkq;
H_bar = M_.params( 36 );
M_.params( 35 ) = call_KbyY_gkq;
K_by_Y = M_.params( 35 );
M_.params( 34 ) = (1-M_.params(1)-M_.params(6)*M_.params(35))*M_.params(11)*M_.params(36)*M_.params(35)^(M_.params(3)/(1-M_.params(3)));
C_bar = M_.params( 34 );
M_.params( 37 ) = M_.params(35)^(M_.params(3)/(1-M_.params(3)))*M_.params(1)*M_.params(11)*M_.params(36);
GSS = M_.params( 37 );
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
M_.Sigma_e(2, 2) = 1;
options_.irf = 60;
options_.nograph = 1;
options_.order = 3;
options_.periods = 0;
options_.replic = 500;
options_.irf_shocks=[];
options_.irf_shocks = 'eps_psi';
var_list_=[];
var_list_ = 'Y';
var_list_ = char(var_list_, 'H');
var_list_ = char(var_list_, 'inv');
var_list_ = char(var_list_, 'c');
var_list_ = char(var_list_, 'k');
var_list_ = char(var_list_, 'r');
var_list_ = char(var_list_, 'D_rate');
var_list_ = char(var_list_, 'E_rate');
var_list_ = char(var_list_, 'spread');
var_list_ = char(var_list_, 'q');
var_list_ = char(var_list_, 'S');
var_list_ = char(var_list_, 'A');
var_list_ = char(var_list_, 'delta');
var_list_ = char(var_list_, 'G');
var_list_ = char(var_list_, 'pi');
info = stoch_simul(var_list_);
save('GKQ_model_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('GKQ_model_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('GKQ_model_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('GKQ_model_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('GKQ_model_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off

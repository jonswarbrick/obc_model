%% Simulation set-up
clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 if calling run from multiple_run.m, 0 otherwise
mult_run = 1;
% 1 = simulation, 2 = irf
sim_type = 2;
% other
number_of_runs = 1;
% dynareOBC options 1: SlowIRF/no cubature, 2: FastIRF/Fast cubature, 3:
% FastIRF/no cubature, 4: FastIRF/default cubature, 5: FastIRF QuasiMC
%opts.dynareOBC_irf_options_1 = ' SlowIRFs FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=60 NoCubature omega=10000 CompileSimulationCode';
%opts.dynareOBC_irf_options_1 = ' FastCubature FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=64 omega=10000 CompileSimulationCode';
%opts.dynareOBC_irf_options_1 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 NoCubature omega=10000 CompileSimulationCode';
opts.dynareOBC_irf_options_1 = ' OrderOverride=3 FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=64 omega=10000 CompileSimulationCode';
%opts.dynareOBC_irf_options_1 = '  QuasiMonteCarloLevel=8 CubatureTolerance=0 FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=64 omega=10000 CompileSimulationCode';

opts.dynareOBC_irf_options_2 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_3 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_4 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_sim_options_2 = ' NoCubature FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=1 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_1 = ' QuasiMonteCarloLevel=8 CubatureTolerance=0 TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=1 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_3 = ' FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=1 omega=10000 CompileSimulationCode Sparse';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if mult_run == 0
models_to_run = 2;
% 1 = non-separable, 2 = additive type 1 , 3 = additive type 2 , 4 =
% non-separable habits on bundles , 5 J-R
utility_type = 5;
% 1 = KQ, 2 = delta
shock_choice = 1;
% 1 = CEE, 2 = Ireland (2003)
adj_type = 2;
% MAT-file names
opts.mat_file_string_1 = '_irfs_test';
opts.mat_file_string_2 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';
opts.mat_file_string_3 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';
opts.mat_file_string_4 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';

% Parameters
parameter_Phi = 2;
parameter_habits_C = .7;
parameter_habits_H = 0;
elseif mult_run == 1
load('mult.mat')
models_to_run = 2;
utility_type = mult.utility_type;
shock_choice = mult.shock_choice;
adj_type = mult.adj_type;
opts.mat_file_string_1 = mult.temp_file_name;
opts.mat_file_string_2 = mult.temp_file_name;
opts.mat_file_string_3 = mult.temp_file_name;
opts.mat_file_string_4 = mult.temp_file_name;
parameter_Phi = mult.parameter_Phi;
parameter_habits_C = mult.parameter_habits_C;
parameter_habits_H = mult.parameter_habits_H;
parameter_xi = mult.parameter_xi;
end

% Parameters
parameter_sigma_g = -.2;
parameter_sigma_a = -.005; 
parameter_sigma_delta = .3;
parameter_sigma_psi = -.1;
parameter_rhoA = 0.7;
parameter_rhoG = 0.95;
parameter_rhodelta = 0.85;

parameter_Theta = 0.80381;
parameter_kappa = 0.05;
parameter_kappa_new = .1;
parameter_nubar = 400;
   
parameter_gamR = 0.9;
parameter_gamPi = 2;
parameter_gamY = 0.4;

parameter_sigma_h = 2.37;
parameter_psi_h = 1;

%% Preparation

[~,numberModels] = size(models_to_run);
if sim_type == 1
opts.dynareOBC_options = {opts.dynareOBC_sim_options_1 ; opts.dynareOBC_sim_options_2 ; opts.dynareOBC_sim_options_3};
elseif sim_type == 2
opts.dynareOBC_options = {opts.dynareOBC_irf_options_1 ; opts.dynareOBC_irf_options_2 ; opts.dynareOBC_irf_options_3 ; opts.dynareOBC_irf_options_4};
end
opts.mat_file_string = {opts.mat_file_string_1 ; opts.mat_file_string_2 ; opts.mat_file_string_3 ; opts.mat_file_string_4};
opts.models = {'rbc' ; 'gkq' ; 'obc' ; 'nk' ; 'nkobc' ; 'newobc' };
opts.sim_type_name = {'sim' ; 'irf'};
if mult_run == 0
opts.loop_num = 1;
opts.total_loops = 1;
elseif mult_run == 1
opts.loop_num = mult.loop_num;
opts.total_loops = mult.total_loops;
end

save('opts.mat')

%%
cd('model');
for jj=1:number_of_runs
for ii=1:numberModels
    save('loop.mat','ii','jj')
    fid_mod = fopen( 'which_model.mod', 'wt' );
    fprintf( fid_mod, '@#define model_type = %d\n', models_to_run(ii));
    fclose(fid_mod);
    fid_sim = fopen( 'sim_type.mod', 'wt' );
    fprintf( fid_sim, '@#define sim_type = %d\n', sim_type);
    fclose(fid_sim);
    fid_util = fopen( 'utility_type.mod', 'wt' );
    fprintf( fid_util, '@#define utility_type = %d\n', utility_type);
    fclose(fid_util);
    fid_shoc = fopen( 'shock_choice.mod', 'wt' );
    fprintf( fid_shoc, '@#define shock_choice = %d\n', shock_choice);
    fclose(fid_shoc);
    fid_adj = fopen( 'adj_type.mod', 'wt' );
    fprintf( fid_adj, '@#define adj_type = %d\n', adj_type);
    fclose(fid_adj);
    
    if models_to_run(ii) == 2
        try
        movefile('GKQ_steadystate.m','all_models_steadystate.m')
        end
    end
    
    try
    eval(strcat('dynareOBC gk ',char(opts.dynareOBC_options(jj,:)),';'));
   %dynare gk
    load('loop.mat')
    fid_log = fopen( '../mutliple_log.txt', 'At' );
    log_txt = strcat(char(strcat('Loop ',num2str(opts.loop_num),'/',num2str(opts.total_loops),': ',opts.models(models_to_run(ii),:),opts.mat_file_string(jj,:))) , ' run ok!\n');
    fprintf( fid_log, log_txt );
    fclose(fid_log);
    catch
    load('loop.mat')
    fid_log = fopen( '../mutliple_log.txt', 'At' );
    log_txt = strcat(char(strcat('Loop ',num2str(opts.loop_num),'/',num2str(opts.total_loops),': ',opts.models(models_to_run(ii),:),opts.mat_file_string(jj,:))) , ' failed!\n');
    fprintf( fid_log, log_txt );
    fclose(fid_log);
    end
    save(char(strcat('../results_',opts.sim_type_name(sim_type,:),'/',opts.models(models_to_run(ii),:),opts.mat_file_string(jj,:),'.mat')));
    save(char(strcat('../results_',opts.sim_type_name(sim_type,:),'/latest_',opts.models(models_to_run(ii),:),'_results_',num2str(jj),'.mat')));

   % if models_to_run(ii) == 2
    %    movefile('all_models_steadystate.m','GKQ_steadystate.m')
   % end
end
end

%% Clean up

warning('off','all')
delete *_IRF_*
delete *Temp*
delete dynareOBCTemp*
delete *_dynamic.m
delete *_static.m
delete *.jnl
delete *.log
delete *.eps
delete *.pdf
delete *_set_auxiliary_variables.m
delete all_models_results.mat
delete all_models.m
warning('on','all')

delete('loop.mat')
cd('../');
delete('opts.mat')


%%

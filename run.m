%% Simulation set-up
clear;% close all;

if exist('results_irf')~=7
    mkdir('results_irf')
end
if exist('results_sim')~=7
    mkdir('results_sim')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 if calling run from multiple_run.m, 0 otherwise
mult_run = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mult_run == 0
% 1 = simulation, 2 = irf
sim_type = 2;
% other
number_of_runs = 1;
% dynareOBC options for each run:
opts.dynareOBC_irf_options_1 = ' slowIRFs OrderOverride=2 FirstOrderConditionalCovariance shockscale=-5000 TimeToEscapeBounds=40 omega=10000 CompileSimulationCode';
opts.dynareOBC_irf_options_1 = ' slowIRFs OrderOverride=2 FirstOrderConditionalCovariance shockscale=5000 TimeToEscapeBounds=40 omega=10000 CompileSimulationCode';
opts.dynareOBC_irf_options_1 = ' slowIRFs OrderOverride=2 FirstOrderConditionalCovariance shockscale=-1 TimeToEscapeBounds=40 omega=10000 CompileSimulationCode';
opts.dynareOBC_irf_options_1 = ' slowIRFs OrderOverride=2 FirstOrderConditionalCovariance shockscale=1 TimeToEscapeBounds=40 omega=10000 CompileSimulationCode';
opts.dynareOBC_irf_options_2 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_3 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_4 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_sim_options_1 = ' NoCubature OrderOverride=2 FirstOrderConditionalCovariance TimeToEscapeBounds=60 TimeToReturnToSteadyState=20 omega=10000 CompileSimulationCode Sparse';
%opts.dynareOBC_sim_options_1 = ' QuasiMonteCarloLevel=8 CubatureTolerance=0 OrderOverride=2 FirstOrderConditionalCovariance TimeToEscapeBounds=60 TimeToReturnToSteadyState=20 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_2 = ' QuasiMonteCarloLevel=8 CubatureTolerance=0 FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=0 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_3 = ' FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=1 omega=10000 CompileSimulationCode Sparse';

% 1 = rbc, 2 = gk, 3 = obc, 4 = nk, 5 = nkobc, 6 = newobc , 7 = gkq
models_to_run = [ 1,7 ];
% 1 = non-separable, 2 = additive type 1 , 3 = additive type 2 , 4 =
% non-separable habits on bundles , 5 J-R
utility_type = 5;
% 1 = KQ, 2 = delta, 4 = epsA  (IRF -- epsA for all for simulation)
shock_choice = 1;
% 1 = CEE, 2 = Ireland (2003)
adj_type = 1;
% MAT-file names
opts.mat_file_string_1 = '_order2';
opts.mat_file_string_1 = '_order2_median_KQ5pc';
opts.mat_file_string_1 = '_order2_median_posKQ5pc';
opts.mat_file_string_1 = '_order2_median_negA1sd';
opts.mat_file_string_1 = '_order2_median_A1sd';
opts.mat_file_string_2 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';
opts.mat_file_string_3 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';
opts.mat_file_string_4 = '_irfs_order3_X3_slow_phi4_shocksPsiA_habitC90_habitH0_sepUtilFrisch';

% 1PC KQ = -1000x, 5PC = -5000x

% Parameters
elseif mult_run == 1
load('mult.mat')
models_to_run = mult.models_to_run;
utility_type = mult.utility_type;
shock_choice = mult.shock_choice;
adj_type = mult.adj_type;
opts.mat_file_string_1 = mult.temp_file_name;
opts.mat_file_string_2 = mult.temp_file_name;
opts.mat_file_string_3 = mult.temp_file_name;
opts.mat_file_string_4 = mult.temp_file_name;
sim_type = mult.sim_type;
number_of_runs = mult.number_of_runs;
opts.dynareOBC_irf_options_1 = mult.dynareOBC_irf_options_1;
opts.dynareOBC_irf_options_2 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_3 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_irf_options_4 = ' FirstOrderConditionalCovariance shockscale=3 TimeToEscapeBounds=40 TimeToReturnToSteadyState=10 NoCubature omega=10000 MLVSimulationMode=2  CompileSimulationCode';
opts.dynareOBC_sim_options_1 = ' NoCubature OrderOverride=2 FirstOrderConditionalCovariance TimeToEscapeBounds=60 TimeToReturnToSteadyState=20 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_2 = ' QuasiMonteCarloLevel=8 CubatureTolerance=0 FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=0 omega=10000 CompileSimulationCode Sparse';
opts.dynareOBC_sim_options_3 = ' FirstOrderConditionalCovariance TimeToEscapeBounds=40 TimeToReturnToSteadyState=20 MLVSimulationMode=1 omega=10000 CompileSimulationCode Sparse';

end

% Parameters
parameter_rho_psi = 0;
parameter_Phi = 2;
parameter_kappa = 0.05;
parameter_kappa_new = .1;
parameter_nubar = 400;
parameter_sigma_c =2;
parameter_gam_jr = 0.001;
parameter_theta_jr = 1.4;

% Unused under baseline environment
parameter_habits_C = 0;
parameter_habits_H = 0;
parameter_sigma_g = 0;
parameter_sigma_delta = .9;
parameter_rhoG = 0.95;
parameter_rhodelta = 0.85;
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
opts.models = {'rbc' ; 'gk' ; 'obc' ; 'nk' ; 'nkobc' ; 'newobc' ; 'gkq' };
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
    % calibrated parameters
    if models_to_run(ii)==1
    parameter_sigma_a = 0.006145871647583;
    parameter_rhoA = 0.95;
    parameter_Theta = 0.9;
    parameter_sigma_psi = 0;
    elseif models_to_run(ii)==2
    parameter_sigma_a = 0.0056549153601962;
    parameter_rhoA = 0.95;
    parameter_Theta = 0.8449834740987057;
    parameter_sigma_psi = 0;
    elseif models_to_run(ii)==6
    parameter_sigma_a = 0.0060626464427887;
    parameter_rhoA = 0.95;
    parameter_Theta = 0.669734151867773;
    parameter_sigma_psi = 0;
    elseif models_to_run(ii)==7
    parameter_sigma_a = 0.005710743233638;
    parameter_rhoA = 0.95;
    parameter_Theta = 0.893717341213108;
    parameter_Theta = 0.84;
    parameter_sigma_psi = 0;
    end
    % compare the same shock size for IRFs
    if sim_type == 2
        if shock_choice == 1
            parameter_sigma_psi = 0.00001;
        elseif shock_choice == 4
            parameter_sigma_a = 0.0060626464427887;
        end
    end
    save('../opts.mat','-append','parameter_sigma_a','parameter_sigma_psi','parameter_rhoA','parameter_Theta')
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
    try
        if models_to_run(ii) == 2
            eval(strcat('dynareOBC gk.mod ',char(opts.dynareOBC_options(jj,:)),';'));
            binding_periods = 1; 
            lev = 1;
        elseif models_to_run(ii) == 7
            eval(strcat('dynareOBC gkq.mod ',char(opts.dynareOBC_options(jj,:)),';'));
            binding_periods = 1;
            lev = ( (oo_.endo_simul(strmatch('Thetax',M_.endo_names,'exact'),:))-((oo_.endo_simul(strmatch('mus',M_.endo_names,'exact'),:))+(oo_.endo_simul(strmatch('mue',M_.endo_names,'exact'),:)).*(oo_.endo_simul(strmatch('E',M_.endo_names,'exact'),:)).*(oo_.endo_simul(strmatch('QE',M_.endo_names,'exact'),:))./(exp((oo_.endo_simul(strmatch('q',M_.endo_names,'exact'),:))).*exp((oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),:))))) ) ./ (oo_.endo_simul(strmatch('nub',M_.endo_names,'exact'),:));
            lev2 = lev + (oo_.endo_simul(strmatch('E',M_.endo_names,'exact'),:)).*(oo_.endo_simul(strmatch('QE',M_.endo_names,'exact'),:))./(exp((oo_.endo_simul(strmatch('q',M_.endo_names,'exact'),:))).*exp((oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),:))));
        elseif models_to_run(ii) == 6
            eval(strcat('dynareOBC obc.mod ',char(opts.dynareOBC_options(jj,:)),';'));
            mv = (oo_.endo_simul(strmatch('mv',M_.endo_names,'exact'),:));
            mv(mv<1e-8) = 0;
            mv(mv>1e-8) = 1;
            binding_periods = mean(mv);
            lev = 1 - exp((oo_.endo_simul(strmatch('b',M_.endo_names,'exact'),:)))./(exp((oo_.endo_simul(strmatch('q',M_.endo_names,'exact'),:))).*exp((oo_.endo_simul(strmatch('k',M_.endo_names,'exact'),:))));
        elseif models_to_run(ii) == 1
            eval(strcat('dynareOBC rbc.mod ',char(opts.dynareOBC_options(jj,:)),';'));
            binding_periods = 0;
        else
            disp('Valid model not chosen!')
        end
        if sim_type==1
        Y = (oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:))./(mean(oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:)));
        I = (oo_.endo_simul(strmatch('inv',M_.endo_names,'exact'),:));
        [~,Y] = hpfilter(Y,1600);
        spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));
        [y_ac,~,~] = autocorr(Y,1);
        disp(horzcat('S.D. Y = ',num2str(std(Y)),'| Target = 0.010563'));
        disp(horzcat('AC(1) Y = ',num2str(y_ac(2)),'| Target = 0.86255'));
        disp(horzcat('mean spread = ',num2str(mean(spread)),'| Target = 0.0057369'));
        disp(horzcat('S.D. spread = ',num2str(std(spread)),'| Target = 0.0017812    '));
        disp(horzcat('Investment skewness: ',num2str(skewness(I))));
        disp(horzcat('Spread skewness: ',num2str(skewness(spread))));
        disp(horzcat('mean bank leverage N/(Q*K): ',num2str(mean(lev))));
        if strcmp(dynareOBC_.BaseFileName,'gkq')
        disp(horzcat('mean bank leverage (N+QE*E)/(Q*K): ',num2str(mean(lev2))));
        end
        disp(horzcat('Constraint binding in ',num2str(100*binding_periods),'% of periods'));
        end
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

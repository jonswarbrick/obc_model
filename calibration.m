clear; close all;

if exist('model/calibrate_log_rbc_order1.txt','file') == 2
delete('model/calibrate_log_rbc_order1.txt')
end
if exist('model/calibrate_log_rbc_order2.txt','file') == 2
delete('model/calibrate_log_rbc_order2.txt')
end
if exist('model/calibrate_log_rbc_order3.txt','file') == 2
delete('model/calibrate_log_rbc_order3.txt')
end
if exist('model/calibrate_log_gk_order1.txt','file') == 2
delete('model/calibrate_log_gk_order1.txt')
end
if exist('model/calibrate_log_gk_order2.txt','file') == 2
delete('model/calibrate_log_gk_order2.txt')
end
if exist('model/calibrate_log_gk_order3.txt','file') == 2
delete('model/calibrate_log_gk_order3.txt')
end
if exist('model/calibrate_log_obc_order1.txt','file') == 2
delete('model/calibrate_log_obc_order1.txt')
end
if exist('model/calibrate_log_obc_order2.txt','file') == 2
delete('model/calibrate_log_obc_order2.txt')
end
if exist('model/calibrate_log_obc_order3.txt','file') == 2
delete('model/calibrate_log_obc_order3.txt')
end

opts.dynareOBC_sim_options = ' NoCubature MILPSolver=cplex OrderOverride=2 FirstOrderConditionalCovariance TimeToEscapeBounds=60 TimeToReturnToSteadyState=20 omega=10000 CompileSimulationCode Sparse';

sim_type = 1;
% 1 = rbc, 2 = gkq, 3 = obc, 4 = nk, 5 = nkobc, 6 = newobc , 7 = gkq
models_to_run = [ 7 ];
num_models = length(models_to_run);
% 1 = non-separable, 2 = additive type 1 , 3 = additive type 2 , 4 =
% non-separable habits on bundles , 5 J-R
utility_type = 5;
% 1 = KQ, 2 = delta
shock_choice = 1;
% 1 = CEE, 2 = Ireland (2003)
adj_type = 1;
% Other
opts.order = 2; % approximation order

[~,numberModels] = size(models_to_run);
if sim_type == 1
opts.dynareOBC_options = opts.dynareOBC_sim_options;
elseif sim_type == 2
opts.dynareOBC_options = opts.dynareOBC_irf_options;
end
opts.models = {'rbc' ; 'gk' ; 'oldobc' ; 'nk' ; 'nkobc' ; 'obc' ; 'gkq'};
opts.sim_type_name = {'sim' ; 'irf'};

save('opts.mat')

cd('model');

for current_model=1:num_models
    save('loop.mat','current_model')
    fid_mod = fopen( 'which_model.mod', 'wt' );
    fprintf( fid_mod, '@#define model_type = %d\n', models_to_run(current_model));
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

    calibrate_parameters;   
    
    load('loop.mat')
 
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


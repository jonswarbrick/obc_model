clear; close all;

if exist('model/calibrate_log_rbc.txt','file') == 2
delete('model/calibrate_log_rbc.txt')
end
if exist('model/calibrate_log_gkq.txt','file') == 2
delete('model/calibrate_log_gkq.txt')
end
if exist('model/calibrate_log_obc.txt','file') == 2
delete('model/calibrate_log_obc.txt')
end

if exist('calibrations')~=7
    mkdir('calibrations')
end

opts.dynareOBC_options = ' OrderOverride=2 TimeToEscapeBounds=60 TimeToReturnToSteadyState=20 omega=10000 CompileSimulationCode Sparse MLVSimulationMode=1';

% 1 = rbc, 2 = gkq, 3 = obc
models_to_run = [ 1 , 2 , 3 ];

num_models = length(models_to_run);
[~,numberModels] = size(models_to_run);
opts.models = {'rbc' ; 'gkq' ; 'obc'};
save('opts.mat')

cd('model');

for current_model=1:num_models
    save('loop.mat','current_model')
    
    
    % Initialization
    options.num_par = 3;
    % Calibration Detail %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sigma_a --> s.d.(Y)                                           %
    % rho_a --> AR(1) of Y                                          %
    % Theta --> s.d.(spread)                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    options.par = {'sigma_a' ; 'rho_a' ; 'Theta'}; 

    if models_to_run(current_model)==1
        options.init_val(1) = 0.006145871647583;
        options.init_val(2) = 0.95;
        options.init_val(3) = 0.9;
    elseif models_to_run(current_model)==2
        options.init_val(1) = 0.005710743233638;
        options.init_val(2) = 0.95;
        options.init_val(3) = 0.893717341213108;
    elseif models_to_run(current_model)==3
        options.init_val(1) = 0.0060626464427887;
        options.init_val(2) = 0.95;
        options.init_val(3) = 0.669734151867773;
    end

    options.target(1) = 0.010563;
    options.target(2) = 0.86255;
    options.target(3) = 0.0017812;

    options.err_tol(1) = 0.000005;
    options.err_tol(2) = 0.15;
    options.err_tol(3) = 0.00002;

    options.init_err_tol = options.err_tol;
    options.curr_val = options.init_val;
    options.prev_val = options.curr_val;

    parameter_sigma_a = options.curr_val(1);
    parameter_rhoA = options.curr_val(2);
    parameter_Theta = options.curr_val(3);

    save('../opts.mat','parameter_sigma_a','parameter_rhoA','parameter_Theta','-append');
    save('settings_file.mat','options');

    %% Initial run
    eval(horzcat('dynareOBC ',char(opts.models(models_to_run(current_model))),'_cal ',char(opts.dynareOBC_options(1,:)),';'));
    Y = dynareOBC_.MLVSimulationWithBounds.y;
    [~,Y] = hpfilter(Y,1600);
    spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));

    load('settings_file.mat');
    load('../opts.mat');
    load('loop.mat');
    [y_ac,~,~] = autocorr(Y,1);
    
    options.prev_realised(1) = std(Y);
    options.prev_realised(2) = y_ac(2);
    options.prev_realised(3) = std(spread);
    options.curr_err(1) = std(Y)-options.target(1);
    options.curr_err(2) = y_ac(2)-options.target(2);
    options.curr_err(3) = std(spread)-options.target(3);
    options.prev_err = options.curr_err;
                if models_to_run(current_model)==1
                    options.curr_err(3)=0;
                end

                err_mess = ['   Errors: ' , num2str(options.curr_err(1)),' ', num2str(options.curr_err(2)),' ', num2str(options.curr_err(3)),' '];
                values_mess = ['   sd(y) = ',num2str(options.prev_realised(1)),' | AC1(y) = ',num2str(options.prev_realised(2)),...
                    ' | sd(spread) = ',num2str(options.prev_realised(3)),];
                calibs_mess = ['   sigma_A = ',num2str(options.curr_val(1)),' | rhoA = ',num2str(options.curr_val(2)),...
                    ' | Theta = ',num2str(options.curr_val(3)),];
                if models_to_run(current_model)==1
                        fid_log = fopen( 'calibrate_log_rbc.txt', 'At' );
                elseif models_to_run(current_model)==2
                        fid_log = fopen( 'calibrate_log_gkq.txt', 'At' );
                elseif models_to_run(current_model)==3
                        fid_log = fopen( 'calibrate_log_obc.txt', 'At' );
                end
                log_txt = strcat(err_mess,'\n',values_mess,'\n',calibs_mess,'\n');
                fprintf( fid_log, log_txt );
                fclose(fid_log);

    %% Calibration
        big_crit = 1;
        options.jj = 1;
    while big_crit
        for ii=1:options.num_par
            options.ii=ii;
            options.iter = 1;
            crit = 1;
            while crit>options.err_tol(ii)
                if abs(options.curr_err(ii))<options.err_tol(ii)
                    break
                end
                if options.iter>1
                    options.step(ii) = options.step(ii) * (  options.prev_err(ii)/(options.prev_err(ii) - options.curr_err(ii))  - 1);
                else
                    options.step(1) = 0.001;
                    options.step(2) = 0.02;
                    options.step(3) = 0.05;
                end
                options.curr_val(ii) = max(0.00001,options.curr_val(ii)+options.step(ii));
                options.curr_val(ii) = min(options.curr_val(ii),.999);
                options.curr_val(2) = min(options.curr_val(2),.98);

                parameter_sigma_a = options.curr_val(1);
                parameter_rhoA = options.curr_val(2);
                parameter_Theta = options.curr_val(3);

                save('../opts.mat','parameter_sigma_a','parameter_rhoA','parameter_Theta','-append');
                save('settings_file.mat','options');
                eval(horzcat('dynareOBC ',char(opts.models(models_to_run(current_model))),'_cal ',char(opts.dynareOBC_options(1,:)),';'));
                Y = dynareOBC_.MLVSimulationWithBounds.y;
                [~,Y] = hpfilter(Y,1600);
                spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));

                load('settings_file.mat');
                load('../opts.mat');
                load('loop.mat');
                ii=options.ii;
                [y_ac,~,~] = autocorr(Y,1);

                options.prev_err(ii) = options.curr_err(ii);
                temp_prev_realised = options.prev_realised;
                options.prev_realised(1) = std(Y);
                options.prev_realised(2) = y_ac(2);
                options.prev_realised(3) = std(spread);
                options.curr_err(1) = std(Y)-options.target(1);
                options.curr_err(2) = y_ac(2)-options.target(2);
                options.curr_err(3) = std(spread)-options.target(3);
                if models_to_run(current_model)==1
                    options.curr_err(3)=0;
                end

                cal_mess = ['>> ',char(opts.models(models_to_run(current_model))),' | ',options.par{ii},' calibration | Big loop ',num2str(options.jj),...
                    ' | Iteration ',num2str(options.iter),' | Step = ',num2str(options.step(ii)),' | ',...
                    options.par{ii},' = ',num2str(options.curr_val(ii),'%.16f'),' | Target = ',num2str(options.target(ii)),' | Realised ',num2str(options.prev_realised(ii)),...
                    ' | Last error = ',num2str(options.prev_err(ii)),' | New error = ',num2str(options.curr_err(ii)),...
                    ' | AC1(y) = ',num2str(options.prev_realised(2))];
                err_mess = ['   Errors: ' , num2str(options.curr_err(1)),' ', num2str(options.curr_err(2)),' ', num2str(options.curr_err(3)),' '];
                disp(cal_mess);
                values_mess = ['   sd(y) = ',num2str(options.prev_realised(1)),' | AC1(y) = ',num2str(options.prev_realised(2)),...
                    ' | mean(spread) = ',num2str(options.prev_realised(3)),];
                calibs_mess = ['   sigma_A = ',num2str(options.curr_val(1)),' | rhoA = ',num2str(options.curr_val(2)),...
                    ' | Theta = ',num2str(options.curr_val(3)),];
                if models_to_run(current_model)==1
                    fid_log = fopen( 'calibrate_log_rbc.txt', 'At' );
                elseif models_to_run(current_model)==2
                    fid_log = fopen( 'calibrate_log_gkq.txt', 'At' );
                elseif models_to_run(current_model)==3
                    fid_log = fopen( 'calibrate_log_obc.txt', 'At' );
                end
                log_txt = strcat(cal_mess,'\n',err_mess,'\n',values_mess,'\n',calibs_mess,'\n');
                fprintf( fid_log, log_txt );
                fclose(fid_log);
                crit = abs(options.curr_err(ii));
                if abs(options.prev_realised-temp_prev_realised)<0.00001
                    crit=0;
                end
                options.iter = options.iter+1;
            end
        end
        if mean(abs(options.curr_err)<10*options.err_tol)==1
            big_crit = 0;
        else
            big_crit = 1;
        end
        options.jj = options.jj+1;
    end
        if models_to_run(current_model)==1
            save( '../calibrations/calibration_rbc.mat' );
        elseif models_to_run(current_model)==2
            save( '../calibrations/calibration_gkq.mat' );
        elseif models_to_run(current_model)==3
            save( '../calibrations/calibration_obc.mat' );
        end

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


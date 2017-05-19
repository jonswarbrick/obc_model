
% Initialization
options.num_par = 3;

% Calibration Detail %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma_a --> s.d.(Y)                                           %
% rho_a --> AR(1) of Y                                          %
% Theta --> s.d.(spread)                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options.par = {'sigma_a' ; 'rho_a' ; 'Theta'}; % order important

if models_to_run(current_model)==1
    options.init_val(1) = 0.0059038615167123 ;
    options.init_val(2) = 0.95;
    options.init_val(3) = 0.9;
elseif models_to_run(current_model)==2
    options.init_val(1) = 0.0068513695401358 ;
    options.init_val(2) = 0.95;
    options.init_val(3) = 0.585227706178260;
elseif models_to_run(current_model)==6
    options.init_val(1) = 0.003500638344089;
    options.init_val(2) = 0.828157065533977;
    options.init_val(3) = 0.569692624538134;
end
    

options.target(1) = 0.010563;
options.target(2) = 0.86255;
options.target(3) = 0.0017812;

options.err_tol(1) = 0.000005;
options.err_tol(2) = 0.15;
options.err_tol(3) = 0.00002;

%options.err_tol(1) = 1;
%options.err_tol(2) = 1;
%options.err_tol(3) = 1;

options.step(1) = 0.001;
options.step(2) = 0.02;
options.step(3) = 0.05;
%options.step = zeros(1,4);

options.init_err_tol = options.err_tol;
options.curr_val = options.init_val;
options.prev_val = options.curr_val;

parameter_sigma_a = options.curr_val(1);
parameter_rhoA = options.curr_val(2);
parameter_Theta = options.curr_val(3);

parameter_sigma_psi = 0;

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
   
save('../opts.mat','parameter_sigma_a','parameter_rhoA','parameter_Theta','parameter_sigma_psi','parameter_Phi','parameter_kappa',...
    'parameter_kappa_new','parameter_nubar','parameter_habits_C','parameter_habits_H','parameter_sigma_g','parameter_sigma_delta',...
    'parameter_rhoG','parameter_rhodelta','parameter_gamR','parameter_gamPi','parameter_gamY','parameter_sigma_h','parameter_psi_h',...
    'parameter_sigma_c','parameter_gam_jr','parameter_theta_jr','parameter_rho_psi','-append');
save('settings_file.mat','options');

%% Initial run

eval(horzcat('dynareOBC ',char(opts.models(models_to_run(current_model))),' ',char(opts.dynareOBC_options(3,:)),';'));
                 
Y = (oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:))./(mean(oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:)));
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
                    fid_log = fopen( 'calibrate_log_rbc_order3.txt', 'At' );
            elseif models_to_run(current_model)==2
                    fid_log = fopen( 'calibrate_log_gk_order3.txt', 'At' );
            elseif models_to_run(current_model)==6
                    fid_log = fopen( 'calibrate_log_obc_order3.txt', 'At' );
            end
            log_txt = strcat(err_mess,'\n',values_mess,'\n',calibs_mess,'\n');
            fprintf( fid_log, log_txt );
            fclose(fid_log);

%% Calibration

for or=1:3;
    options.or=or;
    options.or=3;
    big_crit = 1;
   % options.curr_err = 1;
   % options.prev_err = options.curr_err;
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
%             if options.iter>1
%             options.temp_val = options.curr_val(ii); 
%             options.curr_val(ii) = options.prev_val(ii)*( options.target(ii) / options.prev_realised(ii) );
%             options.step(ii) = options.curr_val(ii)-options.prev_val(ii);
%             options.prev_val(ii) = options.temp_val;
%             end
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
            eval(horzcat('dynareOBC ',char(opts.models(models_to_run(current_model))),' ',char(opts.dynareOBC_options(options.or,:)),';'));
                 
            Y = (oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:))./(mean(oo_.endo_simul(strmatch('Y',M_.endo_names,'exact'),:)));
            [~,Y] = hpfilter(Y,1600);
            spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));

            load('settings_file.mat');
            load('../opts.mat');
            load('loop.mat');
            ii=options.ii;
            or=options.or;

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
                if options.or==1
                    fid_log = fopen( 'calibrate_log_rbc_order1.txt', 'At' );
                elseif options.or==2
                    fid_log = fopen( 'calibrate_log_rbc_order2.txt', 'At' );
                elseif options.or==3
                    fid_log = fopen( 'calibrate_log_rbc_order3.txt', 'At' );
                end
            elseif models_to_run(current_model)==2
                if options.or==1
                    fid_log = fopen( 'calibrate_log_gk_order1.txt', 'At' );
                elseif options.or==2
                    fid_log = fopen( 'calibrate_log_gk_order2.txt', 'At' );
                elseif options.or==3
                    fid_log = fopen( 'calibrate_log_gk_order3.txt', 'At' );
                end
            elseif models_to_run(current_model)==6
                if options.or==1
                    fid_log = fopen( 'calibrate_log_obc_order1.txt', 'At' );
                elseif options.or==2
                    fid_log = fopen( 'calibrate_log_obc_order2.txt', 'At' );
                elseif options.or==3
                    fid_log = fopen( 'calibrate_log_obc_order3.txt', 'At' );
                end
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
    %big_crit = sum(abs(options.curr_err));
    options.jj = options.jj+1;
end
    if models_to_run(current_model)==1
        if options.or==1
            save( '../calibrations/calibration_rbc_order1.mat' );
        elseif options.or==2
            save( '../calibrations/calibration_rbc_order2.mat' );
        elseif options.or==3
            save( '../calibrations/calibration_rbc_order3.mat' );
        end
    elseif models_to_run(current_model)==2
        if options.or==1
            save( '../calibrations/calibration_gk_order1.mat' );
        elseif options.or==2
            save( '../calibrations/calibration_gk_order2.mat' );
        elseif options.or==3
            save( '../calibrations/calibration_gk_order3.mat' );
        end
    elseif models_to_run(current_model)==6
        if options.or==1
            save( '../calibrations/calibration_obc_order1.mat' );
        elseif options.or==2
            save( '../calibrations/calibration_obc_order2.mat' );
        elseif options.or==3
            save( '../calibrations/calibration_obc_order3.mat' );
        end
    end
end

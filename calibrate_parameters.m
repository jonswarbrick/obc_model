clear; close all;
    if exist('calibrate_log_order1.txt','file') == 2
    delete('calibrate_log_order1.txt')
    end
    if exist('calibrate_log_order2.txt','file') == 2
    delete('calibrate_log_order2.txt')
    end
    if exist('calibrate_log_order3.txt','file') == 2
    delete('calibrate_log_order3.txt')
    end

% Initialization
options.num_par = 7;

% Calibration Detail %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma_a --> s.d.(Y)                                           %
% rho_a --> AR(1) of Y                                          %
% Theta --> s.d.(spread)                                        %
% sigma_psi --> s.d.(investment)                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options.par = {'sigma_a' ; 'rho_a' ; 'sigma_risk' ;...
                         'sigma_val' ; 'lambda' ; 'gam_bar' ; 'hatw'}; % order important

options.init_val(1) = 0.002094565089383;
options.init_val(2) = 0.933419530006256;
options.init_val(3) = 0.059919214735920;
options.init_val(4) = 0.004506270901719;
options.init_val(5) = 0.809100000000000;
options.init_val(6) = 0.952212813122282;
options.init_val(7) = 15.8;

options.target(1) = 0.0101462755522287;
options.target(2) = 0.914496282145410;
options.target(3) = -0.141444803174180;
options.target(4) = -0.676263248414584;
options.target(5) = 0.028;
options.target(6) = 0.0099;
options.target(7) = 15.8;

options.step(1) = 0.0005;
options.step(2) = 0.01;
options.step(3) = 0.005;
options.step(4) = 0.001;
options.step(5) = 0.005;
options.step(6) = 0.005;
options.step(7) = 1;

options.err_tol(1) = 0.000001;
options.err_tol(2) = 0.001;
options.err_tol(3) = 0.001;
options.err_tol(4) = 0.001;
options.err_tol(5) = 0.001;
options.err_tol(6) = 0.00001;
options.err_tol(7) = 0.01;

options.init_err_tol = options.err_tol;

for ii=1:options.num_par
options.curr_val(ii) = options.init_val(ii);
end


%% Initial run

for init_run=1:3

param_sigma_a = options.curr_val(1);
param_rho_a = options.curr_val(2);
param_sigma_risk = options.curr_val(3);
param_sigma_val = options.curr_val(4);
param_lambda = options.curr_val(5);
param_gam = options.curr_val(6);
param_hatw = options.curr_val(7);
param_rho_risk = 0.7;
param_rho_eff = 0.7;
param_MC_solver = 1;

options.init_run = init_run;

save('param_sizes.mat','param_sigma_a','param_sigma_val','param_sigma_risk','param_rho_a','param_rho_risk','param_rho_eff','param_gam','param_lambda','param_hatw','param_MC_solver');
save('settings_file.mat','options');

if init_run==1
    dynareOBC adverse_sim LPSolver=cplex OrderOverride=1 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
elseif init_run==2
    dynareOBC adverse_sim LPSolver=cplex OrderOverride=2 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
elseif init_run==3
    dynareOBC adverse_sim LPSolver=cplex OrderOverride=3 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
end

Y = exp(oo_.endo_simul(strmatch('y_a',M_.endo_names,'exact'),:))./mean(exp(oo_.endo_simul(strmatch('y_a',M_.endo_names,'exact'),:)));
I = exp(oo_.endo_simul(strmatch('i_a',M_.endo_names,'exact'),:))./mean(exp(oo_.endo_simul(strmatch('i_a',M_.endo_names,'exact'),:)));
pr = oo_.endo_simul(strmatch('pr_a',M_.endo_names,'exact'),:);
gam = oo_.endo_simul(strmatch('gam',M_.endo_names,'exact'),:);
u = oo_.endo_simul(strmatch('u_a',M_.endo_names,'exact'),:);
e = oo_.endo_simul(strmatch('e_a',M_.endo_names,'exact'),:);
a = oo_.endo_simul(strmatch('a',M_.endo_names,'exact'),:);
xs = oo_.endo_simul(strmatch('xs',M_.endo_names,'exact'),:);
xr = oo_.endo_simul(strmatch('xr',M_.endo_names,'exact'),:);
val = oo_.endo_simul(strmatch('val',M_.endo_names,'exact'),:);
pr_bar = exp(logit_pr_bar)/(exp(logit_pr_bar)+1);
Rs_bar = (param_gam / betta -  param_gam + param_lambda*xs_bar +  (1-param_lambda) )  / (param_lambda*xs_bar + (1-param_lambda)*(pr_bar*omegar_bar - pr_bar*( omegar_bar -1 )* xs_bar ) );
W_bar = (1 - alp)*( (( alp / (Rs_bar - 1 + delta) )^(1/(1-alp))) )^alp;
hhat_bar = (W_bar/chi) / ((param_lambda*xs_bar + (1- param_lambda) )*(( alp / (Rs_bar - 1 + delta) )^(1/(1-alp)))^(alp-1)*(epsilon + kappa)/param_hatw  + ( (param_lambda*xs_bar + (1-param_lambda) )*W_bar/chi -  (param_lambda*xs_bar + (1-param_lambda) )*delta*(epsilon+kappa))/param_hatw +  W_bar/chi );
e_bar = hhat_bar/param_hatw;

load('settings_file.mat');
load('param_sizes.mat');

init_run = options.init_run;

hhat = ones(1,1000);
khat = ones(1,1000);
default_rate = ones(1,1000);
spread = ones(1,1000);
hatw = ones(1,1000);
hatw(1) = (1 - e_bar - u(1))/e_bar;
khat(1) = ( param_lambda*xs_bar + (1-param_lambda))*( kappa +  epsilon )*e_bar;
default_rate(1) =  (1-param_lambda).*(1-pr(1))./param_gam;
spread(1) = ( alp * ( hhat(1) / khat(1) )^(1-alp) + (1-delta) ) * ( exp(val(1))/pr(1) - 1 ) * ( 1 -  xs_bar );
hatw(1) = hhat(1)/(hhat_bar/param_hatw);

for t=2:1000
    khat(t) = ( param_lambda*xs(t-1) + (1-param_lambda)*xr(t-1)*pr(t)*exp(val(t))/pr(t))*( kappa +  epsilon )*e(t-1);
    hhat(t) = 1-e(t-1)-u(t);
    default_rate(t) = (1-param_lambda).*(1-pr(t))./gam(t-1);
    spread(t) = ( alp * exp(a(t)) * ( hhat(t) / khat(t) )^(1-alp) + (1-delta) ) * ( exp(val(t))/pr(t) - 1 ) * ( 1 -  xs(t-1) / xr(t-1) );
    hatw(t) = (1 - e(t-1) - u(t))/e(t-1);
end

[y_ac,~,~] = autocorr(Y,1);
init_err(1) = std(Y)-options.target(1);
init_err(2) = y_ac(2)-options.target(2);
init_err(3) = skewness(Y)-options.target(3);
init_err(4) = skewness(I)-options.target(4);
init_err(5) = mean(default_rate)-options.target(5);
init_err(6) = mean(spread)-options.target(6);
init_err(7) = mean(hatw)-options.target(7);

if init_run==1
    options.init_err_order1 = init_err;
elseif init_run==2
    options.init_err_order2 = init_err;
elseif init_run==3
    options.init_err_order3 = init_err;
end

end

%% Calibration

for or=1:3;
    options.or=or;
    big_crit = 1;
    if or==1
        options.curr_err = options.init_err_order1;
    elseif or==2
        options.curr_err = options.init_err_order2;
    elseif or==3
        options.curr_err = options.init_err_order3;
    end
    options.prev_err = options.curr_err;
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
            end
            options.curr_val(ii) = options.curr_val(ii)+options.step(ii);
            options.curr_val(2) = min( options.curr_val(2) , 0.99 );

            param_sigma_a = options.curr_val(1);
            param_rho_a = options.curr_val(2);
            param_sigma_risk = options.curr_val(3);
            param_sigma_val = options.curr_val(4);
            param_lambda = options.curr_val(5);
            param_gam = options.curr_val(6);
            param_hatw = options.curr_val(7);
            param_rho_risk = 0.7;
            param_rho_eff = 0.7;
            param_MC_solver = 1;

            save('param_sizes.mat','param_sigma_a','param_sigma_val','param_sigma_risk','param_rho_a','param_rho_risk','param_rho_eff','param_gam','param_lambda','param_hatw','param_MC_solver');
            save('settings_file.mat','options');
            if options.or==1
                dynareOBC adverse_sim LPSolver=cplex OrderOverride=1 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
            elseif options.or==2
                dynareOBC adverse_sim LPSolver=cplex OrderOverride=2 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
            elseif options.or==3
                dynareOBC adverse_sim LPSolver=cplex OrderOverride=3 TimeToEscapeBounds=32 NoCubature firstorderconditionalcovariance Omega=10000 Sparse MLVSimulationMode=0 CompileSimulationCode
            end
            Y = exp(oo_.endo_simul(strmatch('y_a',M_.endo_names,'exact'),:))./mean(exp(oo_.endo_simul(strmatch('y_a',M_.endo_names,'exact'),:)));
            I = exp(oo_.endo_simul(strmatch('i_a',M_.endo_names,'exact'),:))./mean(exp(oo_.endo_simul(strmatch('i_a',M_.endo_names,'exact'),:)));
            pr = oo_.endo_simul(strmatch('pr_a',M_.endo_names,'exact'),:);
            gam = oo_.endo_simul(strmatch('gam',M_.endo_names,'exact'),:);
            u = oo_.endo_simul(strmatch('u_a',M_.endo_names,'exact'),:);
            e = oo_.endo_simul(strmatch('e_a',M_.endo_names,'exact'),:);
            a = oo_.endo_simul(strmatch('a',M_.endo_names,'exact'),:);
            xs = oo_.endo_simul(strmatch('xs',M_.endo_names,'exact'),:);
            xr = oo_.endo_simul(strmatch('xr',M_.endo_names,'exact'),:);
            val = oo_.endo_simul(strmatch('val',M_.endo_names,'exact'),:);
            pr_bar = exp(logit_pr_bar)/(exp(logit_pr_bar)+1);
            Rs_bar = (param_gam / betta -  param_gam + param_lambda*xs_bar +  (1-param_lambda) )  / (param_lambda*xs_bar + (1-param_lambda)*(pr_bar*omegar_bar - pr_bar*( omegar_bar -1 )* xs_bar ) );
            W_bar = (1 - alp)*( (( alp / (Rs_bar - 1 + delta) )^(1/(1-alp))) )^alp;
            hhat_bar = (W_bar/chi) / ((param_lambda*xs_bar + (1- param_lambda) )*(( alp / (Rs_bar - 1 + delta) )^(1/(1-alp)))^(alp-1)*(epsilon + kappa)/param_hatw  + ( (param_lambda*xs_bar + (1-param_lambda) )*W_bar/chi -  (param_lambda*xs_bar + (1-param_lambda) )*delta*(epsilon+kappa))/param_hatw +  W_bar/chi );
            e_bar = hhat_bar/param_hatw;

            load('settings_file.mat');
            load('param_sizes.mat');
            ii=options.ii;
            or=options.or;

            hhat = ones(1,1000);
            khat = ones(1,1000);
            default_rate = ones(1,1000);
            spread = ones(1,1000);
            hatw = ones(1,1000);

            hatw(1) = (1 - e_bar - u(1))/e_bar;
            khat(1) = ( param_lambda*xs_bar + (1-param_lambda))*( kappa +  epsilon )*e_bar;
            default_rate(1) =  (1-param_lambda).*(1-pr(1))./param_gam;
            spread(1) = ( alp * ( hhat(1) / khat(1) )^(1-alp) + (1-delta) ) * ( exp(val(1))/pr(1) - 1 ) * ( 1 -  xs_bar );
            hatw(1) = hhat(1)/(hhat_bar/param_hatw);

            for t=2:1000
                khat(t) = ( param_lambda*xs(t-1) + (1-param_lambda)*xr(t-1)*pr(t)*exp(val(t))/pr(t))*( kappa +  epsilon )*e(t-1);
                hhat(t) = 1-e(t-1)-u(t);
                default_rate(t) = (1-param_lambda).*(1-pr(t))./gam(t-1);
                spread(t) = ( alp * exp(a(t)) * ( hhat(t) / khat(t) )^(1-alp) + (1-delta) ) * ( exp(val(t))/pr(t) - 1 ) * ( 1 -  xs(t-1) / xr(t-1) );
                hatw(t) = (1 - e(t-1) - u(t))/e(t-1);
            end

            [y_ac,~,~] = autocorr(Y,1);

            options.prev_err(ii) = options.curr_err(ii);
            options.curr_err(1) = std(Y)-options.target(1);
            options.curr_err(2) = y_ac(2)-options.target(2);
            options.curr_err(3) = skewness(Y)-options.target(3);
            options.curr_err(4) = skewness(I)-options.target(4);
            options.curr_err(5) = mean(default_rate)-options.target(5);
            options.curr_err(6) = mean(spread)-options.target(6);
            options.curr_err(7) = mean(hatw)-options.target(7);


            cal_mess = [options.par{ii},' calibration | Big loop ',num2str(options.jj),' | Iteration ',num2str(options.iter),' | Step = ',num2str(options.step(ii)),' | ',...
                options.par{ii},' = ',num2str(options.curr_val(ii),'%.16f'),' | Last error = ',num2str(options.prev_err(ii)),' | New error = ',num2str(options.curr_err(ii))];
            err_mess = ['Errors: ' , num2str(options.curr_err(1)),' ', num2str(options.curr_err(2)),' ', num2str(options.curr_err(3)),' ', num2str(options.curr_err(4)),' '...
                , num2str(options.curr_err(5)),' ', num2str(options.curr_err(6)),' ', num2str(options.curr_err(7)),' '];
            disp(cal_mess);
            if options.or==1
                fid_log = fopen( 'calibrate_log_order1.txt', 'At' );
            elseif options.or==2
                fid_log = fopen( 'calibrate_log_order2.txt', 'At' );
            elseif options.or==3
                fid_log = fopen( 'calibrate_log_order3.txt', 'At' );
            end
            log_txt = strcat(cal_mess,'\n',err_mess,'\n');
            fprintf( fid_log, log_txt );
            fclose(fid_log);
            crit = abs(options.curr_err(ii));
            options.iter = options.iter+1;
        end
    end
    %options.err_tol = 01*options.init_err_tol;
    if mean(abs(options.curr_err)<10*options.err_tol)==1
        big_crit = 0;
    else
        big_crit = 1;
    end
    %big_crit = sum(abs(options.curr_err));
    options.jj = options.jj+1;
end
    if options.or==1
        save( 'calibration_order1.mat' );
    elseif options.or==2
        save( 'calibration_order2.mat' );
    elseif options.or==3
        save( 'calibration_order3.mat' );
    end
end

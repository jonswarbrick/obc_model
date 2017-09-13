clear;close all;


%% Prepare Data For Plots -- User Settings
opt.epsshock = char('eps_psi');
opt.shock_shortname = char('KQ');
opt.variables = char('Y','inv','H','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('relAroundZero','irfsAroundZero','relAroundZero','irfsAroundZero','irfs','irfs');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 1;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'newobc_order2_median_KQ5pc'...
    ,'rbc_order2_median_KQ5pc'...
    ,'gkq_order2_median_KQ5pc'...
    ...'newobc_order2_slow_KQ5pc'...
    ...,'rbc_order2_QMC_KQ5pc'...
    ...,'gkq_order2_QMC_KQ5pc'...
     };
opt.type_desc = '5pc';
opt.plot_size = [230 250 500 350];
opt.nolegend = 0;
IRF_plotter_paper_new(opt);
clear;
%% Prepare Data For Plots -- User Settings
opt.epsshock = char('eps_psi');
opt.shock_shortname = char('KQ');
opt.variables = char('Y','inv','H','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('relAroundZero','irfsAroundZero','relAroundZero','irfsAroundZero','irfs','irfs');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 2;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'newobc_order2_median_posKQ5pc'...
    ,'rbc_order2_median_posKQ5pc'...
    ,'gkq_order2_median_posKQ5pc'...
    ...'newobc_order2_slow_posKQ5pc'...
    ...,'rbc_order2_QMC_posKQ5pc'...
    ...,'gkq_order2_QMC_posKQ5pc'...
     };
opt.type_desc = 'pos5pc';
opt.plot_size = [230 250 500 350];
opt.nolegend = 0;
IRF_plotter_paper_new(opt);
clear;
%% Prepare Data For Plots -- User Settings
opt.epsshock = char('epsA');
opt.shock_shortname = char('A');
opt.variables = char('Y','inv','H','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('relAroundZero','irfsAroundZero','relAroundZero','irfsAroundZero','irfs','irfs');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 3;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'newobc_order2_median_A1sd'...
    ,'rbc_order2_median_A1sd'...
    ,'gkq_order2_median_A1sd'...
    ...'newobc_order2_slow_A1sd'...
    ...,'rbc_order2_QMC_A1sd'...
    ...,'gkq_order2_QMC_A1sd'...
     };
opt.type_desc = '1sd';
opt.plot_size = [230 250 500 350];
opt.nolegend = 0;
IRF_plotter_paper_new(opt);
clear;

%% Prepare Data For Plots -- User Settings
opt.epsshock = char('epsA');
opt.shock_shortname = char('A');
opt.variables = char('Y','inv','H','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('relAroundZero','irfsAroundZero','relAroundZero','irfsAroundZero','irfs','irfs');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 4;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'newobc_order2_median_negA1sd'...
    ,'rbc_order2_median_negA1sd'...
    ,'gkq_order2_median_negA1sd'...
    ...'newobc_order2_slow_negA1sd'...
    ...,'rbc_order2_QMC_negA1sd'...
    ...,'gkq_order2_QMC_negA1sd'...
     };
opt.type_desc = 'neg1sd';
opt.plot_size = [230 250 500 350];
opt.nolegend = 0;
IRF_plotter_paper_new(opt);
clear;

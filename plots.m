clear;close all;

%% Figure 1
quarterly_baron_data = xlsread('data/baron_bank_data.xlsx',1,'B42:K235');
annual_baron_data = xlsread('data/baron_bank_data.xlsx',1,'B3:K40');
dividend_smoothed = quarterly_baron_data(:,4);
repurchase_smoothed = quarterly_baron_data(:,8);
new_equity_smoothed = quarterly_baron_data(:,10);
new_equity_smoothed_old = [ annual_baron_data(:,10); new_equity_smoothed(1)];
annual_time = [ datetime(1927,10,01) , datetime(1927,10,01) + calyears(1:38)]';
quarterly_time = [ datetime(1965,10,01) , datetime(1965,10,01) + 3*calmonths(1:193)]';

h = figure;
set(h, 'Position', [50 , 50, 500, 300]);
plot(quarterly_time,100*new_equity_smoothed,'color',[0,.4,.8],'LineWidth',1.5);hold on; 
plot(annual_time,100*new_equity_smoothed_old,'color',[0,.4,.8],'LineWidth',1.5);
recessionplot;
ylabel('\% of book equity, annualized','Interpreter','latex')
datetick('x','YYYY')
xmin = datetime(1929,12,03);xmax = datetime(2013,11,01);
xlim([xmin xmax])
ylim([0 50])
set(gca,'TickLabelInterpreter','latex')
img_filename = strcat('plot_files\equity_issuance');
print(h,img_filename,'-r300','-dpng')
print(h,img_filename,'-r300','-depsc')
clear;

%% Figure 3
load('data/rawData.mat');
time = EQTA.EQTA(77:101);
data_1 = 100 * (EQTA.TotalEquitytoTotalAssetsforBanksPercentQuarterlyNotSeasonallyAd(77)./EQTA.TotalEquitytoTotalAssetsforBanksPercentQuarterlyNotSeasonallyAd(77:101) ) - 100;
data_2 = baron.eq_sm4;
data_3 = baron.div_sm4;
load('results_irf/obc_KQ5pc.mat');
obc_1 = [ 0 0 0 0 0 oo_.irfs.lev_eps_psi(1:20) ];
obc_2 = [ 0 0 0 0 0 oo_.irfs.E_rate_eps_psi(1:20) ];
obc_3 = [ 0 0 0 0 0 oo_.irfs.D_rate_eps_psi(1:20) ];
load('results_irf/gkq_KQ5pc.mat');
gkq_1 = [ 0 0 0 0 0 oo_.irfs.lev_eps_psi(1:20) ];
gkq_2 = [ 0 0 0 0 0 oo_.irfs.E_rate_eps_psi(1:20) ];
gkq_3 = [ 0 0 0 0 0 oo_.irfs.D_rate_eps_psi(1:20) ];

h = figure;
plot(time,zeros(1,25),'Color',[.8 .8 .8],'LineWidth',1); hold on;
plot(time,data_1,'Color',[.2 .4 .6],'LineWidth',1.5);
set(gca,'TickLabelInterpreter','latex')
ylabel({'Leverage ratio';'\% deviation from 2007q1'},'Interpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/leverage','-r300','-dpng')
print(h,'plot_files/leverage','-r300','-depsc')

h = figure;
plot(time,100*data_2./4,'Color',[.2 .4 .6],'LineWidth',1.5);
set(gca,'TickLabelInterpreter','latex')
ylabel({'Equity issuance';'(\%)'},'Interpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/equity_issuance','-r300','-dpng')
print(h,'plot_files/equity_issuance','-r300','-depsc')

h = figure;
plot(time,zeros(1,25),'Color',[.8 .8 .8],'LineWidth',1); hold on;
plot(time,100*(data_3-data_3(1))./4,'Color',[.2 .4 .6],'LineWidth',1.5);
set(gca,'TickLabelInterpreter','latex')
ylabel({'Dividend rate';'ppt deviation from 2007q1'},'Interpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/dividends','-r300','-dpng')
print(h,'plot_files/dividends','-r300','-depsc')

h = figure;
yyaxis left
plot(time,zeros(1,25),'Color',[.8 .8 .8],'LineWidth',1); hold on;
plot(time,100*obc_1,'Color',[0 .6 .4],'LineStyle','-','LineWidth',1.5);
ylim([-.5 1]);
set(gca,'YColor',[0 0 0]);
yyaxis right
plot(time,100*gkq_1,'Color',[0.6 0 0],'LineStyle',':','LineWidth',1.5);
ylim([-11 22]);
set(gca,'TickLabelInterpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/leverage_model','-r300','-dpng')
print(h,'plot_files/leverage_model','-r300','-depsc')

h = figure;
yyaxis left
plot(time,zeros(1,25),'Color',[.8 .8 .8],'LineWidth',1); hold on;
plot(time,100*obc_2,'Color',[0 .6 .4],'LineStyle','-','LineWidth',1.5);
ylim([-.5 .5]);
set(gca,'YColor',[0 0 0]);
yyaxis right
plot(time,100*gkq_2,'Color',[0.6 0 0],'LineStyle',':','LineWidth',1.5);
ylim([-3 3]);
set(gca,'TickLabelInterpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/leverage_model','-r300','-dpng')
print(h,'plot_files/leverage_model','-r300','-depsc')

h = figure;
plot(time,zeros(1,25),'Color',[.8 .8 .8],'LineWidth',1); hold on;
plot(time,100*obc_3,'Color',[0 .6 .4],'LineStyle','-','LineWidth',1.5);
plot(time,100*gkq_3,'Color',[0.6 0 0],'LineStyle',':','LineWidth',1.5);
ylim([-0.62 0.08]);
set(gca,'TickLabelInterpreter','latex')
set (h,'Position',[200,200,280,160]);
print(h,'plot_files/leverage_model','-r300','-dpng')
print(h,'plot_files/leverage_model','-r300','-depsc')


%% Figure 4
opt.epsshock = char('eps_psi');
opt.shock_shortname = char('KQ');
opt.variables = char('y','inv','h','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 1;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'obc_KQ5pc'...
    ,'rbc_KQ5pc'...
    ,'gkq_KQ5pc'...
     };
opt.num_periods = 30;
opt.type_desc = '5pc';
opt.plot_size = [230 250 500 300];
opt.nolegend = 0;
plotter(opt);
clear;
%% Figure 5
opt.epsshock = char('epsA');
opt.shock_shortname = char('A');
opt.variables = char('y','inv','h','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 4;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'obc_negA1sd'...
    ,'rbc_negA1sd'...
    ,'gkq_negA1sd'...
     };
opt.num_periods = 30;
opt.type_desc = 'neg1sd';
opt.plot_size = [230 250 500 300];
opt.nolegend = 0;
plotter(opt);
clear;

%% Figure E1
opt.epsshock = char('eps_psi');
opt.shock_shortname = char('KQ');
opt.variables = char('y','inv','h','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 2;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'obc_posKQ5pc'...
    ,'rbc_posKQ5pc'...
    ,'gkq_posKQ5pc'...
     };
opt.num_periods = 30;
opt.type_desc = 'pos5pc';
opt.plot_size = [230 250 500 300];
opt.nolegend = 0;
plotter(opt);
clear;
%% Figure E2
opt.epsshock = char('epsA');
opt.shock_shortname = char('A');
opt.variables = char('y','inv','h','spread','D_rate','E_rate'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero');
opt.names = char('$Y$','$I$','$H$','$\Delta$','$D/N$','$\quad E/N$');
opt.no_rows_sub_plots = 2;
opt.no_cols_sub_plots = 3;
opt.plot_num = 3;

% Choice of data
opt.data_files_desc = {'Our Model','RBC','GK'};
opt.num_datasets = 3;
opt.data_files = {...
    'obc_A1sd'...
    ,'rbc_A1sd'...
    ,'gkq_A1sd'...
     };
opt.num_periods = 30;
opt.type_desc = '1sd';
opt.plot_size = [230 250 500 300];
opt.nolegend = 0;
plotter(opt);
clear;


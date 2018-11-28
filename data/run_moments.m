clear; close all;
warning('off','all')

%% Dividend and Equity moments
load('rawData.mat')
quarterly_baron_data = xlsread('baron_bank_data.xlsx',1,'B42:K235');
dividend_smoothed = quarterly_baron_data(:,4);
repurchase_smoothed = quarterly_baron_data(:,8);
new_equity_smoothed = quarterly_baron_data(:,10);

indexDate = '01/01/01';
pop = fints( freddata.pop.Data(:,1) , freddata.pop.Data(:,2) , 'pop'); 
pop_index = fts2mat(pop(indexDate));
pop = fillts(toquarterly(fints( freddata.pop.Data(:,1) , freddata.pop.Data(:,2)/pop_index , 'pop')),'c');
pop = pop('::07/01/2016');
y = fints(freddata.gdp.Data(:,1),freddata.gdp.Data(:,2),'gdp');
y = y('04/01/1948::');
def = fints(freddata.P_gdp.Data(:,1),freddata.P_gdp.Data(:,2),'gdpdef');
def = def('04/01/1948::');
output = log ( ( fts2mat(y.gdp) ./  fts2mat(def.gdpdef) ) ./ fts2mat(pop) ).*100;
[~,hp.y] = hpfilter(output(142:end),1600);

Y = hp.y(1:end-10);
E = new_equity_smoothed(72:end);
D = dividend_smoothed(72:end)+repurchase_smoothed(72:end);

% Cross correlations
corr_coeff = corrcoef([ Y , D , E ]);
disp('**-------------------------------------------------------------------------------**');
disp('Banking variables correlation coefficients');
disp(table(corr_coeff(:,1),corr_coeff(:,2),corr_coeff(:,3),'VariableNames',{'Y','D','E'},'RowNames',{'Y','D','E'}));

% Moments
moments = [ mean([ D , E ])' std([ D , E ])' skewness([ D , E ])' kurtosis([ D , E ])' ];
disp('**-------------------------- -----------------------------------------------------**');
disp('Banking variables moments');
disp(table(moments(:,1),moments(:,2),moments(:,3),moments(:,4),'VariableNames',{'Mean','StdDev','Skewness','Kurtosis'},'RowNames',{'D','E'}));
disp('**-------------------------------------------------------------------------------**');

%% Y, I, C and spread
load('rawData.mat')
indexDate = '01/01/01';
pop = fints( freddata.pop.Data(:,1) , freddata.pop.Data(:,2) , 'pop'); 
pop_index = fts2mat(pop(indexDate));
pop = fillts(toquarterly(fints( freddata.pop.Data(:,1) , freddata.pop.Data(:,2)/pop_index , 'pop')),'c');
pop = pop('::07/01/2016');
y = fints(freddata.gdp.Data(:,1),freddata.gdp.Data(:,2),'gdp');
y = y('04/01/1948::');
c = fints(freddata.C.Data(:,1),freddata.C.Data(:,2),'c');
c = c('04/01/1948::');
invest = fints(freddata.inv.Data(:,1),freddata.inv.Data(:,2),'invest');
invest = invest('04/01/1948::');
def = fints(freddata.P_gdp.Data(:,1),freddata.P_gdp.Data(:,2),'gdpdef');
def = def('04/01/1948::');

spread = fints(freddata.spread.Data(:,1),freddata.spread.Data(:,2),'spread');
output = log ( ( fts2mat(y.gdp) ./  fts2mat(def.gdpdef) ) ./ fts2mat(pop) ).*100;
cons = log ( ( fts2mat(c.c) ./  fts2mat(def.gdpdef) ) ./ fts2mat(pop) ).*100;
invest = log ( ( fts2mat(invest.invest) ./  fts2mat(def.gdpdef) ) ./ fts2mat(pop) ).*100;

spread = spread('07/01/1983::07/01/2016');
[~,hp.y] = hpfilter(output(142:end),1600);
[~,hp.c] = hpfilter(cons(142:end),1600);
[~,hp.i] = hpfilter(invest(142:end),1600);
S_qtly = 100*((1+0.01*((fts2mat(toquarterly(spread,'CalcMethod','v21x'))))).^(1/4)-1);
S = 100*((1+0.01*((fts2mat(spread)))).^(1/4)-1);

all_data = [ hp.y(1:end-1) , hp.c(1:end-1) , hp.i(1:end-1) , S_qtly(1:end-1)  ]; 

mom.Y = [ mean(hp.y) std(hp.y) skewness(hp.y) kurtosis(hp.y) ];
mom.I = [ mean(hp.i) std(hp.i) skewness(hp.i) kurtosis(hp.i) ];
mom.C = [ mean(hp.c) std(hp.c) skewness(hp.c) kurtosis(hp.c) ];
mom.S = [ mean(S) std(S) skewness(S) kurtosis(S) ];

all_mom = [ mom.Y ; mom.I ; mom.C ; mom.S ];
corr_coeff = corrcoef(all_data);
variables = {'Y','I','C','spread'};

paper_data = [ ...
    corr_coeff(1,1) , corr_coeff(1,2) , corr_coeff(1,3) ,...
    corr_coeff(1,4) ; ...
    ( all_mom(:,2) )' ; ( all_mom(:,3) )' ]; 

tab_paper = table(paper_data(:,1),paper_data(:,2),...
    paper_data(:,3),paper_data(:,4),...
    'VariableNames',variables,'RowNames',{'corr','sd','skew'});

disp('**-------------------------------------------------------------------------------**');
disp('Moments for table 1');
disp(tab_paper);
disp('**-------------------------------------------------------------------------------**');

warning('on','all')

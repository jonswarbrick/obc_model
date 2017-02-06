clear;close all;



variables = {'Y','I','C','D','E','Q','spread'};

% RBC
%load('results_sim/rbc_paper.mat')
load('results_sim/newobc_order3_QMC_default.mat')
spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
rbc.Y =  100*oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
rbc.I =  100*oo_.endo_simul(I_element,:);
rbc.C =  100*oo_.endo_simul(C_element,:);
rbc.Q =  100*oo_.endo_simul(Q_element,:);
rbc.spread =  100*oo_.endo_simul(spread_element,:);
none = zeros(1,length(rbc.Y));

data.rbc = [rbc.Y' rbc.I' rbc.C' none' none' rbc.Q' rbc.spread'];
sd.rbc = std(data.rbc)';
%sd.rbc = sd.rbc./sd.rbc(1);
sk.rbc = skewness(data.rbc)';
corr_coeff = corrcoef(data.rbc);
corrY.rbc = corr_coeff(:,1);
av.rbc = mean(data.rbc)';

% GKQ
clearvars -except sd sk av corrY none variables;

%load('results_sim/gkq_paper.mat')
load('results_sim/gk_order3_QMC_default.mat')
spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
D_element = strmatch('D_rate',M_.endo_names,'exact');
GKQ.Y =  100*oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
GKQ.I =  100*oo_.endo_simul(I_element,:);
GKQ.C =  100*oo_.endo_simul(C_element,:);
GKQ.Q =  100*oo_.endo_simul(Q_element,:);
GKQ.spread =  100*oo_.endo_simul(spread_element,:);
GKQ.D =  100*oo_.endo_simul(D_element,:);
GKQ.E = none;

data.GKQ = [GKQ.Y' GKQ.I' GKQ.C' GKQ.D' GKQ.E' GKQ.Q' GKQ.spread'];
sd.GKQ = std(data.GKQ)';
sk.GKQ = skewness(data.GKQ)';
corr_coeff = corrcoef(data.GKQ);
corrY.GKQ = corr_coeff(:,1);
av.GKQ = mean(data.GKQ)';

% OBC
% clearvars -except sd sk corrY variables;
% load('data/macroTimeSeries.mat')
% 
% timePeriod = '03/31/1986::12/31/2015';
% variables = {'Y','I','C','D','E','Q','spread'};
% 
% % de-mean spread
% raw.spread = fints(raw.spread.dates , fts2mat(raw.spread)-mean(fts2mat(raw.spread)) ,'spread' );
% compm.D =  fints(compm.D.dates , 100*(fts2mat(compm.D)-mean(fts2mat(compm.D))) ,'D' );
% compm.E =  fints(compm.E.dates , 100*(fts2mat(compm.E)-mean(fts2mat(compm.E))) ,'E' );
% data = merge( hp.y , hp.i , hp.c , compm.D, compm.E ,  hp.q , raw.spread );
% data = [ fts2mat(data.y(timePeriod)) , fts2mat(data.i(timePeriod)) ,...
%     fts2mat(data.c(timePeriod)), fts2mat(data.D(timePeriod)) ,...
%     fts2mat(data.E(timePeriod)) , fts2mat(data.q(timePeriod)) , fts2mat(data.spread(timePeriod))];

data = [ none' , none' ,none' , none' ,none' ,none' ,none' ];


% Correlations
corr_coeff = corrcoef(data);
corrY.data = corr_coeff(:,1);

% Relative Standard deviations
std_dev = std(data);
rel_std_dev = (std(data)./std_dev(:,1)); rel_std_dev(4)=NaN;
sd.data = std_dev';
av.data = mean(data);

% Moments
moments = [ mean(data)' std(data)' skewness(data)' kurtosis(data)' ];
sk.data = moments(:,3);

clearvars -except sd av sk corrY variables;

load('results_sim/newobc_order3_QMC_default.mat')
%load('results_sim/obc_phi0.mat')
spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
D_element = strmatch('D_rate',M_.endo_names,'exact');
E_element = strmatch('E_rate',M_.endo_names,'exact');
obc.Y =  100*oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
obc.I =  100*oo_.endo_simul(I_element,:);
obc.C =  100*oo_.endo_simul(C_element,:);
obc.Q =  100*oo_.endo_simul(Q_element,:);
obc.spread =  100*oo_.endo_simul(spread_element,:);
obc.D =  100*oo_.endo_simul(D_element,:);
obc.E =  100*oo_.endo_simul(E_element,:);

data.obc = [obc.Y' obc.I' obc.C' obc.D' obc.E' obc.Q' obc.spread'];
sd.obc = std(data.obc)';
sk.obc = skewness(data.obc)';
corr_coeff = corrcoef(data.obc);
corrY.obc = corr_coeff(:,1);
av.obc = mean(data.obc);



%
% clearvars -except sd sk corrY variables;
% column_titles = {'corr_YOBC','corr_YRBC','corr_YGKQ','SD_YOBC','SD_YRBC','SD_YGKQ',...
%     'sk_YOBC','sk_YRBC','sk_YGKQ'};
% table_1 = table(corrY.obc,corrY.rbc,corrY.GKQ,sd.obc,sd.rbc,sd.GKQ,sk.obc,sk.rbc,sk.GKQ,...
%     'VariableNames',column_titles,'RowNames',variables);

clearvars -except sd av sk corrY variables;
column_titles = {'corr_Ydata','corr_YOBC','corr_YRBC','corr_YGKQ','SD_data','SD_OBC','SD_RBC','SD_GKQ',...
    'sk_data','sk_OBC','sk_RBC','sk_GKQ'};
table_1 = table(corrY.data,corrY.obc,corrY.rbc,corrY.GKQ,sd.data,sd.obc,...
    sd.rbc,sd.GKQ,sk.data,sk.obc,sk.rbc,sk.GKQ,...
    'VariableNames',column_titles,'RowNames',variables);

disp('Moments and Cross Correlations');
disp(table_1);
av


%clear;
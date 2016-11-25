
load('results_sim/latest_rbc_results_1.mat')

spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
rbc.Y = oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
rbc.I = oo_.endo_simul(I_element,:)-(mean(oo_.endo_simul(I_element,:)));
rbc.C = oo_.endo_simul(C_element,:)-(mean(oo_.endo_simul(C_element,:)));
rbc.Q = oo_.endo_simul(Q_element,:)-(mean(oo_.endo_simul(Q_element,:)));
rbc.spread = oo_.endo_simul(spread_element,:);

none = zeros(1,length(rbc.Y));

variables = {'Y','I','C','D','E','Q','spread'};
data.rbc = [rbc.Y' rbc.I' rbc.C' none' none' rbc.Q' rbc.spread'];
moments.rbc = [ mean(data.rbc)' std(data.rbc)' skewness(data.rbc)' kurtosis(data.rbc)' ];
corr_coeff.rbc = corrcoef(data.rbc);

rbc_table_1 = table(moments.rbc(:,1),moments.rbc(:,2),moments.rbc(:,3),...
    moments.rbc(:,4),...
    'VariableNames',{'Mean','StdDev','Skewness','Kurtosis'},'RowNames',variables);

rbc_table_2 = table(corr_coeff.rbc(:,1),corr_coeff.rbc(:,2),corr_coeff.rbc(:,3),...
    corr_coeff.rbc(:,4),corr_coeff.rbc(:,5),corr_coeff.rbc(:,6),corr_coeff.rbc(:,7),...
    'VariableNames',variables,'RowNames',variables);

disp('**-------------------------------------------------------------------------------**');
disp('RBC Moments');
disp(rbc_table_1);
disp('RBC cross correlations');
disp(rbc_table_2);
disp('**-------------------------------------------------------------------------------**');

save('moments_data.mat','rbc');

%% GKQ
clear;

load('results_sim/latest_gkq_results_1.mat')

spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
D_element = strmatch('D_rate',M_.endo_names,'exact');
%E_element = strmatch('E_rate',M_.endo_names,'exact');
GKQ.Y = oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
GKQ.I = oo_.endo_simul(I_element,:)-(mean(oo_.endo_simul(I_element,:)));
GKQ.C = oo_.endo_simul(C_element,:)-(mean(oo_.endo_simul(C_element,:)));
GKQ.Q = oo_.endo_simul(Q_element,:)-(mean(oo_.endo_simul(Q_element,:)));
GKQ.D = oo_.endo_simul(D_element,:);
%GKQ.E = oo_.endo_simul(E_element,:);
GKQ.spread = oo_.endo_simul(spread_element,:);


none = zeros(1,length(GKQ.Y));


variables = {'Y','I','C','D','E','Q','spread'};
data.GKQ = [GKQ.Y' GKQ.I' GKQ.C' GKQ.D' none' GKQ.Q' GKQ.spread'];
moments.GKQ = [ mean(data.GKQ)' std(data.GKQ)' skewness(data.GKQ)' kurtosis(data.GKQ)' ];
corr_coeff.GKQ = corrcoef(data.GKQ);

GKQ_table_1 = table(moments.GKQ(:,1),moments.GKQ(:,2),moments.GKQ(:,3),...
    moments.GKQ(:,4),...
    'VariableNames',{'Mean','StdDev','Skewness','Kurtosis'},'RowNames',variables);

GKQ_table_2 = table(corr_coeff.GKQ(:,1),corr_coeff.GKQ(:,2),corr_coeff.GKQ(:,3),...
    corr_coeff.GKQ(:,4),corr_coeff.GKQ(:,5),corr_coeff.GKQ(:,6),corr_coeff.GKQ(:,7),...
    'VariableNames',variables,'RowNames',variables);
disp('**-------------------------------------------------------------------------------**');
disp('GK Moments');
disp(GKQ_table_1);
disp('GK cross correlations');
disp(GKQ_table_2);
disp('**-------------------------------------------------------------------------------**');

save('moments_data.mat','GKQ','-append');

%% OBC
clear;

load('results_sim/latest_newobc_results_1.mat')

spread_element = strmatch('spread',M_.endo_names,'exact');
Y_element = strmatch('Y',M_.endo_names,'exact');
I_element = strmatch('inv',M_.endo_names,'exact');
C_element = strmatch('c',M_.endo_names,'exact');
Q_element = strmatch('q',M_.endo_names,'exact');
D_element = strmatch('D_rate',M_.endo_names,'exact');
E_element = strmatch('E_rate',M_.endo_names,'exact');
obc.Y = oo_.endo_simul(Y_element,:)./(mean(oo_.endo_simul(Y_element,:)));
obc.I = oo_.endo_simul(I_element,:)-(mean(oo_.endo_simul(I_element,:)));
obc.C = oo_.endo_simul(C_element,:)-(mean(oo_.endo_simul(C_element,:)));
obc.Q = oo_.endo_simul(Q_element,:)-(mean(oo_.endo_simul(Q_element,:)));
obc.D = oo_.endo_simul(D_element,:);
obc.E = oo_.endo_simul(E_element,:);
obc.spread = oo_.endo_simul(spread_element,:);

variables = {'Y','I','C','D','E','Q','spread'};
data.obc = [obc.Y' obc.I' obc.C' obc.D' obc.E' obc.Q'  obc.spread'];
moments.obc = [ mean(data.obc)' std(data.obc)' skewness(data.obc)' kurtosis(data.obc)' ];
corr_coeff.obc = corrcoef(data.obc);

obc_table_1 = table(moments.obc(:,1),moments.obc(:,2),moments.obc(:,3),...
    moments.obc(:,4),...
    'VariableNames',{'Mean','StdDev','Skewness','Kurtosis'},'RowNames',variables);

obc_table_2 = table(corr_coeff.obc(:,1),corr_coeff.obc(:,2),corr_coeff.obc(:,3),...
    corr_coeff.obc(:,4),corr_coeff.obc(:,5),corr_coeff.obc(:,6),corr_coeff.obc(:,7),...
    'VariableNames',variables,'RowNames',variables);

disp('**-------------------------------------------------------------------------------**');
disp('OBC Moments');
disp(obc_table_1);
disp('OBC cross correlations');
disp(obc_table_2);
disp('**-------------------------------------------------------------------------------**');

save('moments_data.mat','obc','-append');

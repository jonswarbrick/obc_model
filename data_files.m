%% data files
% habits:               habitsC0, habitsC40, habitsC90, habitsH0, habitsH90
% capital adjustment:   phi0, phi1, phi4
% utility:              nonsepUtil, sepUtil, sepUtilFrisch, nonsepUtilBundles
% shocks:               shocksPsiA, shocksDeltaA , shocksConfA
% adj:                  CEE, Ireland

% 0 = rbc v newobc   , 1 = rbc v new obc x2 , 2 - rbc v new obc x3,  
% 3 = NK v NK obc x2 , 4 - NK v NK obc x3   , 5 - rbc v gkq v new obc
three_models = 5;

opt.data_files_desc = {...
    %'Baseline (non sep util)',...
    %'obc (non sep util)',...
    'Baseline (J-R util)',...
    'obc (J-R util)',...
    %'Baseline (non-sep util - bundles)',...
    %'obc (non-sep util - bundles)'...
    'gkq'
    };

if exist('mult_plots.mat')==2
load('mult_plots.mat');
data_phi = mult.data_phi;
data_shock = char(mult.data_shock);
data_habitC = mult.data_habitC;
data_habitH = mult.data_habitH;
data_stick = mult.data_stick;
data_utility = 'JRUtil';
data_adj = char(mult.data_adj);
data_phi_2 = data_phi;
data_shock_2 = data_shock;
data_habitC_2 = data_habitC;
data_habitH_2 = data_habitH;
data_stick_2 = mult.data_stick;
data_utility_2 = 'JRUtil';
data_adj_2 = data_adj;
data_phi_3 = data_phi;
data_shock_3 = data_shock;
data_habitC_3 = data_habitC;
data_habitH_3 = data_habitH;
data_stick_3 = mult.data_stick;
data_utility_3 = 'nonsepUtilBundles';
data_adj_3 = data_adj;
shockname = mult.shockname;
else

data_phi = 2;
data_shock = 'Psi';
data_habitC = 0;
data_habitH = 0;
data_utility = 'nonsepUtil';
data_adj = 'CEE';
data_phi_2 = data_phi;
data_shock_2 = data_shock;
data_habitC_2 = 90;
data_habitH_2 = 90;
data_utility_2 = 'nonsepUtil';
data_adj_2 = data_adj;
data_phi_3 = data_phi;
data_shock_3 = data_shock;
data_habitC_3 = 90;
data_habitH_3 = 90;
data_utility_3 = 'nonsepUtilBundles';
data_adj_3 = data_adj;
plot_num = 1;
shockname = 'Psi, phi 2 (CEE style)';

end
%%

if strcmp(data_shock,'Psi')
opt.epsshock = char('eps_psi');
elseif strcmp(data_shock,'Delta')
opt.epsshock = char('epsdelta');
elseif strcmp(data_shock,'Conf')
opt.epsshock = char('epsconf');
elseif strcmp(data_shock,'M')
opt.epsshock = char('epsM');
end
if data_habitC > 0
    habitC_desc = strcat(' - habitsC = ',num2str(data_habitC));
else
    habitC_desc = '';
end
if data_habitH > 0
    habitH_desc = strcat(' - habitsH = ',num2str(data_habitH));
else
    habitH_desc = '';
end
   
%%
%_irfs_order3_X3_slowIRF_NoCubature_phi
%_irfs_order3_X3_fastIRF_NoCubature_phi
%_irfs_order3_X3_fastIRF_FastCubature_phi
%_irfs_order3_X3_fastIRF_DefaultCubature_phi
%_irfs_order3_X3_fastIRF_QuasiMCCubature_phi

common_data_string = '_irfs_order3_X3_fastIRF_NoCubature_phi';

data_files_rbc = char(strcat('rbc',common_data_string,num2str(data_phi),'_',data_adj,'_shocks',data_shock,'A_habitC',num2str(data_habitC),'_habitH',num2str(data_habitH),'_',data_utility,'stick',num2str(100*data_stick)));
data_files_obc = char(strcat('obc',common_data_string,num2str(data_phi),'_',data_adj,'_shocks',data_shock,'A_habitC',num2str(data_habitC),'_habitH',num2str(data_habitH),'_',data_utility,'stick',num2str(100*data_stick)));
data_files_newobc = char(strcat('newobc',common_data_string,num2str(data_phi),'_',data_adj,'_shocks',data_shock,'A_habitC',num2str(data_habitC),'_habitH',num2str(data_habitH),'_',data_utility,'stick',num2str(100*data_stick)));
data_files_NK = char(strcat('nk',common_data_string,num2str(data_phi),'_',data_adj,'_shocks',data_shock,'A_habitC',num2str(data_habitC),'_habitH',num2str(data_habitH),'_',data_utility,'stick',num2str(100*data_stick)));
data_files_NKobc = char(strcat('nkobc',common_data_string,num2str(data_phi),'_',data_adj,'_shocks',data_shock,'A_habitC',num2str(data_habitC),'_habitH',num2str(data_habitH),'_',data_utility,'stick',num2str(100*data_stick)));
data_files_rbc_2 = char(strcat('rbc',common_data_string,num2str(data_phi_2),'_',data_adj_2,'_shocks',data_shock,'A_habitC',num2str(data_habitC_2),'_habitH',num2str(data_habitH_2),'_',data_utility_2,'stick',num2str(100*data_stick_2)));
data_files_newobc_2 = char(strcat('newobc',common_data_string,num2str(data_phi_2),'_',data_adj_2,'_shocks',data_shock,'A_habitC',num2str(data_habitC_2),'_habitH',num2str(data_habitH_2),'_',data_utility_2,'stick',num2str(100*data_stick_2)));
data_files_NK_2 = char(strcat('nk',common_data_string,num2str(data_phi_2),'_',data_adj_2,'_shocks',data_shock,'A_habitC',num2str(data_habitC_2),'_habitH',num2str(data_habitH_2),'_',data_utility_2,'stick',num2str(100*data_stick_2)));
data_files_NKobc_2 = char(strcat('nkobc',common_data_string,num2str(data_phi_2),'_',data_adj_2,'_shocks',data_shock,'A_habitC',num2str(data_habitC_2),'_habitH',num2str(data_habitH_2),'_',data_utility_2,'stick',num2str(100*data_stick_2)));
data_files_rbc_3 = char(strcat('rbc',common_data_string,num2str(data_phi_3),'_',data_adj_3,'_shocks',data_shock,'A_habitC',num2str(data_habitC_3),'_habitH',num2str(data_habitH_3),'_',data_utility_3,'stick',num2str(100*data_stick_3)));
data_files_newobc_3 = char(strcat('newobc',common_data_string,num2str(data_phi_3),'_',data_adj_3,'_shocks',data_shock,'A_habitC',num2str(data_habitC_3),'_habitH',num2str(data_habitH_3),'_',data_utility_3,'stick',num2str(100*data_stick_3)));
data_files_NK_3 = char(strcat('nk',common_data_string,num2str(data_phi_3),'_',data_adj_3,'_shocks',data_shock,'A_habitC',num2str(data_habitC_3),'_habitH',num2str(data_habitH_3),'_',data_utility_3,'stick',num2str(100*data_stick_3)));
data_files_NKobc_3 = char(strcat('nkobc',common_data_string,num2str(data_phi_3),'_',data_adj_3,'_shocks',data_shock,'A_habitC',num2str(data_habitC_3),'_habitH',num2str(data_habitH_3),'_',data_utility_3,'stick',num2str(100*data_stick_3)));



if three_models == 0
opt.num_datasets = 2;
opt.shockname = shockname;%char(strcat(strrep(strrep(opt.epsshock,'eps_psi','KQ: '),'epsdelta','Delta: '),' phi = ',num2str(data_phi),data_adj,habitC_desc,habitH_desc,' / ',data_utility));
opt.data_files = {data_files_rbc,data_files_newobc};
elseif three_models == 1
opt.num_datasets = 4;
opt.shockname = shockname;%char(strcat(strrep(strrep(opt.epsshock,'eps_psi','KQ: '),'epsdelta','Delta: '),data_adj_2,habitC_desc,habitH_desc,' / ',data_utility));
opt.data_files = {data_files_rbc,data_files_newobc,data_files_rbc_2,data_files_newobc_2};
elseif three_models == 2
opt.num_datasets = 6;
opt.shockname = shockname;%char(strcat(strrep(strrep(opt.epsshock,'eps_psi','KQ: '),'epsdelta','Delta: '),data_adj_2,habitC_desc,habitH_desc,' / ',data_utility));
opt.data_files = {data_files_rbc,data_files_newobc,data_files_rbc_2,data_files_newobc_2,data_files_rbc_3,data_files_newobc_3};
elseif three_models == 3
opt.num_datasets = 4;
opt.shockname = shockname;%char(strcat(strrep(strrep(opt.epsshock,'eps_psi','KQ: '),'epsdelta','Delta: '),data_adj_2,habitC_desc,habitH_desc,' / ',data_utility));
opt.data_files = {data_files_NK,data_files_NKobc,data_files_NK_2,data_files_NKobc_2};
elseif three_models == 4
opt.num_datasets = 6;
opt.shockname = shockname;%char(strcat(strrep(strrep(opt.epsshock,'eps_psi','KQ: '),'epsdelta','Delta: '),data_adj_2,habitC_desc,habitH_desc,' / ',data_utility));
opt.data_files = {data_files_NK,data_files_NKobc,data_files_NK_2,data_files_NKobc_2,data_files_NK_3,data_files_NKobc_3};
end
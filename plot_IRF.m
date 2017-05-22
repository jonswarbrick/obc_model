clear;
%close all;
%% Prepare Data For Plots -- User Settings
opt.variables = char('Y','inv','c','k','r','D_rate','E_rate','q','spread','kappa','delta','H','pi'); % rel = % deviation; irfs = level deviation, ; irfsAroundZero = level around 0
opt.var_paper_plot = char('rel','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfsAroundZero','irfs','irfs','irfsAroundZero','irfsAroundZero','irfs','irfs','rel','irfs');
opt.names = char('Output','Investment','Consumption','Capital','Deposit Rate','Dividend Payment','Equity Issuance','Cost of Capital','Spread','kappa','delta','Hours','Inflation');
opt.no_rows_sub_plots = 4;
opt.no_cols_sub_plots = 4;

%% Choice of data
% For large comparison (many configs)
%data_files;

% For small comparison (direct entry)
opt.epsshock = char('eps_psi');
opt.shockname = char('KQ');
% opt.epsshock = char('epsdelta');
% opt.shockname = char('Depreciation Shock');
opt.data_files_desc = {'rbc','gk','obc'};
opt.num_datasets = 3;
opt.data_files = {...
    'latest_rbc_results_1'...
    ,'latest_gk_results_1'...
    ,'latest_newobc_results_1'...
    ...,'newobc_irfs_order3_X3_fast_phi2_shocksNegPsiA_habitC70_habitH0_sepUtil'...
    ...'rbc_irfs_order3_X3_fast_phi2_shocksDeltaA_habitC70_habitH0_sepUtil'...
    ...,'newobc_irfs_order3_X3_fast_phi2_shocksDeltaA_habitC70_habitH0_sepUtil'...
    ...'rbc_irfs_order3_X5_fast_phi2_shocksDeltaA_habitC70_habitH0_sepUtil'...
    ...'rbc_irf_epsA'...
    ...,'gkq_irf_epsA'...
    ...,'newobc_irf_epsA'...
%     ,''...
%     ,''...
    };


%% Data Prep
opt.size_var = size(opt.variables);
opt.num_var = opt.size_var(1);
opt.size_epsshocks = size(opt.epsshock);
opt.number_shocks = opt.size_epsshocks(1);

index = 0;
for ii=1:opt.num_datasets
    clearvars -except ii opt model index
    index = index+1;
    eval(['load ',cell2mat(strcat('results_irf/',(opt.data_files(:,index))))]);
    for jj = 1:opt.number_shocks
    curr_shock = strtrim(opt.epsshock(jj,:));
    for kk=1:opt.num_var
        try
          eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = dynareOBC_.IRFOffsets.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
          eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = oo_.irfs.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
        catch
            try
              eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = dynareOBC_.IRFOffsets.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
              eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = oo_.irfs.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
            catch
              eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = zeros(1,60);'));
              eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = zeros(1,60);'));
            end
        end
    end
    end
    jj = 1;
    for jj = 1:opt.num_var
        eval( strcat('irfs.',strtrim(opt.variables(jj,:)),' = offset.',strtrim(opt.variables(jj,:)),'+irfsAroundZero.',strtrim(opt.variables(jj,:)),';'));
         if strcmp(strtrim(opt.var_paper_plot(jj,:)),'rel')
         eval( strcat('rel.',strtrim(opt.variables(jj,:)),' = irfs.',strtrim(opt.variables(jj,:)),'./offset.',strtrim(opt.variables(jj,:)),';'));
         end
        curr_var = strtrim(opt.variables(jj,:));
        plot_type = strtrim(opt.var_paper_plot(jj,:));
        model.abs(:,:,index,jj) = eval(['irfs.',curr_var,'(:,1:60)']);
        model.aroundZero(:,:,index,jj) = eval(['irfsAroundZero.',curr_var,'(:,1:60)']);
        model.paper(:,:,index,jj) = eval([plot_type,'.',curr_var,'(:,1:60)']);
    end
end
save('data.mat','model')
clearvars -except opt

%% Plots
load('data.mat')

for shocks = 1:opt.number_shocks  
    h = figure;
    for vars = 1:opt.num_var
        subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,vars),
        plot(model.paper(shocks,1:60,1,vars),'Color','r','LineStyle','-','LineWidth',1); hold on;
        if opt.num_datasets >1
        plot(model.paper(shocks,1:60,2,vars), 'b','LineStyle','-', 'LineWidth',1);
        end
        if opt.num_datasets >2
        plot(model.paper(shocks,1:60,3,vars),'Color','k','LineStyle','--','LineWidth',1); hold on;
        end
        if opt.num_datasets >3
        plot(model.paper(shocks,1:60,4,vars), 'b','LineStyle','--', 'LineWidth',1);
        end
        if opt.num_datasets >4
        plot(model.paper(shocks,1:60,5,vars),'Color','r','LineStyle',':','LineWidth',1); hold on;
        end
        if opt.num_datasets >5
        plot(model.paper(shocks,1:60,6,vars), 'b','LineStyle',':', 'LineWidth',1);
        end
        %xlabel('Quarters');
        %ylabel('% dev from SS');
        grid off
        titlename=strtrim(opt.names(vars,:));
        title(titlename,'FontSize',10)
     end
     axis tight;
     legend(opt.data_files_desc,'Location', 'best', 'Orientation', 'horizontal')
     %legend(series1.name,series2.name,series3.name,'Location', 'best', 'Orientation', 'horizontal')
     %legend(series1.name,series2.name,'Location', 'best', 'Orientation', 'horizontal')
     set(findall(gcf,'-property','FontSize'),'FontSize',12)
     [ax4,h3]=suplabel(strtrim(opt.shockname(shocks,:))  ,'t');
     set(h3,'FontSize',14)
end

set(h, 'WindowStyle', 'docked')

delete('data.mat')

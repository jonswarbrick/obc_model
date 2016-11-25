clear; close all;
%% Prepare Data For Plots -- User Settings

opt.variables = char('Y','I','C','B','D','E','spread','delta','a','G'); 
opt.names = char('Output','Investment','Consumption','Deposits','Dividend Payment','Equity Issuance','Spread','delta','A','G');
opt.sim_start = 450;
opt.sim_end = 550;

opt.num_datasets = 2;
% type: 1 = rbc, 2 = GK, 3 = GKQ , 4 = obc
opt.data_file(1,:) = cellstr('rbc_last_sim_results_1.mat');
opt.data_file(2,:) = cellstr('obc_last_sim_results_1.mat');
opt.data_file(3,:) = cellstr('GKQ_last_sim_results_1.mat');
opt.data_file(4,:) = cellstr('rbc_last_sim_results_1.mat');
opt.data_file_desc(1,:) = cellstr('rbc');
opt.data_file_desc(2,:) = cellstr('obc');
opt.data_file_desc(3,:) = cellstr('GK - Equity');
opt.data_file_desc(4,:) = cellstr('GK');
opt.data_file_type(1) = 1;
opt.data_file_type(2) = 4;
opt.data_file_type(3) = 3;
opt.data_file_type(4) = 2;

% Other numbers
opt.no_rows_sub_plots = 4;
opt.no_cols_sub_plots = 3;

%% Data Prep
opt.size_var = size(opt.variables);
opt.num_var = opt.size_var(1);

for ii=1:opt.num_datasets
    clearvars -except ii opt model
    eval(['load ',cell2mat(opt.data_file(ii,:))]);
    
    if opt.data_file_type(ii) == 1
        plot_SIM_rbc;
        model.abs_rbc = abs_rbc_Y;
    elseif opt.data_file_type(ii) == 2
        plot_SIM_GK;
    elseif opt.data_file_type(ii) == 3 
        plot_SIM_GKQ;
    elseif opt.data_file_type(ii) == 4 
        plot_SIM_obc;
        model.abs_obc = abs_obc_Y;
    end
    
    for jj = 1:opt.num_var
        curr_var = strtrim(opt.variables(jj,:));
        model.sim(:,jj,ii) = eval(['sim.',char(curr_var)]);
    end
end

save('data_sim.mat','model')
clearvars -except opt

%% Plots
load('data_sim.mat')
 
    h = figure;
    figure(h);
    set(h, 'Position', [50 , 50, 1000, 800]);
    for vars = 1:opt.num_var
        subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,vars),
        plot(model.sim(opt.sim_start:opt.sim_end,vars,1),'Color','k','LineStyle','-','LineWidth',2); hold on;
        if opt.num_datasets > 1
        plot(model.sim(opt.sim_start:opt.sim_end,vars,2), 'm','LineStyle','--', 'LineWidth',2);
        end
        if opt.num_datasets > 2
        plot(model.sim(opt.sim_start:opt.sim_end,vars,3),'Color','b','LineStyle','--','LineWidth',2); hold on;
        end
        if opt.num_datasets > 3
        plot(model.sim(opt.sim_start:opt.sim_end,vars,4), '-m','LineStyle','--', 'LineWidth',2);
        end
        %xlabel('Quarters');
        %ylabel('% dev from SS');
        grid off
        titlename=strtrim(opt.variables(vars,:));
        title(titlename,'FontSize',10)
    end
    legend('RBC','GK','GK - equity','OBC','Location', 'best', 'Orientation', 'horizontal')
    subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,vars+1) 
    difference = (model.abs_rbc - model.abs_obc)/mean(model.abs_rbc);
    plot(difference(opt.sim_start:opt.sim_end),'Color','r','LineStyle','-','LineWidth',2);
    grid off
    title('Y(RBC) - Y(OBC)','FontSize',10)
    axis tight;
    % legend(strtrim(opt.data_file_desc(1,:)),strtrim(opt.data_file_desc(2,:)),strtrim(opt.data_file_desc(3,:)),strtrim(opt.data_file_desc(4,:)),'Location', 'best', 'Orientation', 'horizontal')
     %legend(series1.name,series2.name,series3.name,'Location', 'best', 'Orientation', 'horizontal')
     %legend(series1.name,series2.name,'Location', 'best', 'Orientation', 'horizontal')
     set(findall(gcf,'-property','FontSize'),'FontSize',12)
     %[ax4,h3]=suplabel(strtrim(opt.shockname(shocks,:))  ,'t');
     %set(h3,'FontSize',14)

%      difference = (model.abs_rbc - model.abs_obc)/mean(model.abs_rbc);
%         h = figure;
%         figure(h);
%         set(h, 'Position', [50 , 50, 1000, 800]);
%         plot(difference(opt.sim_start:opt.sim_end),'Color','r','LineStyle','-','LineWidth',2);
%         grid off
%         title('Y (rbc) - Y (obc) (percent of RBC mean Y)','FontSize',10)
%      axis tight;
%      set(findall(gcf,'-property','FontSize'),'FontSize',12)
%      
clear;
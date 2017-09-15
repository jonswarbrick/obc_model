function [ ] = plotter( opt )
%IRF_PLOTTER plots IRFs for appendix

% Data Prep
opt.size_var = size(opt.variables);
opt.num_var = opt.size_var(1);
opt.size_epsshocks = size(opt.epsshock);
opt.number_shocks = opt.size_epsshocks(1);
opt.num_periods = 40;

index = 0;
for ii=1:opt.num_datasets
    clearvars -except ii opt model index
    index = index+1;
    eval(['load ',cell2mat(strcat('results_irf/',(opt.data_files(:,index))))]);
    for jj = 1:opt.number_shocks
        curr_shock = strtrim(opt.epsshock(jj,:));
        % Y
        eval( strcat('offset.Y(jj,:) = dynareOBC_.IRFOffsets.Y_',curr_shock,';') );
        if max(abs(eval( strcat('oo_.irfs.Y_',curr_shock) )))>1e-9
            eval( strcat('irfsAroundZero.Y(jj,:) = oo_.irfs.Y_',curr_shock,';') );
        else
            irfsAroundZero.Y(jj,:) = zeros(1,60);
        end
        % I
        eval( strcat('offset.inv(jj,:) = dynareOBC_.IRFOffsets.inv_',curr_shock,';') );
        if max(abs(eval( strcat('oo_.irfs.inv_',curr_shock) )))>1e-9
            eval( strcat('irfsAroundZero.inv(jj,:) = oo_.irfs.inv_',curr_shock,';') );
        else
            irfsAroundZero.inv(jj,:) = zeros(1,60);
        end
        % H
        eval( strcat('offset.H(jj,:) = dynareOBC_.IRFOffsets.H_',curr_shock,';') );
        if max(abs(eval( strcat('oo_.irfs.H_',curr_shock) )))>1e-9
            eval( strcat('irfsAroundZero.H(jj,:) = oo_.irfs.H_',curr_shock,';') );
        else
            irfsAroundZero.H(jj,:) = zeros(1,60);
        end
        % spread
        eval( strcat('offset.spread(jj,:) = dynareOBC_.IRFOffsets.spread_',curr_shock,';') );
        if max(abs(eval( strcat('oo_.irfs.spread_',curr_shock) )))>1e-9
            eval( strcat('irfsAroundZero.spread(jj,:) = oo_.irfs.spread_',curr_shock,';') );
        else
            irfsAroundZero.spread(jj,:) = zeros(1,60);
        end
        % D
        if strcmp(dynareOBC_.BaseFileName,'rbc_psi') || strcmp(dynareOBC_.BaseFileName,'rbc_a')
            offset.D_rate(jj,:) = zeros(1,60);
            irfsAroundZero.D_rate(jj,:) = zeros(1,60);
        elseif strcmp(dynareOBC_.BaseFileName,'obc_psi') || strcmp(dynareOBC_.BaseFileName,'obc_a')
            eval( strcat('temp_Q = exp(dynareOBC_.IRFOffsets.q_',curr_shock,') + exp(oo_.irfs.q_',curr_shock,');'));
            eval( strcat('temp_K = exp(dynareOBC_.IRFOffsets.k_',curr_shock,') + exp(oo_.irfs.k_',curr_shock,');'));
            eval( strcat('temp_B = exp(dynareOBC_.IRFOffsets.b_',curr_shock,') + exp(oo_.irfs.b_',curr_shock,');'));
            eval( strcat('temp_D = (dynareOBC_.IRFOffsets.D_',curr_shock,') + (oo_.irfs.D_',curr_shock,');'));
            eval( strcat('offset.D_rate(jj,:) = dynareOBC_.IRFOffsets.D_',curr_shock,'./( (exp(dynareOBC_.IRFOffsets.q_',curr_shock,')).*(exp(dynareOBC_.IRFOffsets.k_',curr_shock,'))- exp(dynareOBC_.IRFOffsets.b_',curr_shock,') );'));
            irfsAroundZero.D_rate(jj,:) = max(0,temp_D ./ (temp_Q.*temp_K - temp_B)) - offset.D_rate(jj,:);
        else
            eval( strcat('offset.D_rate(jj,:) = dynareOBC_.IRFOffsets.D_rate_',curr_shock,';') );
            if max(abs(eval( strcat('oo_.irfs.D_rate_',curr_shock) )))>1e-9
                eval( strcat('irfsAroundZero.D_rate(jj,:) = oo_.irfs.D_rate_',curr_shock,';') );
            else
                irfsAroundZero.D_rate(jj,:) = zeros(1,60);
            end
        end
        % E
        if strcmp(dynareOBC_.BaseFileName,'rbc_psi') || strcmp(dynareOBC_.BaseFileName,'rbc_a')
            offset.E_rate(jj,:) = zeros(1,60);
            irfsAroundZero.E_rate(jj,:) = zeros(1,60);
        else
            eval( strcat('offset.E_rate(jj,:) = dynareOBC_.IRFOffsets.E_rate_',curr_shock,';') );
            if max(abs(eval( strcat('oo_.irfs.E_rate_',curr_shock) )))>1e-9
                eval( strcat('irfsAroundZero.E_rate(jj,:) = oo_.irfs.E_rate_',curr_shock,';') );
            else
                irfsAroundZero.E_rate(jj,:) = zeros(1,60);
            end
        end
    end
    jj = 1;
    for jj = 1:opt.num_var
        eval( strcat('irfs.',strtrim(opt.variables(jj,:)),' = offset.',strtrim(opt.variables(jj,:)),'+irfsAroundZero.',strtrim(opt.variables(jj,:)),';'));
         if strcmp(strtrim(opt.var_paper_plot(jj,:)),'rel') || strcmp(strtrim(opt.var_paper_plot(jj,:)),'relAroundZero')
         eval( strcat('rel.',strtrim(opt.variables(jj,:)),' = irfs.',strtrim(opt.variables(jj,:)),'./offset.',strtrim(opt.variables(jj,:)),';'));
         eval( strcat('relAroundZero.',strtrim(opt.variables(jj,:)),' = rel.',strtrim(opt.variables(jj,:)),'-1;'));
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

h = figure;

subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,1),
    plot(100*model.paper(1,1:opt.num_periods,1,1),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1); hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,1),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,1),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1); 
    gap = [' '; ' '];
    xlabel(gap);
    set(gca,'TickLabelInterpreter','latex')
    grid off
    titlename=strtrim(opt.names(1,:));
    title(titlename,'Interpreter','latex')
subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,2),
    plot(100*model.paper(1,1:opt.num_periods,1,2),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1); hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,2),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,2),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1); 
    gap = [' '; ' '];
    xlabel(gap);
    set(gca,'TickLabelInterpreter','latex')
    grid off
    titlename=strtrim(opt.names(2,:));
    title(titlename,'Interpreter','latex')
subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,3),
    plot(100*model.paper(1,1:opt.num_periods,1,3),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1); hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,3),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,3),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1); 
    gap = [' '; ' '];
    xlabel(gap);
    if opt.plot_num==1
    ylim([-2.1 -1]);
    elseif opt.plot_num==2
    ylim([1 2.3]);
    elseif opt.plot_num==3
    ylim([0 .65]);
    elseif opt.plot_num==4
    ylim([-.65 0]);
    end
    set(gca,'TickLabelInterpreter','latex')
    grid off
    titlename=strtrim(opt.names(3,:));
    title(titlename,'Interpreter','latex')
subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,4),
    plot(100*model.paper(1,1:opt.num_periods,1,4),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1); hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,4),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,4),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1);
    gap = [' '; ' '];
    xlabel(gap);
    if opt.plot_num==1
    ylim([-.05 .5]);
    elseif opt.plot_num==2
    ylim([-.4 .04]);
    elseif opt.plot_num==3
    ylim([-.1 .02]);
    elseif opt.plot_num==4
    ylim([-.02   .1]);
    end
    set(gca,'TickLabelInterpreter','latex')
    grid off
    titlename=strtrim(opt.names(4,:));
    title(titlename,'Interpreter','latex')
subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,5),
    yyaxis left
    plot(100*model.paper(1,1:opt.num_periods,1,5),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1);  hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,5),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    if opt.plot_num==1
    ylim([-.005 .1]);
    elseif opt.plot_num==2
    ylim([-.05 1]);
    elseif opt.plot_num==3
    ylim([-.005 .2]);
    elseif opt.plot_num==4
    ylim([-.005 .2]);
    end
    set(gca,'YColor',[0 0 0]);
    yyaxis right
    plot(100*model.paper(1,1:opt.num_periods,3,5),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1); 
    if opt.plot_num==1
    ylim([-.05 2]);
    elseif opt.plot_num==2  
    ylim([-.05 2]);
    elseif opt.plot_num==3
    ylim([-.05 2]);
    elseif opt.plot_num==4
    ylim([-.05 2]);
    end
    gap = [' '; ' '];
    xlabel(gap);
    set(gca,'TickLabelInterpreter','latex')
    grid off
    titlename=strtrim(opt.names(5,:));
    title(titlename,'Interpreter','latex')
subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,6),
    if opt.plot_num==1
    yyaxis left
    plot(100*model.paper(1,1:opt.num_periods,1,6),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1);  hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,6),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    ylim([-.1 .1]);
    set(gca,'YColor',[0 0 0]);
    yyaxis right
    plot(100*model.paper(1,1:opt.num_periods,3,6),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1);    
    elseif opt.plot_num==2
    plot(100*model.paper(1,1:opt.num_periods,1,6),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1);  hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,6),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,6),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1);    
    ylim([-.1 1.5]);
    elseif opt.plot_num==3
    plot(100*model.paper(1,1:opt.num_periods,1,6),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1);  hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,6),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    plot(100*model.paper(1,1:opt.num_periods,3,6),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1);    
    ylim([-.1 .25]);
    elseif opt.plot_num==4
    yyaxis left
    plot(100*model.paper(1,1:opt.num_periods,1,6),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1);  hold on;
    plot(100*model.paper(1,1:opt.num_periods,2,6),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
    ylim([-.01 .01]);
    set(gca,'YColor',[0 0 0]);
    yyaxis right
    plot(100*model.paper(1,1:opt.num_periods,3,6),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1);    
    end
    gap = [' '; ' '];
    xlabel(gap);
    set(gca,'TickLabelInterpreter','latex');
    grid off
    titlename=strtrim(opt.names(6,:));
    title(titlename,'Interpreter','latex')
 %axis tight;
 set (h,'Position',opt.plot_size)
 if ~opt.nolegend
 legend(opt.data_files_desc,'position', [0.275 0.05 0.45 0.05], 'Orientation', 'horizontal','Interpreter','latex')
 end
 img_filename = strcat('plot_files\IRF_',opt.shock_shortname(1,:),'_',opt.type_desc);
 if exist('plot_files')~=7
    mkdir('plot_files')
 end
 print(h,img_filename,'-r300','-dpng')
 print(h,img_filename,'-r300','-depsc')

delete('data.mat')

end


function [ ] = IRF_plotter_paper_nogap( opt )
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
    for kk=1:opt.num_var
        %try
          %eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = dynareOBC_.IRFOffsets.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
         % eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = oo_.irfs.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
        %catch
        try
            eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = dynareOBC_.IRFOffsets.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
            if max(abs(eval( ( strcat('oo_.irfs.',strtrim(opt.variables(kk,:)),'_',curr_shock)))))>1e-9
              eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = oo_.irfs.',strtrim(opt.variables(kk,:)),'_',curr_shock,';'));
            else
              eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = zeros(1,60);'));
            end
        catch
            eval( strcat('offset.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = zeros(1,60);'));
            eval( strcat('irfsAroundZero.',strtrim(opt.variables(kk,:)),'(',num2str(jj),',:) = zeros(1,60);'));
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

for shocks = 1:opt.number_shocks  
    h = figure;
    for vars = 1:opt.num_var
        subplot(opt.no_rows_sub_plots,opt.no_cols_sub_plots,vars),
        plot(model.paper(shocks,1:opt.num_periods,1,vars),'Color',[0 .6 .4],'LineStyle','-','LineWidth',1); hold on;
        if opt.num_datasets >1
        plot(model.paper(shocks,1:opt.num_periods,2,vars),'color',[.2 .4 .6],'LineStyle','--', 'LineWidth',1);
        end
        if opt.num_datasets >2
        plot(model.paper(shocks,1:opt.num_periods,3,vars),'Color',[0.6 0 0],'LineStyle',':','LineWidth',1); hold on;
        end
        if opt.num_datasets >3
        plot(model.paper(shocks,1:opt.num_periods,4,vars), 'color',[0 .6 .4],'LineStyle','--','LineWidth',1);
        end
        if opt.num_datasets >4
        plot(model.paper(shocks,1:opt.num_periods,5,vars), 'color',[.2 .4 .6],'LineStyle','--','LineWidth',1); hold on;
        end
        if opt.num_datasets >5
        plot(model.paper(shocks,1:opt.num_periods,6,vars), 'color',[0.6 0 0],'LineStyle','--','LineWidth',1);
        end
        %ylabel('% dev from SS');
        grid off
        titlename=strtrim(opt.names(vars,:));
        title(titlename,'FontSize',10,'Interpreter','latex')
     end
     axis tight;
     %[~,h3]=suplabel(strtrim(opt.shockname(shocks,:))  ,'t');
     %set(h3,'FontSize',14)
     set(findall(gcf,'-property','FontSize'),'FontSize',12,'Fontname','Palatino Linotype');
     set (h,'Position',opt.plot_size)
     if ~opt.nolegend
     legend(opt.data_files_desc,'position', [0.275 0.03 0.45 0.05], 'Orientation', 'horizontal')
     end
     img_filename = strcat('C:\Users\Jonathan\Dropbox\PhD\Thesis\OBC\OBC\images\IRF_',opt.shock_shortname(shocks,:),'_',opt.type_desc);
     %h.PaperUnits = 'points';
     %h.PaperPosition = opt.plot_size; 
     print(h,img_filename,'-r300','-dpng')
     print(h,img_filename,'-r300','-deps')
end

delete('data.mat')

end


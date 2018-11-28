clear; close all;


for mod=1:3
    if mod==1
        load('results_sim/obc.mat') 
        disp('**- Our model -------------------------------------------------------------------**');
    elseif mod==2
        load('results_sim/rbc.mat') 
        disp('**- RBC -------------------------------------------------------------------------**');
    elseif mod==3
        load('results_sim/gkq.mat') 
        disp('**- GK --------------------------------------------------------------------------**');
    end

    y = dynareOBC_.MLVSimulationWithBounds.y;
    inv = (oo_.endo_simul(strmatch('inv',M_.endo_names,'exact'),:));
    c = (oo_.endo_simul(strmatch('c',M_.endo_names,'exact'),:));
    spread = (oo_.endo_simul(strmatch('spread',M_.endo_names,'exact'),:));
    D_rate = dynareOBC_.MLVSimulationWithBounds.D_rate;
    E_rate = dynareOBC_.MLVSimulationWithBounds.E_rate;
    lev = dynareOBC_.MLVSimulationWithBounds.lev;
    [~,hp.y] = hpfilter(y,1600);
    [~,hp.inv] = hpfilter(inv,1600);
    [~,hp.c] = hpfilter(c,1600);
    [y_ac,~,~] = autocorr(hp.y,1);
    lev_corr = corrcoef(lev',hp.y);
    disp(horzcat('S.D. Y = ',num2str(std(hp.y)),'| Target = 0.010563'));
    disp(horzcat('AC(1) Y = ',num2str(y_ac(2)),'| Target = 0.86255'));
    disp(horzcat('mean spread = ',num2str(mean(spread)),'| Target = 0.0057369'));
    disp(horzcat('S.D. spread = ',num2str(std(spread)),'| Target = 0.0017812'));
    disp(horzcat('Investment skewness: ',num2str(skewness(hp.inv))));
    disp(horzcat('Spread skewness: ',num2str(skewness(spread))));
    disp(horzcat('Asset / Equity (corr with Y): ',num2str(lev_corr(1,2))));

    % Table for paper
    variables = {'Y','I','C','D','E','spread'};
    data = [ hp.y.*100 , hp.inv.*100 , hp.c.*100 , D_rate'.*100 , E_rate'.*100 , spread'.*100  ];
    corr_coeff = corrcoef(data);
    moments = [ mean(data)' std(data)' skewness(data)' kurtosis(data)' ];

    paper_data = [ ...
        corr_coeff(1,1) , corr_coeff(1,2) , corr_coeff(1,3) ,...
        corr_coeff(1,4) , corr_coeff(1,5) , corr_coeff(1,6) ; ...
        ( moments(:,2) )' ; ( moments(:,3) )' ]; 

    tab_paper = table(paper_data(:,1),paper_data(:,2),paper_data(:,3),...
        paper_data(:,4),paper_data(:,5),paper_data(:,6),...
        'VariableNames',variables,'RowNames',{'corr','sd','skew'});

    disp('**-------------------------------------------------------------------------------**');
    disp(tab_paper);
    disp('**-------------------------------------------------------------------------------**');

    if strcmp(dynareOBC_.BaseFileName,'obc_sim')
                mv = (oo_.endo_simul(strmatch('mv',M_.endo_names,'exact'),:));
                mv(mv<1e-8) = 0;
                mv(mv>1e-8) = 1;
                binding_periods = mean(mv);
                disp(horzcat('Constraint binding in ',num2str(100*binding_periods),'% of periods'));
    end

end

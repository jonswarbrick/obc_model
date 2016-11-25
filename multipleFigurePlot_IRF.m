clear;% close all; 

mult.data_shock_vec     = {'Psi'};%{'Psi','Delta'};
mult.data_phi_vec       = [ 0 , 2 ];%[ 0 , 2 , 2 ];
mult.data_adj_vec       = {'CEE'};%{'CEE','CEE','Ireland'};
mult.data_habitC_vec    = 70;%[  0 , 70 ];
mult.data_habitH_vec    = [  0  ];
mult.data_stick_vec     = 0;%[ 0.01 , 0.1 , 0.4 , 0.7];

for mult_i = 1:length(mult.data_shock_vec)
    for mult_j = 1:length(mult.data_phi_vec)
        for mult_k = 1:length(mult.data_habitC_vec)
        	for mult_l = 1:length(mult.data_habitH_vec)
                for mult_m = 1:length(mult.data_stick_vec)
                    mult.data_shock     = mult.data_shock_vec(mult_i);
                    mult.data_phi       = mult.data_phi_vec(mult_j);
                    mult.data_adj       = mult.data_adj_vec(mult_j);
                    mult.data_habitC    = mult.data_habitC_vec(mult_k);
                    mult.data_habitH    = mult.data_habitH_vec(mult_l);
                    mult.data_stick     = mult.data_stick_vec(mult_m);

                    adj_desc = {'phi 0','phi 2 (CEE style)','phi 2 (Ireland style)'};
                    habs_desc = {'no habits','habits in C' ; 'habits in H','habits in C and H' };
                    stick_desc = 'Flexible';%{'[Calvo: 0.01]','[Calvo: 0.1]' ; '[Calvo: 0.4]','[Calvo: 0.7]' };

                    mult.shockname = strcat(char(mult.data_shock),' , ',...
                        char(adj_desc(mult_j)),' , ',char(habs_desc(mult_l,mult_k)),...
                        char(stick_desc(mult_m)));
mult.shockname='Quasi Monte Carlo';
                    save('mult_plots.mat','mult','mult_i','mult_j','mult_k','mult_l','mult_m')

                    plot_IRF;

                    load('mult_plots.mat');
                end
            end
        end
    end
end

delete('mult_plots.mat')
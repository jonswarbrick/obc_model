clear; close all;
fclose('all');
if exist('mutliple_log.txt','file') == 2
delete('mutliple_log.txt')
end
% 1 = rbc, 2 = gkq, 3 = obc, 4 = nk, 5 = nkobc, 6 = newobc
mult.models_to_run = [ 1 , 2 , 6 ];
multopts.parameter_habits_C_vec = [ 0 , .7 ];
multopts.parameter_habits_H_vec = [ 0 ];%, .7 ];
multopts.parameter_Phi_vec = [0];
% 1 = non-separable, 2 = additive type 1 , 3 = additive type 2 , 4 =
% non-separable habits on bundles , 5 J-R
multopts.mult_util_type = [ 1,  5 ];
% 1 = KQ, 2 = delta, 3 = conf, 4 = MP
multopts.mult_shock_choice = [ 1 , 2 ];
% 1 = CEE, 2 = Ireland (2003)
multopts.mult_adj_type = [ 1 ];
multopts.length_mult_i = length(multopts.mult_util_type);  
multopts.length_mult_j = length(multopts.parameter_habits_C_vec);
multopts.length_mult_k = length(multopts.parameter_habits_H_vec);
multopts.length_mult_l = length(multopts.mult_shock_choice);
multopts.length_mult_m = length(multopts.parameter_Phi_vec);
multopts.length_mult_n = length(multopts.mult_adj_type);
multopts.total_loops = multopts.length_mult_i*multopts.length_mult_j*multopts.length_mult_k*multopts.length_mult_l*multopts.length_mult_m*multopts.length_mult_n;           

% Parameters
multopts.parameter_xi = 0;%[0.01 , 0.1 , 0.4 , 0.7];

multopts.length_mult_par = length(multopts.parameter_xi);

multopts.loop_num = 1;

for mult_i = 1:multopts.length_mult_i
  for mult_j = 1:multopts.length_mult_j 
    for mult_k = 1:multopts.length_mult_k
      for mult_l = 1:multopts.length_mult_l
          for mult_m = 1:multopts.length_mult_m
            for mult_n = 1:multopts.length_mult_n
                for mult_par = 1:multopts.length_mult_par
              
save('mult_loop.mat','mult_j','mult_i','mult_k','mult_l','mult_m','mult_n','mult_par','multopts');
 
mult.utility_type = multopts.mult_util_type(mult_i);   
mult.shock_choice = multopts.mult_shock_choice(mult_l);   
mult.adj_type = multopts.mult_adj_type(mult_n);     
mult.parameter_habits_C = multopts.parameter_habits_C_vec(mult_j);
mult.parameter_habits_H = multopts.parameter_habits_H_vec(mult_k);
mult.parameter_Phi = multopts.parameter_Phi_vec(mult_m);
mult.parameter_xi = multopts.parameter_xi(mult_par);
mult.total_loops = multopts.total_loops;
mult.loop_num = multopts.loop_num;

str_util = {'nonsepUtil';'sepUtil';'sepUtilFrisch';'nonsepUtilBundles';'JRutil'};
str_shocks = {'Psi';'Delta';'Conf';'M'};
str_adj = {'CEE';'Ireland'};
str_phi = num2str(mult.parameter_Phi);
str_habC = num2str(100*mult.parameter_habits_C);
str_habH = num2str(100*mult.parameter_habits_H);
str_xi = num2str(100*mult.parameter_xi);

%_irfs_order3_X3_slowIRF_NoCubature_phi
%_irfs_order3_X3_fastIRF_NoCubature_phi
%_irfs_order3_X3_fastIRF_FastCubature_phi
%_irfs_order3_X3_fastIRF_DefaultCubature_phi
%_irfs_order3_X3_fastIRF_QuasiMCCubature_phi

mult.temp_file_name = char(strcat('_irfs_order2_X3_slowIRF_nocubature_phi',str_phi,'_',str_adj(mult.adj_type),'_shocks',str_shocks(mult.shock_choice),'A_habitC',str_habC,'_habitH',str_habH,'_',str_util(mult.utility_type),'stick',str_xi));
save('mult.mat','mult');

run;

load('mult_loop.mat');
multopts.loop_num = multopts.loop_num+1;
    
                end
            end
          end
      end
    end
  end
end

%% Clean up
clear;
delete('mult.mat')
delete('mult_loop.mat')

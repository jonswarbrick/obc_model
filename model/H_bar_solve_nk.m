function F = H_bar_solve_nk(H_bar)
global M_
 
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  if strcmp(paramname , 'H_bar')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end

F = H_bar*varrho*(1-epsilonC)*( 1/Disp_bar - gySS/Disp_bar - deltaSS * MC_bar*alp/Z_bar )/(MC_bar*(1-alp))-((1-epsilonH)*(1-H_bar))^(sigma_h);
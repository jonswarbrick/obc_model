function F = call_csolve_obc(~)
global M_
 
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                    %    Get the name of parameter i. 
  if strcmp(paramname , 'H_bar')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end
C_over_Y = (1 - gySS - deltaSS * alp/( 1 / ( 1 - gam ) / betta - ( 1 -  deltaSS ) / ( 1 - gam ) + gam * (  Theta * ( 1 - deltaSS ) ) ));
H_bar = csolve(@(H_bar) H_bar_solve(C_over_Y,H_bar),0.5,[],1e-8,200);
F = H_bar;
end
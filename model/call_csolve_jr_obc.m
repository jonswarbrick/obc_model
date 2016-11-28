function F = call_csolve_jr_obc(~)
global M_
 
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                    %    Get the name of parameter i. 
  if strcmp(paramname , 'H_bar')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end
K_over_Y = (alp / ( ( 1/((1-gam)*betta*( ( 1 - gam*(1- ( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ) )*(1-(1-gam)*(1-Theta)) )/(1-gam) ))  ) -1+deltaSS ));
H_bar = csolve(@(H_bar) H_bar_jr_solve(K_over_Y,H_bar),0.5,[],1e-8,200);
F = H_bar;
end
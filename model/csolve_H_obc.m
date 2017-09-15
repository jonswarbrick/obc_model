function F = csolve_H_obc(~)
global M_
 
NumberOfParameters = M_.param_nbr;         
for i = 1:NumberOfParameters                   
  paramname = deblank(M_.param_names(i,:));       
  if strcmp(paramname , 'H_bar')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);  
end
K_over_Y = (alp / ( ( 1/((1-gam)*betta*( ( 1 - gam*(1- ( (1-gam)*betta*MD1_bar/(1+(1-gam)*betta*MD1_bar) ) )*(1-(1-gam)*(1-Theta)) )/(1-gam) ))  ) -1+deltaSS ));
H_bar = csolve(@(H_bar) H_bar_solve(K_over_Y,H_bar),0.5,[],1e-8,200);
F = H_bar;
end
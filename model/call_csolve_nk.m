function F = call_csolve_nk(~)
global M_
 
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                    %    Get the name of parameter i. 
  if strcmp(paramname , 'H_bar')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end
H_bar = csolve(@(H_bar) H_bar_solve_nk(H_bar),0.5,[],1e-8,200);
F = H_bar;
end
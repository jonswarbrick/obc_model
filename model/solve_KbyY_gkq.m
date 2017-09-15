function K_by_Y = solve_KbyY_gkq(~)
global M_

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i.                    %    Get the name of parameter i.
  if strcmp(paramname , 'H_bar')
      continue
  end
  if strcmp(paramname , 'K_by_Y')
      continue
  end
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=[1.1 , log(0.3)];
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxIter', 5000, 'MaxFunEvals', 5000, 'Display','off' );
[x,fval] =fsolve(@solve_gkq,x0,options);

K_by_Y = exp( x(1) );
H_bar = exp( x(2) );

end

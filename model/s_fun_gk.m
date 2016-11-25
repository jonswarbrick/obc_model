function y=s_fun_gk(x)
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

load('../opts.mat');

K_by_Y = exp( x(1) );
H_bar = exp( x(2) );

H = H_bar;
R = 1 / betta;
Z = alp/K_by_Y;
I_by_Y = deltaSS * K_by_Y;
C_by_Y = 1 - gySS - I_by_Y;
Y = ( Ass * H ) * K_by_Y ^ ( alp / ( 1 - alp ) );
K = K_by_Y * Y;
RK = (Z + (1 - deltaSS));
S = K;
N = ((sigmaB + xiB)*RK - R*sigmaB)*S/(1-R*sigmaB);
phi_GK = S/N;
m = phi_GK*Theta;
Omega = 1 - sigmaB + sigmaB*m;	
C = C_by_Y*Y;

if utility_type == 1
    UH = -varrho*( (C*(1 - epsilonC))^((1-varrho)*(1-sigma_c)))*(((1-H)*(1- epsilonH))^(varrho*(1-sigma_c)-1));
    UC = (1-varrho)*( (C*(1 - epsilonC))^((1-varrho)*(1-sigma_c)-1))*(((1-H)*(1- epsilonH))^(varrho*(1-sigma_c)));
elseif utility_type == 2
    UH = -varrho*((1-H)*(1- epsilonH))^(-sigma_h);
    UC = (C*(1 - epsilonC))^(-1);
elseif utility_type == 3
    UH = -(H*(1 - epsilonH))^(psi_h);
    UC = (C*(1 - epsilonC))^(-1);
elseif utility_type == 4
    UH = -varrho*((1-H)/C)^(varrho-1)*( (1-epsilonC)*C^(1-varrho)*(1-H)^varrho )^(-sigma_c);
    UC = (1-varrho)*((1-H)/C)^varrho*( (1-epsilonC)*C^(1-varrho)*(1-H)^varrho )^(-sigma_c);
elseif utility_type == 5
    habits = epsilonC*( C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)));
    UH = - (theta_jr+1-gam_jr)*varrho*C^gam_jr*H^(theta_jr-gam_jr)*(C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)) - habits )^(-sigma_c);
    UC = (1 - gam_jr*varrho*H^theta_jr*(H/C)^(1-gam_jr))*(C - varrho*H^theta_jr*(C^gam_jr*H^(1-gam_jr)) - habits )^(-sigma_c);
end

y = [m*(1-betta*Omega*(RK-R)/Theta)-betta*Omega*R;
    UH/UC+(1-alp)*Y/H];




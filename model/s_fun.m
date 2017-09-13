function y=s_fun(x)
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
C_by_Y = 1 - I_by_Y;
Y = ( Ass * H ) * K_by_Y ^ ( alp / ( 1 - alp ) );
K = K_by_Y * Y;
RK = (Z + (1 - deltaSS));
S = K;
spread = RK - R;
xE = -epsilon/kappa_GK;
Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);
N = ((sigmaB + xiB)*RK - R*sigmaB)*S/(1-R*sigmaB);
phi_GK = S/N;
Omega = 1 - sigmaB + sigmaB*Thetax*phi_GK;
mus = betta*Omega*spread;
nub = Omega;
C = C_by_Y*Y;

if utility_type == 1
    UH = -varrho*( (C*(1 - epsilonC))^((1-varrho)*(1-sigma_c)))*(((1-H)*(1- epsilonH))^(varrho*(1-sigma_c)-1));
    lambdaC = (1-varrho)*( (C*(1 - epsilonC))^((1-varrho)*(1-sigma_c)-1))*(((1-H)*(1- epsilonH))^(varrho*(1-sigma_c)));
elseif utility_type == 2
    UH = -varrho*((1-H)*(1- epsilonH))^(-sigma_h);
    lambdaC = (C*(1 - epsilonC))^(-1);
elseif utility_type == 3
    UH = -(H*(1 - epsilonH))^(psi_h);
    lambdaC = (C*(1 - epsilonC))^(-1);
elseif utility_type == 4
    UH = -varrho*((1-H)/C)^(varrho-1)*( (1-epsilonC)*C^(1-varrho)*(1-H)^varrho )^(-sigma_c);
    lambdaC = (1-varrho)*((1-H)/C)^varrho*( (1-epsilonC)*C^(1-varrho)*(1-H)^varrho )^(-sigma_c);
elseif utility_type == 5
    Xjr = C;
    UH = - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * theta_jr*varrho*Xjr*H^(theta_jr-1);
    UC = (C - varrho*H^theta_jr*Xjr)^(-sigma_c);
    UX =  - (C - varrho*H^theta_jr*Xjr)^(-sigma_c) * varrho * H^(theta_jr);
    lambdaX = UX / ( 1 - betta*( (1-gam_jr) ) );
    lambdaC = UC + lambdaX*gam_jr;
end

y = [phi_GK-nub/(Thetax-mus)
    UH/lambdaC+(1-alp)*Y/H];

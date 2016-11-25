function [ys,check] = GKQ_steadystate(ys,exe)

global M_
 
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;

x0=[log(1.5) log(0.6)];
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxIter', 100, 'MaxFunEvals', 1000, 'Display','iter' );
[x,fval] =fsolve(@s_fun,x0,options);

load('../opts.mat');

disp = 0;
pi = 0;
mc = 0;
q = 0;
Q = 1;
R = 1 / betta;
r = log( R );
a = log( Ass );
logit_delta = logit_deltaSS;

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
spread = RK - R;
xE = -epsilon/kappa_GK;
Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);
N = ((sigmaB + xiB)*RK - R*sigmaB)*S/(1-R*sigmaB);
phi_GK = S/N;
Omega = 1 - sigmaB + sigmaB*Thetax*phi_GK;	
mus = betta*Omega*spread;	
nub = Omega;
C = C_by_Y*Y;

g = log ( GSS );
k = log( K );
c = log( C_by_Y * Y );
inv = log( I_by_Y * Y );
I = I_by_Y * Y;

RK = (Z + (1 - deltaSS));
RE = R;
QE =Z/(RE  - (1 - deltaSS));
S = K;
mue = 0;
spread = RK - R;
xE = -epsilon/kappa_GK;
Thetax = Theta*(1 + epsilon*xE + kappa_GK*xE^2/2);
E = xE*S/QE;
N = ((sigmaB + xiB)*RK - R*sigmaB)*S/(1-R*sigmaB);
phi_GK = S/N;
Omega = 1 - sigmaB + sigmaB*Thetax*phi_GK;	
mus = betta*Omega*spread;	
nub = Omega;
D = (1-sigmaB)*(RK*S-R*(S-N-QE*E) - RE*QE*E); 
D_rate = D/N;
E_rate = E/N;
psi = 1;


NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
for i = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(i) ') = ' varname ';']);                %    Get the steady state vZNue of this variable.
end                                                           % End of the loop.
%%


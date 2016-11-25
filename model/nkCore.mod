# lag_Disp = exp( disp(-1) );
# y = log(Y);

1 = xi*Pi^(zzeta-1) + (1-xi)*(X1/X2)^(1-zzeta);
Disp =  lag_Disp*xi*Pi^zzeta + (1-xi)*(X1/X2)^(-zzeta);
X1 = (zzeta/(zzeta-1))*Y*MC + xi*X1(+1)*lead_Lambda*lead_Pi^zzeta;
X2 = Y + xi*X2(+1)*lead_Lambda*lead_Pi^(zzeta-1);
r = max( 0 , gamR*r(-1) + (1-gamR)*( r_bar + gamPi*(pi - pi_bar) + gamY*( y - y_bar) ) + sigmaM*epsM );


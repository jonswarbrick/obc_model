//required
var mue mus nub E Thetax QE;   // GKQ
//var m;   //GK
//for analysis
var D_rate  E_rate;

parameters C_bar K_by_Y H_bar GSS;

H_bar  = call_Hbar_gkq;
K_by_Y  = call_KbyY_gkq;

C_bar = ( 1 - gySS - deltaSS * K_by_Y )*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );
GSS = gySS*( Ass * H_bar ) * K_by_Y ^ ( alp / ( 1 - alp ) );

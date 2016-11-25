@#include "which_model.mod"
@#include "sim_type.mod"
@#include "utility_type.mod"
@#include "shock_choice.mod"
@#include "adj_type.mod"

@#include "commonVarsAndPars.mod"

@#if model_type == 1
    @#include "rbcVars.mod"
@#endif
@#if model_type == 2
    @#include "gkqVars.mod"
@#endif
@#if model_type == 3
    @#include "obcVars.mod"
@#endif
@#if model_type == 4
    @#include "commonNKVarsAndPars.mod"
    @#include "nkVars.mod"
@#endif
@#if model_type == 5
    @#include "commonNKVarsAndPars.mod"
    @#include "nkobcVars.mod"
@#endif
@#if model_type == 6
    @#include "new_obcVars.mod"
@#endif

model;
@#include "rbcCoreMlvs.mod"

@#if model_type == 1
    @#include "rbcEquations.mod"
@#endif
@#if model_type == 2
    @#include "gkqEquations.mod"
@#endif
@#if model_type == 3
    @#include "obcEquations.mod"
    @#include "obcCore.mod"
@#endif
@#if model_type == 4
    @#include "nkCore.mod"
    @#include "nkEquations.mod"
@#endif
@#if model_type == 5
    @#include "nkCore.mod"
    @#include "new_obcCore.mod"
    @#include "obc_shock_processes.mod"
@#endif
@#if model_type == 6
    @#include "obcEquations.mod"
    @#include "new_obcCore.mod"
    @#include "obc_shock_processes.mod"
@#endif

@#include "rbcCore.mod"
@#include "shock_processes.mod"

end;

@#if model_type == 1
    @#include "rbcSS.mod"
@#endif
@#if model_type == 2
   % @#include "gkqSS.mod"
@#endif
@#if model_type == 3
    @#include "obcSS.mod"
@#endif
@#if model_type == 4
    @#include "nkSS.mod"
@#endif
@#if model_type == 5
    @#include "nkobcSS.mod"
@#endif
@#if model_type == 6
    @#include "new_obcSS.mod"
@#endif

steady;
check;

shocks;
@#include "shocks_2.mod"
end;


@#if sim_type == 1
    stoch_simul( order = 2, irf = 0, periods = 1000 );
@#endif
@#if sim_type == 2
    stoch_simul( nograph , replic = 500, order = 3, irf = 60, periods = 0 , irf_shocks = ( @#include "shocks_3.mod" ) ); 
@#endif



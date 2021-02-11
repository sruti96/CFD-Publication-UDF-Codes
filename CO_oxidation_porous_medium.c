#include <udf.h>
#define CO 2

DEFINE_SR_RATE(user_rate,c,cthread,r,mole_weight,species_mf,rr)
{
    if (FLUID_THREAD_P(cthread) && THREAD_VAR(cthread).fluid.porous)
    {
        real T = C_T(c, cthread);
	    real CO_CONC = C_R(c, cthread)*species_mf[CO]/mole_weight[CO];
        *rr = r->A*exp(-(r->E)/(UNIVERSAL_GAS_CONSTANT*T))*CO_CONC;
    }
        
    else
        *rr = 0;
}
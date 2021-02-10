#include <udf.h>
#define NH3 0
#define NH3_S 1
#define O2 2
#define NO 3
#define NO2 4
#define N2 5
#define H2O 6
#define R 8.314
DEFINE_SR_RATE(LH_Rate,f,t,r,mw,yi,rr)
{
    
   	real C_NH3, theta_NH3, C_O2, C_NO, C_NO2, C_tot;
	real k1, k2, k3, k4, k5, k6, k7, k8, K;

if (FLUID_THREAD_P(t) && THREAD_VAR(t).fluid.porous)
    {
		real T = F_T(f,t);
        C_tot = C_R(f,t);
        
        C_NH3 = C_tot*yi[NH3]/mw[NH3]*10^-3; //concentrations in mol/m3
        C_O2 = C_tot*yi[O2]/mw[O2]*10^-3;
        C_NO = C_tot*yi[NO]/mw[NO]*10^-3;
        C_NO2 = C_tot*yi[NO2]/mw[NO2]*10^-3;
        theta_NH3 = yi[NH3_S]; //site coverage of NH3
        
        k1 = 4.5;
        k2 = 2.49E5*exp(-97500*(1-0.38*theta_NH3)/(R*T));
        k3 = 1.39E6*exp(-63800*(1-theta_NH3)/(R*T));
        k4 = 3.63*exp(-32100*(1-theta_NH3)/(R*T));
        k5 = 3.18E8*exp(-88000*(1-theta_NH3)/(R*T));
        k6 = 2.33E7*exp(-32100*(1-theta_NH3)/(R*T));
        k7 = 4.24E5*exp(-58300*(1-theta_NH3)/(R*T));
        k8 = 3.07E4*exp(-48200*(1-theta_NH3)/(R*T));
		KEQ = exp(5.0462 + (6343.4/T) - 2.2973*log(T) + 3.0315e-3*T - 8.2812e-7*T*T + 1.1412e-10*T*T*T);
        K = 10^-3*120*(0.001/4)//K is the conversion factor from mole/mol site. s to kmole/m2. s
        // storage capacity of NH3 = 120 mol/m3
        // Diameter of monolith channel = 0.001 m

		if (STREQ(r->name,"reaction-1"))
			*rr = k1*C_NH3*(1-theta_NH3)*K;
		if (STREQ(r->name,"reaction-2"))
			*rr = k2*theta_NH3*K;
		if (STREQ(r->name,"reaction-3"))
			*rr = k3*C_O2*C_NH3^2*K;
		if (STREQ(r->name,"reaction-4"))
			*rr = k4*(C_O2^0.5*C_NO-(C_NO2/KEQ))*K;
        if (STREQ(r->name,"reaction-5"))
            *rr = k5*C_NO*theta_NH3^2*K;
        if (STREQ(r->name,"reaction-6"))
            *rr = k6*C_NO*C_NO2*theta_NH3*K;
        if (STREQ(r->name,"reaction-7"))
            *rr = k7*C_NO2*theta_NH3*K;
        if (STREQ(r->name,"reaction-8"))
            *rr = k8*C_NO2*theta_NH3*K;
}

}

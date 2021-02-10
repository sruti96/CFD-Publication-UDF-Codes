#include <udf.h>
#define CO 2
#define O2 0
#define C3H6 8
#define NO 6
#define H2 5
#define NO2 7
#define CO2 1
#define N2 4
#define H2O 3
DEFINE_SR_RATE(LH_Rate,f0,t0,r,mw,yi,rr)
{
    Thread *t = THREAD_T0(t0);
	cell_t f = F_C0(f0,t0);
	real YCO, YO2, YC3H6, YNO, YH2, YNO2, YCO2, YN2, YH2O, Y_T;
	real XCO, XO2, XC3H6, XNO, XH2, XNO2, XCO2, XN2, XH2O;
	real G1, G2, G3, G4, G;
	real K1, K2, K3, K4, K5, KNO, KEQ;
	real k1, k2, k3;
	real P = 101300;

		real T = F_T(f,t);
		YCO = C_R(f,t)*yi[CO]/mw[CO];
		YO2 = C_R(f,t)*yi[O2]/mw[O2];
		YC3H6 = C_R(f,t)*yi[C3H6]/mw[C3H6];
		YNO = C_R(f,t)*yi[NO]/mw[NO];
		YH2 = C_R(f,t)*yi[H2]/mw[H2];
		YNO2 = C_R(f,t)*yi[NO2]/mw[NO2];
		YCO2 = C_R(f,t)*yi[CO2]/mw[CO2];
		YN2 = C_R(f,t)*yi[N2]/mw[N2];
		YH2O = C_R(f,t)*yi[H2O]/mw[H2O];
		Y_T = YCO+YO2+YC3H6+YNO+YH2+YNO2+YCO2+YN2+YH2O;
		XCO = YCO/Y_T; XO2 = YO2/Y_T; XC3H6 = YC3H6/Y_T; XNO = YNO/Y_T; 
		XH2 = YH2/Y_T; XNO2 = YNO2/Y_T; XCO2 = YCO2/Y_T; XN2 = YN2/Y_T;

		k1 = 2.03e17*exp(-6241/T)*10e-8;
		k2 = 2.8e19*exp(-10825/T)*10e-8;
		k3 = 30e8*exp(-2567/T)*10e-8;

		
		K1 = 111.91*exp(790.7/T);
		K2 = 644.32*exp(1591/T);
		K3 = 1.3e9*exp(4811.16/T);
		K4 = 10955.5*exp(539.1/T);
		K5 = 2.3e8*exp(-8083.6/T);
		KNO = 1600;
		KEQ = exp(5.0462 + (6343.4/T) - 2.2973*log(T) + 3.0315e-3*T - 8.2812e-7*T*T + 1.1412e-10*T*T*T);

		G1 = pow(1 + K1*XCO + K2*XC3H6, 2);
		G2 = 1 + K3*pow(XCO*XC3H6,2);
		G3 = 1 + K4*XNO;
		G4 = 1 + K5*XO2;
		G=T*G1*G2*G3;

		if (STREQ(r->name,"reaction-1"))
			*rr = k1*XCO*XO2/G;
		if (STREQ(r->name,"reaction-2"))
			*rr = k2*XC3H6*XO2*XNO*KNO/(G*G4);
		if (STREQ(r->name,"reaction-3"))
			*rr = k2*XC3H6*XO2/G;
		if (STREQ(r->name,"reaction-4"))
			*rr = k3*((XNO*pow(XO2,0.5))-(300*XNO2/(KEQ*pow(P,0.5))))/G;
		if (STREQ(r->name,"reaction-5"))
			*rr = k1*XH2*XO2/G;       
    
}
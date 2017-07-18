#include "profiles.h"
#include <complex.h>

#define RT2LN2 1.1774100225154747     // sqrt(2ln2)
#define RTPI 1.7724538509055159       // sqrt(pi)
#define RT2PI 2.5066282746310002      // sqrt(2pi)
#define RT2 1.4142135623730951        // sqrt(pi)

typedef std::complex<double> dcomp;

double humlic(double x, double y)
{
    dcomp T(y,- x);
    double S = fabs(x) + y;
   // printf("x: %16.4e y: %16.4e\n",x,y);
    if (S >= 15) {
        // Region I
        dcomp humlic1 = T*0.5641896/(0.5+T*T);
        //fprintf(stdout, "I");
        return humlic1.real();
    }
    if (S >= 5.5) {
        // Region II
        dcomp U = T * T;
        dcomp humlic1 = T * (1.410474 + U*.5641896)/(.75 + U*(3.+U));
        //fprintf(stdout, "II");
        return humlic1.real();
    }
    if (y >= 0.195 * fabs(x) - 0.176) {
        // Region III
        dcomp humlic1 = (16.4955+T*(20.20933+T*(11.96482
                +T*(3.778987+T*.5642236)))) / (16.4955+T*(38.82363
                +T*(39.27121+T*(21.69274+T*(6.699398+T)))));
        //fprintf(stdout, "III");
        return humlic1.real();
	}
	else {
		// Region IV
		dcomp U = T * T;
		dcomp humlic1;
		humlic1 = exp(U) - T*(36183.31 - U*(3321.9905 - U*(1540.787 - U*(219.0313 - U*
			(35.76683 - U*(1.320522 - U*.56419)))))) / (32066.6 - U*(24322.84 - U*
			(9022.228 - U*(2186.181 - U*(364.2191 - U*(61.57037 - U*(1.841439 - U)))))));
		//fprintf(stdout, "IV");
		return humlic1.real();
	}
	
	
}

#include <math.h>
#include <complex.h>

#define RT2LN2 1.1774100225154747     // sqrt(2ln2)
#define RTPI 1.7724538509055159       // sqrt(pi)
#define RT2PI 2.5066282746310002      // sqrt(2pi)
#define RT2 1.4142135623730951        // sqrt(pi)

double humlic(double x, double y)
{
    double complex T = y - x*I;
    double S = fabs(x) + y;
    if (S >= 15) {
        // Region I
        double complex humlic1 = T*0.5641896/(0.5+T*T);
        //fprintf(stdout, "I");
        return creal(humlic1);
    }
    if (S >= 5.5) {
        // Region II
        double complex U = T * T;
        double complex humlic1 = T * (1.410474 + U*.5641896)/(.75 + U*(3.+U));
        //fprintf(stdout, "II");
        return creal(humlic1);
    }
    if (y >= 0.195 * fabs(x) - 0.176) {
        // Region III
        double complex humlic1 = (16.4955+T*(20.20933+T*(11.96482
                +T*(3.778987+T*.5642236)))) / (16.4955+T*(38.82363
                +T*(39.27121+T*(21.69274+T*(6.699398+T)))));
        //fprintf(stdout, "III");
        return creal(humlic1);
    }
    // Region IV
    double complex U = T * T;
    double complex humlic1;
    humlic1 = cexp(U)-T*(36183.31-U*(3321.9905-U*(1540.787-U*(219.0313-U*
       (35.76683-U*(1.320522-U*.56419))))))/(32066.6-U*(24322.84-U*
       (9022.228-U*(2186.181-U*(364.2191-U*(61.57037-U*(1.841439-U)))))));
    //fprintf(stdout, "IV");
    return creal(humlic1);
}

double voigt_humlicek(double nu, double nu0, double alpha, double gamma)
{
    double sigma = alpha / RT2LN2;
    double x = (nu0 - nu) / sigma / RT2;
    double y = gamma / sigma / RT2;
    return humlic(x, y) / sigma / RT2PI;
}



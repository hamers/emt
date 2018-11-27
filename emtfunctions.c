#include <stdio.h>
#include <math.h>



/* Functions that appear in the right-hand-sides of the
 * orbit-averaged equations of motion for the `emt' model. 
 * Interfaced to Python using ctypes (uses wrapperemtlibrary.py)
 * Adrian Hamers
 * November 2018 */

#define c_1div3 (1.0/3.0)

double fm(double e, double x, double E0, double Etau);
double fa(double e, double x, double E0, double Etau);
double fe(double e, double x, double E0, double Etau);
double fomega(double e, double x, double E0, double Etau);
double ga(double e, double x, double E0);
double ge(double e, double x, double E0);
double ha(double e, double x, double E0);
double he(double e, double x, double E0);
double compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity);
double XL0_q(double q);

double R_Lc_div_a(double q);

/* Roche-lobe radius for circular orbit */
double R_Lc_div_a(double q)
/* 1983ApJ...268..368E
    q is defined as m_primary/m_secondary */
{
    double q_pow_one_third = pow(q,c_1div3);
    double q_pow_two_third = q_pow_one_third*q_pow_one_third;
    return 0.49*q_pow_two_third/(0.6*q_pow_two_third + log(1.0 + q_pow_one_third));
}

double fm(double e, double x, double E0, double Etau)
{
    return -(-192*E0 + 576*E0*x - 576*E0*pow(x,2) - 288*pow(e,2)*E0*pow(x,2) + 192*E0*pow(x,3) + 288*pow(e,2)*E0*pow(x,3) + 72*pow(e,2)*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) - 
      96*e*(-1 + x)*(2 - 4*x + (2 + 3*pow(e,2))*pow(x,2))*sin(E0) + 6*pow(e,4)*pow(x,3)*sin(2*E0 - 3*Etau) - 8*pow(e,3)*pow(x,3)*sin(3*E0 - 3*Etau) + 
      3*pow(e,4)*pow(x,3)*sin(4*E0 - 3*Etau) + 72*pow(e,3)*pow(x,2)*sin(E0 - 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - 2*Etau) - 72*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) + 
      72*pow(e,2)*pow(x,3)*sin(2*E0 - 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 - 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 - 2*Etau) - 288*e*x*sin(E0 - Etau) + 
      576*e*pow(x,2)*sin(E0 - Etau) - 288*e*pow(x,3)*sin(E0 - Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - Etau) + 72*pow(e,2)*x*sin(2*E0 - Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 - Etau) + 
      72*pow(e,2)*pow(x,3)*sin(2*E0 - Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 - Etau) - 288*e*x*sin(E0 + Etau) + 576*e*pow(x,2)*sin(E0 + Etau) - 288*e*pow(x,3)*sin(E0 + Etau) - 
      72*pow(e,3)*pow(x,3)*sin(E0 + Etau) - 72*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) + 72*pow(e,2)*pow(x,3)*sin(2*(E0 + Etau)) - 8*pow(e,3)*pow(x,3)*sin(3*(E0 + Etau)) + 
      72*pow(e,2)*x*sin(2*E0 + Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 + Etau) + 72*pow(e,2)*pow(x,3)*sin(2*E0 + Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 + Etau) + 
      72*pow(e,3)*pow(x,2)*sin(E0 + 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 + 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 + 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 + 2*Etau) + 
      6*pow(e,4)*pow(x,3)*sin(2*E0 + 3*Etau) + 3*pow(e,4)*pow(x,3)*sin(4*E0 + 3*Etau))/(192.*M_PI);
}

double fa(double e, double x, double E0, double Etau)
{
    return (192*E0 - 576*E0*x + 576*E0*pow(x,2) + 288*pow(e,2)*E0*pow(x,2) - 192*E0*pow(x,3) - 288*pow(e,2)*E0*pow(x,3) + 72*pow(e,2)*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) - 
     96*e*(-1 + x)*(2 - 4*x + (2 + 3*pow(e,2))*pow(x,2))*sin(E0) + 6*pow(e,4)*pow(x,3)*sin(2*E0 - 3*Etau) + 8*pow(e,3)*pow(x,3)*sin(3*E0 - 3*Etau) + 
     3*pow(e,4)*pow(x,3)*sin(4*E0 - 3*Etau) + 72*pow(e,3)*pow(x,2)*sin(E0 - 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 - 2*Etau) + 72*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) - 
     72*pow(e,2)*pow(x,3)*sin(2*E0 - 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 - 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 - 2*Etau) + 288*e*x*sin(E0 - Etau) - 576*e*pow(x,2)*sin(E0 - Etau) + 
     288*e*pow(x,3)*sin(E0 - Etau) + 72*pow(e,3)*pow(x,3)*sin(E0 - Etau) + 72*pow(e,2)*x*sin(2*E0 - Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 - Etau) + 
     72*pow(e,2)*pow(x,3)*sin(2*E0 - Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 - Etau) + 288*e*x*sin(E0 + Etau) - 576*e*pow(x,2)*sin(E0 + Etau) + 288*e*pow(x,3)*sin(E0 + Etau) + 
     72*pow(e,3)*pow(x,3)*sin(E0 + Etau) + 72*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) - 72*pow(e,2)*pow(x,3)*sin(2*(E0 + Etau)) + 8*pow(e,3)*pow(x,3)*sin(3*(E0 + Etau)) + 
     72*pow(e,2)*x*sin(2*E0 + Etau) - 144*pow(e,2)*pow(x,2)*sin(2*E0 + Etau) + 72*pow(e,2)*pow(x,3)*sin(2*E0 + Etau) + 18*pow(e,4)*pow(x,3)*sin(2*E0 + Etau) + 
     72*pow(e,3)*pow(x,2)*sin(E0 + 2*Etau) - 72*pow(e,3)*pow(x,3)*sin(E0 + 2*Etau) + 24*pow(e,3)*pow(x,2)*sin(3*E0 + 2*Etau) - 24*pow(e,3)*pow(x,3)*sin(3*E0 + 2*Etau) + 
     6*pow(e,4)*pow(x,3)*sin(2*E0 + 3*Etau) + 3*pow(e,4)*pow(x,3)*sin(4*E0 + 3*Etau))/(192.*M_PI);
}

double fe(double e, double x, double E0, double Etau)
{
    return -((-1 + pow(e,2))*(12*e*E0*x*(4 - 8*x + (4 + pow(e,2))*pow(x,2))*cos(Etau) + 
        (32 - 96*x + 96*pow(x,2) + 48*pow(e,2)*pow(x,2) - 32*pow(x,3) - 48*pow(e,2)*pow(x,3) + 3*pow(e,3)*pow(x,3)*cos(E0 - 3*Etau) + pow(e,3)*pow(x,3)*cos(3*E0 - 3*Etau) + 
           8*pow(e,2)*pow(x,2)*cos(2*E0 - 2*Etau) - 8*pow(e,2)*pow(x,3)*cos(2*E0 - 2*Etau) + 24*e*x*cos(E0 - Etau) - 48*e*pow(x,2)*cos(E0 - Etau) + 24*e*pow(x,3)*cos(E0 - Etau) + 
           6*pow(e,3)*pow(x,3)*cos(E0 - Etau) + 32*pow(e,2)*pow(x,2)*cos(2*Etau) - 32*pow(e,2)*pow(x,3)*cos(2*Etau) + 24*e*x*cos(E0 + Etau) - 48*e*pow(x,2)*cos(E0 + Etau) + 
           24*e*pow(x,3)*cos(E0 + Etau) + 6*pow(e,3)*pow(x,3)*cos(E0 + Etau) + 8*pow(e,2)*pow(x,2)*cos(2*(E0 + Etau)) - 8*pow(e,2)*pow(x,3)*cos(2*(E0 + Etau)) + 
           pow(e,3)*pow(x,3)*cos(3*(E0 + Etau)) + 3*pow(e,3)*pow(x,3)*cos(E0 + 3*Etau))*sin(E0)))/(32.*M_PI);
}
double fomega(double e, double x, double E0, double Etau)
{
    return -(sqrt(1 - pow(e,2))*x*sin(Etau)*(-48*E0 + 96*E0*x - 48*E0*pow(x,2) - 12*pow(e,2)*E0*pow(x,2) + 4*(6 - 12*x + (6 + pow(e,2))*pow(x,2))*sin(2*E0) + pow(e,2)*pow(x,2)*sin(4*E0) - 
        2*pow(e,2)*pow(x,2)*sin(2*E0 - 2*Etau) + pow(e,2)*pow(x,2)*sin(4*E0 - 2*Etau) - 24*e*x*sin(E0 - Etau) + 24*e*pow(x,2)*sin(E0 - Etau) + 8*e*x*sin(3*E0 - Etau) - 
        8*e*pow(x,2)*sin(3*E0 - Etau) - 24*e*x*sin(E0 + Etau) + 24*e*pow(x,2)*sin(E0 + Etau) - 2*pow(e,2)*pow(x,2)*sin(2*(E0 + Etau)) + 8*e*x*sin(3*E0 + Etau) - 
        8*e*pow(x,2)*sin(3*E0 + Etau) + pow(e,2)*pow(x,2)*sin(4*E0 + 2*Etau)))/(32.*M_PI);
}

double ga(double e, double x, double E0)
{
    return (4*E0*x*(-8*(3 + (-3 + x)*x) + pow(e,2)*(12 + (-8 + pow(e,2))*pow(x,2))) + 64*sqrt(1 - pow(e,2))*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) + 
     e*x*(-16*(6 + pow(e,2)*(-3 + x) - 4*x)*x*sin(E0) + e*(8*(3 + x*(-6 + (2 + pow(e,2))*x))*sin(2*E0) + e*x*(-16*(-1 + x)*sin(3*E0) + 3*e*x*sin(4*E0)))))/(32.*M_PI);
}
double ge(double e, double x, double E0)
{
    return -((-1 + pow(e,2))*(12*E0*(-2 + pow(e,2)*x*(6 + x*(-9 + (4 + pow(e,2))*x))) + 48*sqrt(1 - pow(e,2))*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) - 
        6*e*(-4 + x*(24 + 8*(-3 + x)*x + pow(e,2)*x*(-15 + 14*x)))*sin(E0) + pow(e,2)*x*(6*(6 + x*(-15 + 2*(4 + pow(e,2))*x))*sin(2*E0) + e*x*(2*(9 - 10*x)*sin(3*E0) + 3*e*x*sin(4*E0))))
      )/(48.*e*M_PI);
}
double ha(double e, double x, double E0)
{
    return ((8*atan(sqrt(-((1 + e)/(-1 + e)))*tan(E0/2.)))/sqrt(1 - pow(e,2)) + e*(-12*x - (-4 + pow(e,2))*pow(x,3) - 4/(-1 + e*cos(E0)))*sin(E0) + 
     x*(-2*E0*(6 + x*(-6 + (2 + pow(e,2))*x)) - pow(e,2)*x*(-3*(-2 + x)*sin(2*E0) + e*x*sin(3*E0))))/(4.*M_PI);
}
double he(double e, double x, double E0)
{
    return (144*pow(1 - pow(e,2),1.5)*x*atan(sqrt(-((1 + e)/(-1 + e)))*tan(E0/2.)) + 48*sqrt(1 - pow(e,2))*(1 + 3*(-1 + pow(e,2))*x)*atan(((1 + e)*tan(E0/2.))/sqrt(1 - pow(e,2))) + 
     ((-1 + pow(e,2))*(-24*E0 - 36*pow(e,2)*E0*pow(x,2) + 24*pow(e,2)*E0*pow(x,3) - 12*e*E0*(-2 + pow(e,2)*pow(x,2)*(-3 + 2*x))*cos(E0) - 
          3*e*(-8 + 48*x - 3*(16 + 3*pow(e,2))*pow(x,2) + 2*(8 + 7*pow(e,2))*pow(x,3))*sin(E0) + 72*pow(e,2)*x*sin(2*E0) - 126*pow(e,2)*pow(x,2)*sin(2*E0) + 
          60*pow(e,2)*pow(x,3)*sin(2*E0) + 16*pow(e,4)*pow(x,3)*sin(2*E0) + 27*pow(e,3)*pow(x,2)*sin(3*E0) - 26*pow(e,3)*pow(x,3)*sin(3*E0) + 4*pow(e,4)*pow(x,3)*sin(4*E0))
        )/(-1 + e*cos(E0)))/(48.*e*M_PI);
}


/* Solve Kepler equation to get E_\tau as a function of \tau */
double compute_eccentric_anomaly_from_mean_anomaly(double mean_anomaly, double eccentricity)
{
    double eccentric_anomaly;
    double eccentric_anomaly_next = mean_anomaly;
    double epsilon = 1e-10;
    double error = 2.0*epsilon;
    int j = 0;
    while (error > epsilon || j < 15)
    {
        j += 1;
        eccentric_anomaly = eccentric_anomaly_next;
        eccentric_anomaly_next = eccentric_anomaly - (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly)/(1.0 - eccentricity*cos(eccentric_anomaly));
        error = fabs(eccentric_anomaly_next - eccentric_anomaly);
    }
    return eccentric_anomaly;
}

/* XL0(q) */
double XL0_q(double q)
{
    double Aplus = pow(q*(54 + pow(q,2) + 6*sqrt(3)*sqrt(27 + pow(q,2))),0.3333333333333333);
    double Aminus = pow(q*(54 + pow(q,2) - 6*sqrt(3)*sqrt(27 + pow(q,2))),0.3333333333333333);
    return (3 + sqrt(3)*sqrt(3 + Aminus + Aplus - 2*q) - sqrt(3)*sqrt(6 - Aplus - 4*q - pow(q,2)/Aplus + (6*sqrt(3)*(1 + q))/sqrt(3 + Aminus + Aplus - 2*q)))/6.;
}

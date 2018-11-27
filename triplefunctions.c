#include <stdio.h>
#include <math.h>

#define c_1div3 (1.0/3.0)
#define c_1div16 (1.0/16.0)
#define c_15div16 (15.0/16.0)
#define c_1div4 (1.0/4.0)
#define c_5div2 (5.0/2.0)

void triple_EOM_(double CONST_G, double CONST_C, double e_in, double e_out, double g_in, double g_out, double cositot, \
    double a_in, double a_out, double m1, double m2, double m3, \
    int include_quadrupole_terms, int include_octupole_terms, int include_1PN_terms, \
    double *de_in_dt, double *de_out_dt, double *dg_in_dt, double *dg_out_dt, double *dcositot_dt);


void triple_EOM_(double CONST_G, double CONST_C, double e_in, double e_out, double g_in, double g_out, double cositot, \
    double a_in, double a_out, double m1, double m2, double m3, \
    int include_quadrupole_terms, int include_octupole_terms, int include_1PN_terms, \
    double *de_in_dt, double *de_out_dt, double *dg_in_dt, double *dg_out_dt, double *dcositot_dt)
{
    
     
	/*	mass quantities	*/
	double m1_div_m2 = m1/m2;
	double m2_div_m1 = 1.0/m1_div_m2;    
    double m1_plus_m2 = m1+m2;
    double m1_plus_m2_plus_m3 = m1_plus_m2+m3;
    double m1_times_m2 = m1*m2;
    double m1_times_m2_times_m3 = m1_times_m2*m3;
    double m1_plus_m2_div_m3 = m1_plus_m2/m3;
    double m1_times_m1 = m1*m1;
    double m2_times_m2 = m2*m2;
        
	/*	eccentricity quantities	*/
	double e_in_p2 = e_in*e_in;
	double l_in_p2 = 1.0 - e_in_p2;
    double l_in = sqrt(l_in_p2);
    double l_in_p3 = l_in_p2*l_in;
    double l_in_p4 = l_in_p3*l_in;
    double l_in_p5 = l_in_p4*l_in;
    double l_in_p6 = l_in_p5*l_in;
    double l_in_p7 = l_in_p6*l_in;
    
	double e_out_p2 = e_out*e_out;
	double l_out_p2 = 1.0 - e_out_p2;
    double l_out = sqrt(l_out_p2);
    double l_out_p3 = l_out_p2*l_out;
    double l_out_p4 = l_out_p3*l_out;
    double l_out_p5 = l_out_p4*l_out;    
    double l_out_p6 = l_out_p5*l_out;
    double l_out_p7 = l_out_p6*l_out;
    
	/*	triple secular gravitational dynamics quantities */
    /* 2000ApJ...535..385F */
	double L_in = m1_times_m2*sqrt(CONST_G*a_in/(m1_plus_m2));
	double L_out = (m1_plus_m2)*m3*sqrt(CONST_G*a_out/(m1_plus_m2_plus_m3));
	double G_in = L_in*sqrt(l_in_p2);
	double G_out = L_out*sqrt(l_out_p2);
    double G_tot = sqrt( G_in*G_in + G_out*G_out + 2.0*G_in*G_out*cositot );

	double a_in_div_a_out = a_in/a_out;
    double C2,C3;
    if (include_quadrupole_terms == 0)
    {
        C2 = 0.0;
    }
    else
    {
        C2 = CONST_G*c_1div16*(m1_times_m2_times_m3/m1_plus_m2)*pow(l_out,-3.0)*a_in_div_a_out*a_in_div_a_out/a_out;
    }
    if (include_octupole_terms == 0)
    {
        C3 = 0.0;
    }
    else
    {
        C3 = -CONST_G*c_15div16*c_1div4*(m1_times_m2_times_m3/(m1_plus_m2*m1_plus_m2))*(m1-m2)*pow(l_out,-5.0)*a_in_div_a_out*a_in_div_a_out*a_in_div_a_out/a_out;
    }

    
	if (cositot > 1.0)
	{
		cositot = 2.0 - cositot;
	}

	if (cositot < -1.0)
	{
		cositot = -2.0 - cositot;
	}
    //printf("test %g %g %g\n",C2,C3,cositot);

	double cositot_p2 = cositot*cositot;
	double sinitot = sqrt(1.0 - cositot_p2); // NOTE: 0 < itot < PI, so sinitot > 0 always
	double sinitot_p2 = sinitot*sinitot;

	double sin_g_in = sin(g_in);
	double sin_2g_in = sin(2.0*g_in);
	double sin_g_out = sin(g_out);
	double cos_g_in = cos(g_in);
	double cos_2g_in = cos(2.0*g_in);
	double cos_g_out = cos(g_out);
	
	/*	required for octupole-order terms	*/
	double B = 2.0 + 5.0*e_in_p2 - 7.0*e_in_p2*cos_2g_in;
	double A = 4.0 + 3.0*e_in_p2 - c_5div2*B*sinitot_p2;
	double cosphi = -cos_g_in*cos_g_out - cositot*sin_g_in*sin_g_out;


    /************************************************
     * the calculations of the ODE right-hand-sides *
     * **********************************************/
     	

    /*******************************
     * e_in_dot                    *
     * *****************************/
    
	double e_in_dot_newtonian = 0.0;

    /* Newtonian point particle -- up and including octupole order */    
     e_in_dot_newtonian = C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in) \
            + C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in \
                - 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2 \
                - A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out));
    
    *de_in_dt = e_in_dot_newtonian;


    /*******************************
     * e_out_dot                       *
     * *****************************/

    double e_out_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    e_out_dot_newtonian = -C3*e_in*(l_out_p2/G_out)*(10.0*cositot*sinitot_p2*l_in_p2*sin_g_in*cos_g_out \
        + A*(cos_g_in*sin_g_out - cositot*sin_g_in*cos_g_out));
            
    /* combined */
	*de_out_dt = e_out_dot_newtonian;



    /*******************************
     * g_in_dot                    *
     * *****************************/

	double g_in_dot_newtonian = 0.0;
    double g_in_dot_GR_1PN_in = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */

    g_in_dot_newtonian = 6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) \
                + (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) \
            - C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2 \
                - 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) \
                - (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) \
                + cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));
    if (include_1PN_terms == 1)
    {
        g_in_dot_GR_1PN_in = (3.0/(CONST_C*CONST_C*a_in*l_in_p2))*pow(CONST_G*m1_plus_m2/a_in,3.0/2.0);
    }
    
    /* combined */  
    *dg_in_dt = g_in_dot_newtonian + g_in_dot_GR_1PN_in;


    /*******************************
     * g_out_dot                   *
     * *****************************/
     
    double g_out_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    g_out_dot_newtonian = 3.0*C2*((2.0*cositot/G_in)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)) \
            + (1.0/G_out)*(4.0 + 6.0*e_in_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)))) \
        + C3*e_in*(sin_g_in*sin_g_out*(((4.0*e_out_p2 + 1.0)/(e_out*G_out))*10.0*cositot*sinitot_p2*l_in_p2 \
            - e_out*((1.0/G_in) + (cositot/G_out))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l_in_p2)) \
            + cosphi*(5.0*B*cositot*e_out*((1.0/G_in) + (cositot/G_out)) + ((4.0*e_out_p2 + 1.0)/(e_out*G_out))*A));
    
    /* combined */
    *dg_out_dt = g_out_dot_newtonian;
	


    /**********************************************
     * cositot_dot                                *
     * due to dynamical triple interaction only!  *
     * ********************************************/

    //double cositot_dot = 0.0;
    double G_in_dot = -G_in*e_in*e_in_dot_newtonian/l_in_p2;
    double G_out_dot = -G_out*e_out*e_out_dot_newtonian/l_out_p2;
    
    *dcositot_dt = (-1.0/(G_in*G_out))*(G_in_dot*(G_in + G_out*cositot) + G_out_dot*(G_out + G_in*cositot));

    //printf("dcos %g %g\n",*dcositot_dt,e_in_dot_newtonian);
}

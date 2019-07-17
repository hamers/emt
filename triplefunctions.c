#include <stdio.h>
#include <math.h>

#define c_1div3 (1.0/3.0)
#define c_1div4 (1.0/4.0)
#define c_1div8 (1.0/8.0)
#define c_1div16 (1.0/16.0)
#define c_3div2 (3.0/2.0)
#define c_3div8 (3.0/8.0)
#define c_5div2 (5.0/2.0)
#define c_5div16 (5.0/16.0)
#define c_5div64 (5.0/64.0)
#define c_11div18 (11.0/18.0)
#define c_15div2 (15.0/2.0)
#define c_15div4 (15.0/4.0)
#define c_15div8 (15.0/8.0)
#define c_15div16 (15.0/16.0)
#define c_25div64 (25.0/64.0)
#define c_31div2 (31.0/2.0)
#define c_45div8 (45.0/8.0)
#define c_185div16 (185.0/16.0)
#define c_255div8 (255.0/8.0)



void triple_EOM_(double CONST_G, double CONST_C, double e_in, double e_out, double g_in, double g_out, double cositot, \
    double a_in, double a_out, double m1, double m2, double m3, \
    double spin_angular_frequency1, double spin_angular_frequency2, double R1, double R2, double k_div_T_tides_star1, double k_div_T_tides_star2, \
    double gyration_radius_star1, double gyration_radius_star2, double apsidal_motion_constant_star1, double apsidal_motion_constant_star2, \
    int include_quadrupole_terms, int include_octupole_terms, int include_1PN_terms, int include_inner_tidal_terms, \
    double *da_in_dt, double *de_in_dt, double *de_out_dt, double *dg_in_dt, double *dg_out_dt, double *dcositot_dt, double *spin_angular_frequency1_dot, double *spin_angular_frequency2_dot);
double f_tides1(double e_p2);
double f_tides2(double e_p2);
double f_tides3(double e_p2);
double f_tides4(double e_p2);
double f_tides5(double e_p2);

void triple_EOM_(double CONST_G, double CONST_C, double e_in, double e_out, double g_in, double g_out, double cositot, \
    double a_in, double a_out, double m1, double m2, double m3, \
    double spin_angular_frequency1, double spin_angular_frequency2, double R1, double R2, double k_div_T_tides_star1, double k_div_T_tides_star2, \
    double gyration_radius_star1, double gyration_radius_star2, double apsidal_motion_constant_star1, double apsidal_motion_constant_star2, \
    int include_quadrupole_terms, int include_octupole_terms, int include_1PN_terms, int include_inner_tidal_terms, \
    double *da_in_dt, double *de_in_dt, double *de_out_dt, double *dg_in_dt, double *dg_out_dt, double *dcositot_dt, double *spin_angular_frequency1_dot, double *spin_angular_frequency2_dot)
{
   
    //printf("=========\n");
    //printf("test pars %g %g %g %g %g %g %g\n",CONST_G,CONST_C,e_in,e_out,g_in,g_out,cositot);
    //printf("test pars %g %g %g %g %g\n",a_in,a_out,m1,m2,m3);
    //printf("test pars %g %g %g %g %g %g\n",spin_angular_frequency1,spin_angular_frequency2,R1,R2,k_div_T_tides_star1,k_div_T_tides_star2);
    //printf("test pars %g %g %g %g \n",gyration_radius_star1,gyration_radius_star2,apsidal_motion_constant_star1,apsidal_motion_constant_star2);
    //printf("test pars %d %d %d %d \n",include_quadrupole_terms,include_octupole_terms,include_1PN_terms,include_inner_tidal_terms);
    
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

    int ignore_tertiary = 0;
    if (include_quadrupole_terms == 0 && include_octupole_terms == 0)
    {
        ignore_tertiary = 1;
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


    /* tides quantities */
	double R1_div_a_in = R1/a_in;
	double R1_div_a_in_p2 = R1_div_a_in*R1_div_a_in;
	double R1_div_a_in_p5 = R1_div_a_in_p2*R1_div_a_in_p2*R1_div_a_in;
	double R1_div_a_in_p6 = R1_div_a_in*R1_div_a_in_p5;
	double R2_div_a_in = R2/a_in;
	double R2_div_a_in_p2 = R2_div_a_in*R2_div_a_in;
	double R2_div_a_in_p5 = R2_div_a_in_p2*R2_div_a_in_p2*R2_div_a_in;
	double R2_div_a_in_p6 = R2_div_a_in*R2_div_a_in_p5;
	
	double m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1 = m2_div_m1*R1_div_a_in_p6*k_div_T_tides_star1;
	double m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2 = m1_div_m2*R2_div_a_in_p6*k_div_T_tides_star2;
    
    double a_in_p2 = a_in*a_in;
    double a_in_p3 = a_in_p2*a_in;
	double n_in = sqrt(CONST_G*m1_plus_m2/a_in_p3); // mean orbital angular speed
    double spin_angular_frequency1_div_n_in = spin_angular_frequency1/n_in;
    double spin_angular_frequency2_div_n_in = spin_angular_frequency2/n_in;

    double f_tides1_in = f_tides1(e_in_p2);
    double f_tides2_in = f_tides2(e_in_p2);
    double f_tides3_in = f_tides3(e_in_p2);
    double f_tides4_in = f_tides4(e_in_p2);
    double f_tides5_in = f_tides5(e_in_p2);


    /************************************************
     * the calculations of the ODE right-hand-sides *
     * **********************************************/
     	
    /*******************************
     * a_in_dot                    *
     * *****************************/
     
    double a_in_dot_tides = 0.0;
    /* tides */
    if (include_inner_tidal_terms == 1)
    {
        double a_in_dot_tides_star1 = -6.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency1_div_n_in);
        double a_in_dot_tides_star2 = -6.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*a_in*pow(l_in,-15.0)*(f_tides1_in \
            - l_in_p3*f_tides2_in*spin_angular_frequency2_div_n_in);            
        a_in_dot_tides = a_in_dot_tides_star1 + a_in_dot_tides_star2;
    }
    *da_in_dt = a_in_dot_tides; 


    /*******************************
     * e_in_dot                    *
     * *****************************/
    
	double e_in_dot_newtonian = 0.0;
    double e_in_dot_tides = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */    
    if (ignore_tertiary == 0)
    {
        e_in_dot_newtonian = C2*(l_in_p2/G_in)*(30.0*e_in*sinitot_p2*sin_2g_in) \
            + C3*e_out*(l_in_p2/G_in)*(35.0*cosphi*sinitot_p2*e_in_p2*sin_2g_in \
                - 10.0*cositot*sinitot_p2*cos_g_in*sin_g_out*l_in_p2 \
                - A*(sin_g_in*cos_g_out - cositot*cos_g_in*sin_g_out));
    }
    
    /* tides */
    if (include_inner_tidal_terms == 1)
    {
        double e_in_dot_tides_star1 = -27.0*(1.0+m2_div_m1)*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*R1_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency1_div_n_in);
        double e_in_dot_tides_star2 = -27.0*(1.0+m1_div_m2)*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*R2_div_a_in_p2*e_in*pow(l_in,-13.0)*(f_tides3_in \
            - c_11div18*l_in_p3*f_tides4_in*spin_angular_frequency2_div_n_in);            

        e_in_dot_tides = e_in_dot_tides_star1 + e_in_dot_tides_star2;
    }

    /* combined */
    *de_in_dt = e_in_dot_newtonian + e_in_dot_tides;

    /*******************************
     * e_out_dot                       *
     * *****************************/

    double e_out_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    if (ignore_tertiary == 0)
    {
        e_out_dot_newtonian = -C3*e_in*(l_out_p2/G_out)*(10.0*cositot*sinitot_p2*l_in_p2*sin_g_in*cos_g_out \
        + A*(cos_g_in*sin_g_out - cositot*sin_g_in*cos_g_out));
    }
    
    /* combined */
	*de_out_dt = e_out_dot_newtonian;



    /*******************************
     * g_in_dot                    *
     * *****************************/

	double g_in_dot_newtonian = 0.0;
    double g_in_dot_GR_1PN_in = 0.0;
    double g_in_dot_tides = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    if (ignore_tertiary == 0)
    {
        g_in_dot_newtonian = 6.0*C2*((1.0/G_in)*(4.0*cositot_p2 + (5.0*cos_2g_in - 1.0)*(l_in_p2 - cositot_p2)) \
                + (cositot/G_out)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in))) \
            - C3*e_out*(e_in*((1.0/G_out) + (cositot/G_in))*(sin_g_in*sin_g_out*(10.0*(3.0*cositot_p2 \
                - 1.0)*(1.0 - e_in_p2) + A) - 5.0*B*cositot*cosphi) \
                - (l_in_p2/(e_in*G_in))*(sin_g_in*sin_g_out*10.0*cositot*sinitot_p2*(1.0 - 3.0*e_in_p2) \
                + cosphi*(3.0*A - 10.0*cositot_p2 + 2.0)));
    }
    
    /* GR */
    if (include_1PN_terms == 1)
    {
        g_in_dot_GR_1PN_in = (3.0/(CONST_C*CONST_C*a_in*l_in_p2))*pow(CONST_G*m1_plus_m2/a_in,3.0/2.0);
    }

    /* tides/rotation */
    if (include_inner_tidal_terms == 1)
    {    
        double g_in_dot_tides_star1 = R1_div_a_in_p5*apsidal_motion_constant_star1*(n_in/(l_in_p4))*(15.0*m2_div_m1*f_tides4_in/l_in_p6 \
		+ (1.0+m2_div_m1)*spin_angular_frequency1_div_n_in*spin_angular_frequency1_div_n_in);
        double g_in_dot_tides_star2 = R2_div_a_in_p5*apsidal_motion_constant_star2*(n_in/(l_in_p4))*(15.0*m1_div_m2*f_tides4_in/l_in_p6 \
		+ (1.0+m1_div_m2)*spin_angular_frequency2_div_n_in*spin_angular_frequency2_div_n_in);
        g_in_dot_tides = g_in_dot_tides_star1 + g_in_dot_tides_star2;
    }
    

    /* combined */  
    *dg_in_dt = g_in_dot_newtonian + g_in_dot_GR_1PN_in + g_in_dot_tides;


    /*******************************
     * g_out_dot                   *
     * *****************************/
     
    double g_out_dot_newtonian = 0.0;
    
    /* Newtonian point particle -- up and including octupole order */
    if (ignore_tertiary == 0)
    {
        g_out_dot_newtonian = 3.0*C2*((2.0*cositot/G_in)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)) \
            + (1.0/G_out)*(4.0 + 6.0*e_in_p2 + (5.0*cositot_p2 - 3.0)*(2.0 + e_in_p2*(3.0 - 5.0*cos_2g_in)))) \
        + C3*e_in*(sin_g_in*sin_g_out*(((4.0*e_out_p2 + 1.0)/(e_out*G_out))*10.0*cositot*sinitot_p2*l_in_p2 \
            - e_out*((1.0/G_in) + (cositot/G_out))*(A + 10.0*(3.0*cositot_p2 - 1.0)*l_in_p2)) \
            + cosphi*(5.0*B*cositot*e_out*((1.0/G_in) + (cositot/G_out)) + ((4.0*e_out_p2 + 1.0)/(e_out*G_out))*A));
    }
    
    /* combined */
    *dg_out_dt = g_out_dot_newtonian;
	


    /**********************************************
     * cositot_dot                                *
     * due to dynamical triple interaction only!  *
     * ********************************************/

    //double cositot_dot = 0.0;
    if (ignore_tertiary == 0)
    {
        double G_in_dot = -G_in*e_in*e_in_dot_newtonian/l_in_p2;
        double G_out_dot = -G_out*e_out*e_out_dot_newtonian/l_out_p2;
    
        *dcositot_dt = (-1.0/(G_in*G_out))*(G_in_dot*(G_in + G_out*cositot) + G_out_dot*(G_out + G_in*cositot));
    }
    else
    {
        *dcositot_dt = 0.0;
    }
    
    //printf("dcos %g %g\n",*dcositot_dt,e_in_dot_newtonian);
    
    
    /************************************
     * spin_angular_frequency1_dot      *
     * **********************************/

	double spin_angular_frequency1_dot_tides = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == 1)
    {
        spin_angular_frequency1_dot_tides = 3.0*m2_div_m1_times_R1_div_a_in_p6_times_k_div_T_tides_star1*(m2_div_m1/(gyration_radius_star1*gyration_radius_star1)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency1_div_n_in);
    }

    *spin_angular_frequency1_dot = spin_angular_frequency1_dot_tides;



    /************************************
     * spin_angular_frequency2_dot      *
     * **********************************/

	double spin_angular_frequency2_dot_tides = 0.0;
    
    /* tides */
    if (include_inner_tidal_terms == 1)
    {
        spin_angular_frequency2_dot_tides = 3.0*m1_div_m2_times_R2_div_a_in_p6_times_k_div_T_tides_star2*(m1_div_m2/(gyration_radius_star2*gyration_radius_star2)) \
            *(n_in/(l_in_p6*l_in_p6))*(f_tides2_in - l_in_p3*f_tides5_in*spin_angular_frequency2_div_n_in);
    }
   
    *spin_angular_frequency2_dot = spin_angular_frequency2_dot_tides;
    
    
    //printf("test out %g %g %g %g %g %g %g %g\n",*da_in_dt,*de_in_dt,*de_out_dt,*dg_in_dt,*dg_out_dt,*dcositot_dt,*spin_angular_frequency1_dot,*spin_angular_frequency2_dot);
}

/* tides (1981A&A....99..126H) */
double f_tides1(double e_p2)
{
    return 1.0 + e_p2*(c_31div2 + e_p2*(c_255div8 + e_p2*(c_185div16 + e_p2*c_25div64)));
}
double f_tides2(double e_p2)
{
    return 1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16));
}
double f_tides3(double e_p2)
{
    return 1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64));
}
double f_tides4(double e_p2)
{
    return 1.0 + e_p2*(c_3div2 + e_p2*c_1div8);
}
double f_tides5(double e_p2)
{
    return 1.0 + e_p2*(3.0 + e_p2*c_3div8);
}

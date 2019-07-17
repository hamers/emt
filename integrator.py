"""
Script to integrate the orbit-averaged equations of motion. 
Requires: included `emtfunctions' library (type `make to install to current directory'), numpy, and scipy.
Testing emtfunctions: use nosetests on testemtlibrary.py: `nosetests emtlibrary.py'
Usage: "python integrator.py". Use command-line arguments to specify parameters,
or, for more flexibility, use the included `applications.py' script to define parameters
in-script.

Adrian Hamers
November 2018
"""

import argparse
import numpy as np

from scipy.integrate import odeint

from wrapperemtlibrary import emtlibrary

### Do not change the numbers below, unless you want to use different units ###
CONST_G = 4.0*np.pi**2
CONST_C = 63239.72638679138
RSun_in_AU = 0.004649130343817401

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def add_bool_arg(parser, name, default=False,help=None):
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true',help="Enable %s"%help)
    group.add_argument('--no-' + name, dest=name, action='store_false',help="Disable %s"%help)
    parser.set_defaults(**{name:default})

def parse_arguments():
    

    parser = argparse.ArgumentParser()

    parser.add_argument("--name",                           type=str,       dest="name",                        default="test01",       help="name")
    parser.add_argument("--M_d",                            type=float,     dest="M_d",                         default=1.0,            help="Donor mass [MSun]")
    parser.add_argument("--M_d_dot_av",                     type=float,     dest="M_d_dot_av",                  default=-1.0e-10,       help="Mass loss rate [MSun/yr]")
    parser.add_argument("--M_a",                            type=float,     dest="M_a",                         default=0.9,            help="Accretor mass [MSun]")    
    parser.add_argument("--a",                              type=float,     dest="a",                           default=1.0,            help="Semimajor axis [AU]")    
    parser.add_argument("--e",                              type=float,     dest="e",                           default=0.95,           help="Eccentricity")    
    parser.add_argument("--R",                              type=float,     dest="R",                           default=10.0,           help="Donor radius [RSun]")    
    parser.add_argument("--R_M_d_power",                    type=float,     dest="R_M_d_power",                 default=0.7,            help="Dimensionless number \alpha in R=RSun (M_d/MSun)^\alpha. To enable this radius set --R_from_M_d")    
    parser.add_argument("--t_end",                          type=float,     dest="t_end",                       default=0.1,            help="End time. Default: scaled in units of M_d/M_d_dot. Use --scaled_t_end to specify an absolute time, in units of yr.")    
    parser.add_argument("--model",                          type=str,       dest="model",                       default="emt",          help="Mass transfer model. Can be `emt' (`emt' model), or `sep' (Sepinsky model).")    
    parser.add_argument("--N_steps",                        type=int,       dest="N_steps",                     default=1000,           help="Number of (output) steps")    
    parser.add_argument("--mxstep",                         type=int,       dest="mxstep",                      default=100000,         help="Maximum number of internal steps taken in the ODE integratino. Increase if ODE integrator gives errors. ")    
    parser.add_argument("--rp_min",                         type=float,     dest="rp_min",                      default=0.01,           help="Minimum periapsis distance, r_p=a(1-e), below which further evolution is halted (i.e,. zero ODE right-hand-sides).")   
    parser.add_argument("--tau",                            type=float,     dest="tau",                         default=0.0,            help="Mass transfer time lag [yr]")    
    parser.add_argument("--tau_mean_anomaly",               type=float,     dest="tau_mean_anomaly",            default=0.0,            help="Mass transfer time lag expressed in mean anomaly shift [rad]. tau_mean_anomaly = n*tau, where n is the (instantaneous) mean anomaly.")    
    parser.add_argument("--ejection_radius_mode",           type=int,       dest="ejection_radius_mode",        default=1,              help="Donor ejection radius mode. 1: low donor spin limit; 2: high mass ratio limit.")    
    parser.add_argument("--donor_normalized_spin",          type=float,     dest="donor_normalized_spin",       default=1.0,            help="Donor spin parameter, \Omega_\hat = \Omega/\omega_{orb,P}: donor spin normalized to orbital angular frequency at periapsis. Only applies when ejection_radius_mode = 1.")    
    parser.add_argument("--accretion_radius",               type=float,     dest="accretion_radius",            default=1.0,            help="Accretion radius (constant) [RSun]")    
    
    parser.add_argument("--M_t",                            type=float,     dest="M_t",                         default=1.0,            help="Tertiary mass [MSun]")    
    parser.add_argument("--e_out",                          type=float,     dest="e_out",                       default=0.5,            help="Outer orbit eccentricity")    
    parser.add_argument("--a_out",                          type=float,     dest="a_out",                       default=1000,           help="Outer orbit semimajor axis [AU}")    
    parser.add_argument("--omega",                          type=float,     dest="omega",                       default=0.1,            help="Argument of periapsis [rad]")    
    parser.add_argument("--omega_out",                      type=float,     dest="omega_out",                   default=0.1,            help="Outer orbit argument of periapsis [rad]")    
    parser.add_argument("--i_rel",                          type=float,     dest="i_rel",                       default=1.0,            help="Inner-outer orbit mutual inclination [rad]")    
    
    parser.add_argument("--spin_angular_frequency_d",       type=float,     dest="spin_angular_frequency_d",    default=0.0,            help="Donor spin angular frequency used for tides [rad/yr]")    
    parser.add_argument("--spin_angular_frequency_a",       type=float,     dest="spin_angular_frequency_a",    default=0.0,            help="Accretor spin angular frequency used for tides [rad/yr]")    
    parser.add_argument("--R_a",                            type=float,     dest="R_a",                         default=1.0,            help="Accretor radius for tides [RSun]")    
    parser.add_argument("--k_div_T_tides_d",                type=float,     dest="k_div_T_tides_d",             default=0.0,            help="Donor spin angular frequency used for tides [rad/yr]")    
    parser.add_argument("--k_div_T_tides_a",                type=float,     dest="k_div_T_tides_a",             default=0.0,            help="Accretor spin angular frequency used for tides [rad/yr]")    
    parser.add_argument("--gyration_radius_d",              type=float,     dest="gyration_radius_d",           default=0.28,           help="Donor gyration radius rg used for tides [rad/yr] (I=rg^2 MR^2) ")    
    parser.add_argument("--gyration_radius_a",              type=float,     dest="gyration_radius_a",           default=0.28,           help="Accretor gyration radius rg used for tides [rad/yr] (I=rg^2 MR^2)")    
    parser.add_argument("--apsidal_motion_constant_d",      type=float,     dest="apsidal_motion_constant_d",   default=0.014,          help="Donor apsidal motion constant used for tides [rad/yr]")    
    parser.add_argument("--apsidal_motion_constant_a",      type=float,     dest="apsidal_motion_constant_a",   default=0.014,          help="Accretor apsidal motion constant used for tides [rad/yr]")    
    

    ### boolean arguments ###
    add_bool_arg(parser, 'verbose',                     default=False,          help="verbose terminal output")
    add_bool_arg(parser, 'plot',                        default=True,           help="plotting") 
    add_bool_arg(parser, 'plot_fancy',                  default=False,          help="using LaTeX for plot labels (slower).")
    add_bool_arg(parser, 'R_from_M_d',                  default=False,          help="variable donor radius according to R=RSun (M_d/MSun)^\alpha. Set the value of \alpha using --R_M_d_power. If disabled, R=constant.")
    add_bool_arg(parser, 'use_constant_tau',            default=False,          help="constant tau during integration.")
    add_bool_arg(parser, 'scaled_t_end',                default=True,           help="using scaled t_end.")
    add_bool_arg(parser, 'include_ejection_radius',     default=False,          help="using nonzero ejection radius. Specify mode with --ejection_radius_mode.")
    add_bool_arg(parser, 'include_accretion_radius',    default=False,          help="using nonzero and constant accretion radius. Specify accretion radius with --accretion_radius.")
    add_bool_arg(parser, 'include_tertiary',            default=False,          help="secular gravitational perturbation by tertiary star.")
    add_bool_arg(parser, 'include_quadrupole_terms',    default=True,           help="inclusion of secular three-body quadrupole-order terms.")
    add_bool_arg(parser, 'include_octupole_terms',      default=True,           help="inclusion of secular three-body octupole-order terms.")
    add_bool_arg(parser, 'include_1PN_terms',           default=False,          help="inclusion of first post-Newtonian terms (only applied if --include_tertiary)")
    add_bool_arg(parser, 'include_inner_tidal_terms',   default=False,          help="inclusion of inner binary equilibrium tides terms (only applied if --include_tertiary)")
    
    args = parser.parse_args()
                       
    return args
    

def determine_E_0(e,x,verbose=False):
    no_RLOF = False
    if x > 1.0/(1.0-e):
        if verbose==True:
            print 'No RLOF','x',x,'1.0/(1.0-e)',1.0/(1.0-e)
        #return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        E_0 = 0.0
        no_RLOF = True

    elif x <= 1.0/(1.0+e):
        E_0 = np.pi
        if verbose==True:
            print 'RLOF all phases','x',x,'1.0/(1.0+e)',1.0/(1.0+e)

    else:
        E_0 = np.arccos( (1.0/e)*(1.0 - 1.0/x) )
        if verbose==True:
            print 'RLOF partial orbit; E_0 = ',E_0
    return no_RLOF, E_0

def R_function(args,M_d):
    R = args.R*RSun_in_AU
    if args.R_from_M_d == True:
        R = pow(M_d,args.R_M_d_power)*RSun_in_AU
    return R

def RHS_function(RHR_vec, t, *ODE_args):
    ### initialization ###
    CONST_G,CONST_C,args = ODE_args
    verbose = args.verbose
    model = args.model

    a = RHR_vec[0]
    e = 1.0-pow(10.0,RHR_vec[1])
    omega = RHR_vec[2]

    M_d = RHR_vec[3]
    M_a = RHR_vec[4]

    e_out = 1.0-pow(10.0,RHR_vec[5])
    omega_out = RHR_vec[6]
    cos_i_rel = RHR_vec[7]
    spin_angular_frequency_d = RHR_vec[8]
    spin_angular_frequency_a = RHR_vec[9]
    
    epsilon=1.0e-12
    if e<epsilon: ### some of the emt functions give nans for exactly zero eccentricity -- replace with tiny number in this case ###
        e=epsilon

    if M_d <= epsilon or M_a <= epsilon:
        print 'No more mass left! No further evolution.'
        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

    if a*(1.0-e) <= args.rp_min:
        print 'Collision! No further evolution'
        return [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

    q = M_d/M_a

    R = R_function(args,M_d)
    R_a = args.R_a*RSun_in_AU
    
    da_dt = de_dt = domega_dt = dM_d_dt = dM_a_dt = de_out_dt = domega_out_dt = dcos_i_rel_dt = dspin_angular_frequency_d_dt = dspin_angular_frequency_a_dt = 0.0

    ### calculate tertiary object EOM ###
    if args.include_tertiary == True:
        da_in_dt_, de_in_dt_, de_out_dt_, domega_dt_, domega_out_dt_, dcos_i_rel_dt_, dspin_angular_frequency_d_dt_, dspin_angular_frequency_a_dt_ = \
            emtlibrary.triple_EOM(CONST_G,CONST_C,e,e_out,omega,omega_out,cos_i_rel, \
            a, args.a_out, M_d, M_a, args.M_t, \
            spin_angular_frequency_d, spin_angular_frequency_a, R, R_a, args.k_div_T_tides_d, args.k_div_T_tides_a, \
            args.gyration_radius_d,args.gyration_radius_a,args.apsidal_motion_constant_d,args.apsidal_motion_constant_a, \
            int(args.include_quadrupole_terms), int(args.include_octupole_terms), int(args.include_1PN_terms), int(args.include_inner_tidal_terms))
        da_dt += da_in_dt_
        de_dt += de_in_dt_
        de_out_dt += de_out_dt_
        domega_dt += domega_dt_
        domega_out_dt += domega_out_dt_
        dcos_i_rel_dt += dcos_i_rel_dt_
        dspin_angular_frequency_d_dt += dspin_angular_frequency_d_dt_
        dspin_angular_frequency_a_dt += dspin_angular_frequency_a_dt_
    
        if verbose==True:
            print 'include_tertiary',da_in_dt_, de_in_dt_, de_out_dt_, domega_dt_, domega_out_dt_, dcos_i_rel_dt_, dspin_angular_frequency_d_dt_, dspin_angular_frequency_a_dt_

        
    
    ### calculate mass transfer EOM ###
    finite_size_term_a = finite_size_term_e = 0.0
    fm=fa=fe=fomega=ga=ge=ha=he=0.0
    
    if model == "sep":
        no_RLOF = False
        fm = 1.0
        fa = np.sqrt(1.0-e**2)
        fe = fa*(1.0-e)
        fomega = 0.0
    elif model == "emt":
        R_Lc = a*emtlibrary.R_Lc_div_a(q)
        x = R_Lc/R
       
        ### determine E_0, i.e., at which parts in the orbit is there RLOF? ###
        no_RLOF, E_0 = determine_E_0(e,x,verbose=args.verbose)

        if no_RLOF == False:
            if args.use_constant_tau == True:
                n = np.sqrt(CONST_G*(M_d+M_a)/(a**3))
                MA_tau = n*args.tau
            else:
                MA_tau = args.tau_mean_anomaly

            E_tau = emtlibrary.compute_eccentric_anomaly_from_mean_anomaly(MA_tau,e)
            if E_tau > np.pi: ### E_tau should not be larger than \pi
                E_tau = np.pi
            
            fm = emtlibrary.fm(e,x,E_0,E_tau)
            fa = emtlibrary.fa(e,x,E_0,E_tau)
            fe = emtlibrary.fe(e,x,E_0,E_tau)
            fomega = emtlibrary.fomega(e,x,E_0,E_tau)

            if np.fabs(fm) <= epsilon:
                fm = epsilon

            if args.include_ejection_radius == True or args.include_accretion_radius == True:
                ha = emtlibrary.ha(e,x,E_0)
                he = emtlibrary.he(e,x,E_0)

                if args.ejection_radius_mode == 1: ### limit of small donor spin
                    ga = emtlibrary.ga(e,x,E_0)
                    ge = emtlibrary.ge(e,x,E_0)
                    XL0 = emtlibrary.XL0_q(q)
                    finite_size_term_a += XL0*ga
                    finite_size_term_e += XL0*ge
                elif args.ejection_radius_mode == 2: ### limit of large mass ratio q
                    XL0 = pow(args.donor_normalized_spin,-2.0/3.0)*(1.0-e)/pow(1.0+e,1.0/3.0)
                    finite_size_term_a += XL0*ha
                    finite_size_term_e += XL0*he
                
                if args.include_accretion_radius == True:
                    finite_size_term_a -= q*(args.accretion_radius*RSun_in_AU/a)*ha
                    finite_size_term_e -= q*(args.accretion_radius*RSun_in_AU/a)*he
            
                if verbose==True:
                    print 'include_ejection_radius',args.include_ejection_radius,'include_accretion_radius',args.include_accretion_radius,'finite_size_term_a',finite_size_term_a,'ga',ga,'ge',ge,'ha',ha,'he',he
    else:
        print 'ERROR: invalid model -- should be sep or emt'
        exit(-1)
    
    if verbose==True and model == "emt":
        print 'x',x,'a',a,'e',e,'fm',fm,'fa',fa,'fe',fe,'fomega',fomega
    

    if no_RLOF == False:
        M_d_dot_av = args.M_d_dot_av
        common_factor = -2.0*(M_d_dot_av/M_d)*(1.0/fm)

        da_dt += common_factor*a*(fa*(1.0-q) + finite_size_term_a)
        de_dt += common_factor*(fe*(1.0-q) + finite_size_term_e)
        domega_dt += common_factor*fomega*(1.0-q)

        dM_d_dt += M_d_dot_av
        dM_a_dt += -dM_d_dt


    
    ### NOTE: X = np.log10(1.0-e); Y = np.log10(1.0-e_out) ###
    dX_dt = -de_dt/(np.log(10.0)*(1.0-e))
    dY_dt = -de_out_dt/(np.log(10.0)*(1.0-e_out))

    RHR_vec_dot = [da_dt,dX_dt,domega_dt,dM_d_dt,dM_a_dt,dY_dt,domega_out_dt,dcos_i_rel_dt,dspin_angular_frequency_d_dt,dspin_angular_frequency_a_dt]

    return RHR_vec_dot

def integrate(args,return_spins=False):
    if args.verbose==True:
        print 'arguments:'
        from pprint import pprint
        pprint(vars(args))

    if args.plot_fancy == True:
        pyplot.rc('text',usetex=True)
        pyplot.rc('legend',fancybox=True)

    ### initial conditions ###   
    M_d = args.M_d
    M_a = args.M_a
    M_d_dot_av = args.M_d_dot_av
    
    a = args.a
    e = args.e
    
    e_out = args.e_out
    a_out = args.a_out
    omega = args.omega
    omega_out = args.omega_out
    i_rel = args.i_rel

    spin_angular_frequency_d = args.spin_angular_frequency_d
    spin_angular_frequency_a = args.spin_angular_frequency_a

    q = M_d/M_a

    if args.scaled_t_end == True:
        t_end = (args.t_end)*np.fabs(M_d/M_d_dot_av)
    else:
        t_end = args.t_end
    
    if args.verbose==True:
        print 't_end/yr',t_end
        
    N_steps = args.N_steps
    times = np.linspace(0.0, t_end, N_steps)

    ODE_args = (CONST_G,CONST_C,args)

    RHR_vec = [a,np.log10(1.0-e),omega,M_d,M_a,np.log10(1.0-e_out),omega_out,np.cos(i_rel),spin_angular_frequency_d,spin_angular_frequency_a]
        
    if args.verbose==True:
        print 'RHR_vec',RHR_vec
    
    ### numerical solution ###
    sol = odeint(RHS_function, RHR_vec, times, args=ODE_args,mxstep=args.mxstep)
    
    a_sol = np.array(sol[:,0])
    e_sol = 1.0 - pow(10.0,np.array(sol[:,1]))
    omega_sol = np.array(sol[:,2])
    M_d_sol = np.array(sol[:,3])
    M_a_sol = np.array(sol[:,4])
    e_out_sol = 1.0 - pow(10.0,np.array(sol[:,5]))
    omega_out_sol = np.array(sol[:,6])
    cos_i_rel_sol = np.array(sol[:,7])
    spin_angular_frequency_d_sol = np.array(sol[:,8])
    spin_angular_frequency_a_sol = np.array(sol[:,9])
        
    E_0_sol = []
    for index,e in enumerate(e_sol):
        q = M_d_sol[index]/M_a_sol[index]
        a = a_sol[index]
        R = R_function(args,M_d)
        R_Lc = a*emtlibrary.R_Lc_div_a(q)
        x = R_Lc/R
        no_RLOF, E_0 = determine_E_0(e,x,verbose=args.verbose)

        E_0_sol.append(E_0)
    E_0_sol = np.array(E_0_sol)
    
    if return_spins==False:
        return times,a_sol,e_sol,omega_sol,M_d_sol,M_a_sol,E_0_sol
    else:
        return times,a_sol,e_sol,omega_sol,M_d_sol,M_a_sol,E_0_sol,spin_angular_frequency_d_sol,spin_angular_frequency_a_sol
    
def plot_function(args,data):
    a = args.a
    e = args.e
    M_d = args.M_d
    M_a = args.M_a
    times,a_sol,e_sol,omega_sol,M_d_sol,M_a_sol,E_0_sol = data
    times*= 1.0e-6

    fontsize=18
    labelsize=12
    
    fig=pyplot.figure(figsize=(8,10))
    plot1=fig.add_subplot(3,1,1,yscale="linear")
    plot2=fig.add_subplot(3,1,2,yscale="linear")
    plot3=fig.add_subplot(3,1,3)
    
    plot1.plot(times,a_sol,color='k')
    #plot1.plot(times,a_sol*(1.0-e_sol),linestyle='dotted')
    #plot1.plot(times,a_sol*(1.0+e_sol),linestyle='dotted')
    #plot1.plot(times,a*(M_d*M_a/(M_d_sol*M_a_sol))**2,color='r',linestyle='dashed')

    plot2.plot(times,e_sol,color='k')
    plot3.plot(times,M_d_sol/M_a_sol,color='k',label="$q$")
    plot3.plot(times,E_0_sol,color='C3',linestyle='dashed',label=r"$\mathcal{E}_0$")

    plots = [plot1,plot2,plot3]
    #labels = [r"$a/\mathrm{AU}; \, r_\mathrm{p}/\mathrm{AU};\, r_\mathrm{a}/\mathrm{AU}$",r"$e$",r"$q; \, \mathcal{E}_0/\mathrm{rad}$"]
    labels = [r"$a/\mathrm{AU}$",r"$e$",r"$q; \, \mathcal{E}_0/\mathrm{rad}$"]
    for index,plot in enumerate(plots):
        if index==2:
            plot.set_xlabel(r"$t/\mathrm{Myr}$",fontsize=fontsize)
        plot.set_ylabel(labels[index],fontsize=fontsize)
                
        plot.tick_params(axis='both', which ='major', labelsize = labelsize)

        if index in [0,1]:
            plot.set_xticklabels([])
    plot3.axhline(y=1.0,linestyle='dotted',color='k')

    handles,labels = plot3.get_legend_handles_labels()
    plot3.legend(handles,labels,loc="upper right",fontsize=0.8*fontsize)
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    filename = 'elements_' + str(args.name) + '.pdf'
    fig.savefig(filename,dpi=200)
    
    pyplot.show()
    
def plot_function_tides(args,data):
    a = args.a
    e = args.e
    M_d = args.M_d
    M_a = args.M_a
    times,a_sol,e_sol,omega_sol,M_d_sol,M_a_sol,E_0_sol,spin_angular_frequency_d_sol,spin_angular_frequency_a_sol = data
    times*= 1.0e-6

    fontsize=18
    labelsize=12
    
    fig=pyplot.figure(figsize=(8,8))
    plot1=fig.add_subplot(3,1,1,yscale="linear")
    plot2=fig.add_subplot(3,1,2,yscale="linear")
    plot3=fig.add_subplot(3,1,3)
    
    plot1.plot(times,a_sol,color='k')
    plot2.plot(times,e_sol,color='k')
    plot3.plot(times,spin_angular_frequency_d_sol,color='k',label=r"$\Omega_1$")
    plot3.plot(times,spin_angular_frequency_a_sol,color='g',label=r"$\Omega_2$")

    af = a*(1.0-e**2) ### expected final semimajor axis (neglecting spins)
    nf = np.sqrt(CONST_G*(M_d+M_a)/af**3)
    plot1.axhline(y=a*af,color='r',linestyle='dashed',label="$\mathrm{Expected\,final\,}a$") 
    plot3.axhline(y=nf,color='r',linestyle='dashed',label="$\mathrm{Expected\,final\,spin}$") 
    
    plots = [plot1,plot2,plot3]
    #labels = [r"$a/\mathrm{AU}; \, r_\mathrm{p}/\mathrm{AU};\, r_\mathrm{a}/\mathrm{AU}$",r"$e$",r"$q; \, \mathcal{E}_0/\mathrm{rad}$"]
    labels = [r"$a/\mathrm{AU}$",r"$e$",r"$\Omega_i$"]
    for index,plot in enumerate(plots):
        if index==2:
            plot.set_xlabel(r"$t/\mathrm{Myr}$",fontsize=fontsize)
        plot.set_ylabel(labels[index],fontsize=fontsize)
                
        plot.tick_params(axis='both', which ='major', labelsize = labelsize)

        if index in [0,1]:
            plot.set_xticklabels([])

    handles,labels = plot1.get_legend_handles_labels()
    plot1.legend(handles,labels,loc="upper right",fontsize=0.8*fontsize)

    handles,labels = plot3.get_legend_handles_labels()
    plot3.legend(handles,labels,loc="upper right",fontsize=0.8*fontsize)
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    filename = 'elements_tides_' + str(args.name) + '.pdf'
    fig.savefig(filename,dpi=200)
    
    pyplot.show()
    
if __name__ == '__main__':
    args = parse_arguments()

    data = integrate(args)
    if args.plot == True:
        if HAS_MATPLOTLIB == False:
            print 'Error importing Matplotlib -- not making plot'
            exit(-1)
        plot_function(args,data)


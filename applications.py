import argparse
import numpy as np

#from scipy.integrate import odeint


import integrator

RSun_in_AU = 0.004649130343817401
MJ_in_MSun = 0.0009546386983890755
G = 4.0*np.pi**2


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

    parser.add_argument("--id",                 type=int,       dest="application_id",        default=0,       help="Id to identify application number. ")

    add_bool_arg(parser, 'plot_fancy',                      default=False,          help="using LaTeX for plot labels (slower).")
    add_bool_arg(parser, 'calc',                            default=True,           help="doing calculations; will be subsequently saved to disk.")
    add_bool_arg(parser, 'verbose',                         default=False,          help="verbose terminal output")
    
    args = parser.parse_args()

    return args
    

   
def plot_function_compare_sep(cmd_args,data1,data2,ylim3=None,\
    extra_A_data1=None,extra_A_data2=None, \
    extra_B_data1=None,extra_B_data2=None, \
    label_A_1=None,label_A_2=None, \
    label_B_1=None,label_B_2=None):
    if cmd_args.plot_fancy == True:
        pyplot.rc('text',usetex=True)
        pyplot.rc('legend',fancybox=True)
    
    a = args.a
    e = args.e
    M_d = args.M_d
    M_a = args.M_a
    times1,a_sol1,e_sol1,omega_sol1,M_d_sol1,M_a_sol1,E_0_sol1 = data1
    times2,a_sol2,e_sol2,omega_sol2,M_d_sol2,M_a_sol2,E_0_sol2 = data2
    times1 *= 1.0e-6
    times2 *= 1.0e-6

    fontsize=22
    labelsize=18
    linewidth=2.0
    linewidth2=3.0
    
    fig=pyplot.figure(figsize=(8,10))
    plot1=fig.add_subplot(3,1,1,yscale="linear")
    plot2=fig.add_subplot(3,1,2)
    plot3=fig.add_subplot(3,1,3)
    
    figB=pyplot.figure(figsize=(8,10))
    plotB=figB.add_subplot(3,1,1,yscale="log")

    label1="$\mathrm{emt}$"
    label2="$\mathrm{Sepinsky}$"
    
    plot1.plot(times1,a_sol1,color='k',label=label1,linewidth=linewidth)
    plot1.plot(times2,a_sol2,color='k',linestyle='dashed',label=label2,linewidth=linewidth2)
    plot1.plot(times1,a*(M_d*M_a/(M_d_sol1*M_a_sol1))**2,color='C3',linestyle='dotted',label="$\mathrm{Analytic\,(circular \,only)}$",linewidth=linewidth)
    
    plot2.plot(times1,e_sol1,color='k',linewidth=linewidth,label=label1)
    plot2.plot(times2,e_sol2,color='k',linestyle='dashed',linewidth=linewidth2,label=label2)
    plot3.plot(times1,M_d_sol1/M_a_sol1,color='k',linewidth=linewidth,label=label1)
    plot3.plot(times2,M_d_sol2/M_a_sol2,color='k',linestyle='dashed',linewidth=linewidth2,label=label2)
    plot3.plot(times1,E_0_sol1,color='C3',linestyle='dashed',label=label1)
    
    if extra_A_data1!=None and extra_A_data2!=None:
        A_times1,A_a_sol1,A_e_sol1,A_omega_sol1,A_M_d_sol1,A_M_a_sol1,A_E_0_sol1 = extra_A_data1
        A_times2,A_a_sol2,A_e_sol2,A_omega_sol2,A_M_d_sol2,A_M_a_sol2,A_E_0_sol2 = extra_A_data2
        A_times1 *= 1.0e-6
        A_times2 *= 1.0e-6
        linewidth=1.0
        plot2.plot(A_times1,A_e_sol1,color='C0',linewidth=linewidth,label=label_A_1)
        plot2.plot(A_times2,A_e_sol2,color='C0',linestyle='dashed',linewidth=linewidth,label=label_A_2)
    if extra_B_data1!=None and extra_B_data2!=None:
        B_times1,B_a_sol1,B_e_sol1,B_omega_sol1,B_M_d_sol1,B_M_a_sol1,B_E_0_sol1 = extra_B_data1
        B_times2,B_a_sol2,B_e_sol2,B_omega_sol2,B_M_d_sol2,B_M_a_sol2,B_E_0_sol2 = extra_B_data2
        B_times1 *= 1.0e-6
        B_times2 *= 1.0e-6
        linewidth=0.5
        plot2.plot(B_times1,B_e_sol1,color='C5',linewidth=linewidth,label=label_B_1)
        plot2.plot(B_times2,B_e_sol2,color='C5',linestyle='dashed',linewidth=linewidth,label=label_B_2)

    
    plotB.plot(times1,M_d_sol1*M_a_sol1*np.sqrt(a_sol1*(1.0-e_sol1**2)),color='k',label=label1,linewidth=linewidth)
    plotB.plot(times2,M_d_sol2*M_a_sol2*np.sqrt(a_sol2*(1.0-e_sol2**2)),color='k',linestyle='dashed',label=label2,linewidth=linewidth2)


    plots = [plot1,plot2,plot3]
    labels = [r"$a/\mathrm{AU}$",r"$e$",r"$q; \mathcal{E}_0/\mathrm{rad}$"]
    for index,plot in enumerate(plots):
        if index==2:
            plot.set_xlabel(r"$t/\mathrm{Myr}$",fontsize=fontsize)
        plot.set_ylabel(labels[index],fontsize=fontsize)
                
        plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

        if index in [0,1]:
            plot.set_xticklabels([])

    if ylim3 != None:
        plot3.set_ylim(0,ylim3)

    beta=0.7
    handles,labels = plot1.get_legend_handles_labels()
    plot1.legend(handles,labels,loc=cmd_args.plot1_loc,fontsize=beta*fontsize)

    handles,labels = plot2.get_legend_handles_labels()
    plot2.legend(handles,labels,loc=cmd_args.plot2_loc,fontsize=beta*fontsize)

    handles,labels = plot3.get_legend_handles_labels()
    plot3.legend(handles,labels,loc=cmd_args.plot3_loc,fontsize=beta*fontsize)

    plot3.axhline(y=1.0,linestyle='dotted',color='k')
    plot3.axhline(y=np.pi,linestyle='dotted',color='C3')
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    filename = cmd_args.filename
    fig.savefig(filename,dpi=200)
    
    pyplot.show()

def plot_function_compare_size(cmd_args,data1,data2,plot_omega=False,ylim3=None):
    if cmd_args.plot_fancy == True:
        pyplot.rc('text',usetex=True)
        pyplot.rc('legend',fancybox=True)
    
    a = args.a
    e = args.e
    M_d = args.M_d
    M_a = args.M_a
    times1,a_sol1,e_sol1,omega_sol1,M_d_sol1,M_a_sol1,E_0_sol1 = data1
    times2,a_sol2,e_sol2,omega_sol2,M_d_sol2,M_a_sol2,E_0_sol2 = data2
    times1 *= 1.0e-6
    times2 *= 1.0e-6

    figB=pyplot.figure(figsize=(8,10))
    plotB=figB.add_subplot(3,1,1,yscale="log")



    fontsize=22
    labelsize=18
    linewidth=2.0
    linewidth2=3.0
    
    fig=pyplot.figure(figsize=(8,10))
    plot1=fig.add_subplot(3,1,1,yscale="log")
    plot2=fig.add_subplot(3,1,2)
    plot3=fig.add_subplot(3,1,3)
    
    label1=cmd_args.label1
    label2=cmd_args.label2
    
    plot1.plot(times1,a_sol1,color='k',label=label1,linewidth=linewidth)
    plot1.plot(times2,a_sol2,color='k',linestyle='dashed',label=label2,linewidth=linewidth2)
    plot1.plot(times1,a*(M_d*M_a/(M_d_sol1*M_a_sol1))**2,color='C3',linestyle='dotted',label="$\mathrm{Analytic\,(circular \,only)}$",linewidth=linewidth)

    plotB.plot(times1,M_d_sol1*M_a_sol1*np.sqrt(a_sol1*(1.0-e_sol1**2)),color='k',label=label1,linewidth=linewidth)
    plotB.plot(times2,M_d_sol2*M_a_sol2*np.sqrt(a_sol2*(1.0-e_sol2**2)),color='k',linestyle='dashed',label=label2,linewidth=linewidth2)


    plot2.plot(times1,e_sol1,color='k',linewidth=linewidth,label=label1)
    plot2.plot(times2,e_sol2,color='k',linestyle='dashed',linewidth=linewidth2,label=label2)
    plot3.plot(times1,M_d_sol1/M_a_sol1,color='k',linewidth=linewidth)
    plot3.plot(times2,M_d_sol2/M_a_sol2,color='k',linestyle='dashed',linewidth=linewidth2)
    plot3.plot(times1,E_0_sol1,color='C3',linestyle='dashed',label=label1)
    plot3.plot(times1,E_0_sol2,color='C3',linestyle='dashed',linewidth=linewidth2,label=label2)

    if plot_omega==True:
        plot3.plot(times2,omega_sol2,color='C0',linestyle='dashed',linewidth=linewidth2,label=r"$\omega\,(\tau\neq0)$")
    
    #pyplot.rc('text.latex', preamble='\usepackage{xcolor}')
    plots = [plot1,plot2,plot3]
    labels = [r"$a/\mathrm{AU}$",r"$e$",r"$q; \,\mathcal{E}_0/\mathrm{rad}$"]
    if plot_omega==True:
        labels = [r"$a/\mathrm{AU}$",r"$e$",r"$q; \,\mathcal{E}_0/\mathrm{rad}; \, \omega/\mathrm{rad}$"]
    for index,plot in enumerate(plots):
        if index==2:
            plot.set_xlabel(r"$t/\mathrm{Myr}$",fontsize=fontsize)
        plot.set_ylabel(labels[index],fontsize=fontsize)
                
        plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

        if index in [0,1]:
            plot.set_xticklabels([])

    if ylim3 != None:
        plot3.set_ylim(0,ylim3)
        
    handles,labels = plot1.get_legend_handles_labels()
    plot1.legend(handles,labels,loc=cmd_args.plot1_loc,fontsize=0.8*fontsize)

    handles,labels = plot2.get_legend_handles_labels()
    plot2.legend(handles,labels,loc=cmd_args.plot2_loc,fontsize=0.8*fontsize)

    handles,labels = plot3.get_legend_handles_labels()
    plot3.legend(handles,labels,loc=cmd_args.plot3_loc,fontsize=0.8*fontsize)

    plot3.axhline(y=1.0,linestyle='dotted',color='k')
    plot3.axhline(y=np.pi,linestyle='dotted',color='C3')
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    filename = cmd_args.filename
    fig.savefig(filename,dpi=200)
    
    pyplot.show()
    
def plot_function_triple_slow(cmd_args,data1,data2):
    if cmd_args.plot_fancy == True:
        pyplot.rc('text',usetex=True)
        pyplot.rc('legend',fancybox=True)
        #pyplot.rc('font', weight='bold')
       # pyplot.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
        
    a = args.a
    e = args.e
    M_d = args.M_d
    M_a = args.M_a
    times1,a_sol1,e_sol1,omega_sol1,M_d_sol1,M_a_sol1,E_0_sol1 = data1
    times2,a_sol2,e_sol2,omega_sol2,M_d_sol2,M_a_sol2,E_0_sol2 = data2
    times1 *= 1.0e-6
    times2 *= 1.0e-6


    fontsize=22
    labelsize=18
    linewidth=2.0
    linewidth2=1.0
    
    fig=pyplot.figure(figsize=(8,10))
    plot1=fig.add_subplot(3,1,1,yscale="log")
    plot2=fig.add_subplot(3,1,2)
    plot3=fig.add_subplot(3,1,3)
    
    label1=cmd_args.label1
    label2=cmd_args.label2
    
    plot1.plot(times2,a_sol2,color='k',linestyle='dashed',label=label2,linewidth=linewidth2)
    plot1.plot(times1,a_sol1,color='k',label=label1,linewidth=linewidth)
    
    plot2.plot(times2,e_sol2,color='k',linestyle='dashed',linewidth=linewidth2,label=label2)
    plot2.plot(times1,e_sol1,color='k',linewidth=linewidth,label=label1)
    
    plot3.plot(times2,M_d_sol2/M_a_sol2,color='k',linestyle='dashed',linewidth=linewidth2)
    plot3.plot(times1,M_d_sol1/M_a_sol1,color='k',linewidth=linewidth)
    
    plot3.plot(times1,E_0_sol2,color='C3',linestyle='dashed',linewidth=linewidth2,label=label2)
    plot3.plot(times1,E_0_sol1,color='C3',linestyle='dashed',label=label1,linewidth=linewidth)
    
    
    #pyplot.rc('text.latex', preamble='\usepackage{xcolor}')
    plots = [plot1,plot2,plot3]
    labels = [r"$a/\mathrm{AU}$",r"$e$",r"$q; \,\mathcal{E}_0/\mathrm{rad}$"]
    for index,plot in enumerate(plots):
        if index==2:
            plot.set_xlabel(r"$t/\mathrm{Myr}$",fontsize=fontsize)
        plot.set_ylabel(labels[index],fontsize=fontsize)
                
        plot.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True)

        if index in [0,1]:
            plot.set_xticklabels([])

    tmax = 17.0
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    plot1_inset = inset_axes(plot1,
    width="25%", # width = 30% of parent_bbox
    height=1.4, # height : 1 inch
    loc="center right")

    plot1_inset.plot(times1,a_sol1,color='k',linewidth=linewidth)
    plot1_inset.set_xlim(10,tmax)
    plot1_inset.set_ylim(0.5,1.0)
    plot1_inset.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True) #labelcolor='C3'


    plot2_inset = inset_axes(plot2,
    width="25%", # width = 30% of parent_bbox
    height=1.4, # height : 1 inch
    loc="center right")

    plot2_inset.plot(times1,e_sol1,color='k',linewidth=linewidth)
    plot2_inset.set_xlim(10,tmax)
    plot2_inset.set_ylim(0.5,1.0)
    plot2_inset.tick_params(axis='both', which ='major', labelsize = 1.0*labelsize,bottom=True, top=True, left=True, right=True,labelcolor='k')


    plot3_inset = inset_axes(plot3,
    width="25%", # width = 30% of parent_bbox
    height=1.4, # height : 1 inch
    loc="center right")

    plot3_inset.plot(times1,E_0_sol1,color='C3',linewidth=linewidth)
    plot3_inset.set_xlim(10,tmax)
    plot3_inset.set_ylim(0.0,0.15)
    plot3_inset.tick_params(axis='both', which ='major', labelsize = labelsize,bottom=True, top=True, left=True, right=True) #labelcolor='C3'
    
    handles,labels = plot1.get_legend_handles_labels()
    plot1.legend(handles,labels,loc=cmd_args.plot1_loc,fontsize=0.8*fontsize)

    handles,labels = plot2.get_legend_handles_labels()
    plot2.legend(handles,labels,loc=cmd_args.plot2_loc,fontsize=0.8*fontsize)

    handles,labels = plot3.get_legend_handles_labels()
    plot3.legend(handles,labels,loc=cmd_args.plot3_loc,fontsize=0.8*fontsize)

    plot3.axhline(y=1.0,linestyle='dotted',color='k')
    plot3.axhline(y=np.pi,linestyle='dotted',color='C3')
    
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    filename = cmd_args.filename
    fig.savefig(filename,dpi=200)
    
    pyplot.show()
    
class args:
    pass

if __name__ == '__main__':
    cmd_args = parse_arguments()
    application_id = cmd_args.application_id

    args = args()
    args.plot_fancy = cmd_args.plot_fancy
    args.verbose = cmd_args.verbose
    args.R_from_M_d = False
    args.R_M_d_power = 0.0
    args.scaled_t_end = True
    args.N_steps = 1000
    args.rp_min = 0.01
    args.include_ejection_radius = False
    args.include_accretion_radius = False
    args.ejection_radius_mode = 1
    args.donor_normalized_spin = 1.0
    args.accretion_radius = 1.0
    args.use_constant_tau = False
    args.tau = 0.0
    args.tau_mean_anomaly = 0.0
    args.include_tertiary = False
    args.include_quadrupole_terms = False
    args.include_octupole_terms = False
    args.include_1PN_terms = False
    args.include_inner_tidal_terms = False
    args.omega = 0.1
    args.omega_out = 0.1
    args.e_out = 0.1
    args.a_out = 100.0
    args.i_rel = 0.0
    args.mxstep = 10000
    
    args.spin_angular_frequency_d = 1.0
    args.spin_angular_frequency_a = 1.0
    args.R_a = 1.0
    args.k_div_T_tides_d = 0.0
    args.k_div_T_tides_a = 0.0
    args.gyration_radius_d = 0.28
    args.gyration_radius_a = 0.28
    args.apsidal_motion_constant_d = 0.014
    args.apsidal_motion_constant_a = 0.014
    
    
    if application_id==-1: ### test triple implementation with Naoz reference system ###
        args.e = 0.001
        args.e_out = 0.6
        args.a = 6.0
        args.a_out = 100.0
        args.M_d = 1.0
        args.M_a = 1.0*MJ_in_MSun
        args.M_t = 40.0*MJ_in_MSun
        args.R = 1.0
        args.i_rel = 65.0*np.pi/180.0
        args.rp_min = 1.0e-10
        args.omega = 45.0*np.pi/180.0
        args.include_tertiary = True
        args.include_quadrupole_terms = True
        args.include_octupole_terms = True
        args.include_1PN_terms = True

        args.scaled_t_end = False
        args.t_end = 30.0e6
        args.M_d_dot_av = -1.0e-20

        args.model = "emt"
        args.name = "test_Naoz"
        
        data = integrator.integrate(args)
        integrator.plot_function(args,data)
    if application_id==-2: ### test equilibrium tides implementation
        args.e = 0.9
        args.e_out = 0.6
        args.a = 1.0
        args.a_out = 100.0
        args.M_d = 1.0
        args.M_a = 1.0
        args.M_t = 0.0
        args.R = 10.0
        args.R_a = 5.0
        args.i_rel = 65.0*np.pi/180.0
        args.rp_min = 1.0e-10
        args.omega = 45.0*np.pi/180.0
        args.include_tertiary = True ### should be True!
        args.include_quadrupole_terms = False
        args.include_octupole_terms = False
        args.include_1PN_terms = False
        args.include_inner_tidal_terms = True

        args.k_div_T_tides_d = 1.0e-6
        args.k_div_T_tides_a = 1.0e-6

        args.scaled_t_end = False
        args.t_end = 1.0e10
        args.M_d_dot_av = 0.0

        args.model = "emt"
        args.name = "test_tides"
        
        data = integrator.integrate(args,return_spins=True)
        integrator.plot_function_tides(args,data)
    if application_id==0: ### circular binary
        args.M_d = 1.0
        args.M_d_dot_av = -1.0e-8

        args.M_a = 0.8
        args.e = 0.0
        args.a = 1.0
        args.R = 1000.0
        args.t_end = 0.7

        args.model = "emt"
        data_emt = integrator.integrate(args)

        args.model = "sep"

        data_sep = integrator.integrate(args)

        cmd_args.filename = "binary_circular.pdf"
        cmd_args.plot1_loc = cmd_args.plot2_loc = "upper left"
        cmd_args.plot2_loc = "center left"
        cmd_args.plot3_loc = "center right"        

        plot_function_compare_sep(cmd_args,data_emt,data_sep)

    elif application_id==10: ### circular binary; initial mass ratio <1
        args.M_d = 1.5
        args.M_d_dot_av = -1.0e-8

        args.M_a = 2.0
        args.e = 0.0
        args.a = 1.0
        args.R = 1000.0
        args.t_end = 0.7

        args.model = "emt"
        data_emt = integrator.integrate(args)

        args.model = "sep"

        data_sep = integrator.integrate(args)

        cmd_args.filename = "binary_circular_low_q.pdf"
        cmd_args.plot1_loc = cmd_args.plot2_loc = "upper left"
        cmd_args.plot2_loc = "center left"
        cmd_args.plot3_loc = "center right"        

        plot_function_compare_sep(cmd_args,data_emt,data_sep)

    elif application_id==1: ### eccentric binary -- emt vs sep
        args.M_d = 8.0
        args.M_d_dot_av = -1.0e-8

        args.rp_min = 10.0*RSun_in_AU
        args.M_a = 1.4
        
        args.a = 1.0
        args.R = 10.0
        args.t_end = 0.5
        
        args.e = 0.92
        args.model = "emt"
        data_emt = integrator.integrate(args)

        args.model = "sep"
        data_sep = integrator.integrate(args)

        args.e = 0.45
        args.a = 0.1
        args.model = "emt"
        extra_A_data1 = integrator.integrate(args)

        args.model = "sep"
        extra_A_data2 = integrator.integrate(args)

        args.e = 0.2
        args.a = 0.08
        args.model = "emt"
        extra_B_data1 = integrator.integrate(args)

        args.model = "sep"
        extra_B_data2 = integrator.integrate(args)

        ### change args.a back to 1.0 for plotting purposes ###
        args.a = 1.0

        cmd_args.filename = "binary_eccentric.pdf"
        cmd_args.plot1_loc = "upper right"
        cmd_args.plot2_loc = "upper right"
        cmd_args.plot3_loc = "upper right"        
        plot_function_compare_sep(cmd_args,data_emt,data_sep,ylim3=2.0*np.pi,\
            extra_A_data1=extra_A_data1,extra_A_data2=extra_A_data2, \
            extra_B_data1=extra_B_data1,extra_B_data2=extra_B_data2, \
            label_A_1 = r"$\mathrm{emt}\,(a=0.1\,\mathrm{AU}, \, e=0.45)$", \
            label_A_2 = r"$\mathrm{Sepinsky}\,(a=0.1\,\mathrm{AU}, \, e=0.45)$", \
            label_B_1 = r"$\mathrm{emt}\,(a=0.08\,\mathrm{AU}, \, e=0.2$)", \
            label_B_2 = r"$\mathrm{Sepinsky}\,(a=0.08\,\mathrm{AU}, \, e=0.2$)")

    elif application_id==2: ### eccentric binary -- zero ejection/accretion radii vs nonzero
        args.M_d = 8.0
        args.M_d_dot_av = -1.0e-8

        args.rp_min = 10.0*RSun_in_AU
        args.M_a = 1.4
        args.e = 0.92
        args.a = 1.0
        args.R = 10.0
        args.t_end = 0.5
        
        args.model = "emt"
        data_zero = integrator.integrate(args)

        args.include_ejection_radius = True
        args.include_accretion_radius = True

        args.ejection_radius_mode = 1
        args.donor_normalized_spin = 1.0
        args.accretion_radius = 0.01


        data_nonzero = integrator.integrate(args)

        cmd_args.plot1_loc = "upper right"
        cmd_args.plot2_loc = "upper right"        
        cmd_args.plot3_loc = "upper right"        

        cmd_args.label1 = r"$\mathbf{r}_\mathrm{A} = \mathbf{0}$"
        cmd_args.label2 = r"$\mathbf{r}_\mathrm{A} \neq \mathbf{0}$"
        cmd_args.filename = "binary_eccentric_nonzero.pdf"
        plot_function_compare_size(cmd_args,data_zero,data_nonzero,ylim3=2.0*np.pi)

    elif application_id==3: ### eccentric binary -- ejection mode 1 vs mode 2
        args.M_d = 8.0
        args.M_d_dot_av = -1.0e-8

        args.rp_min = 10.0*RSun_in_AU
        args.M_a = 1.4
        args.e = 0.92
        args.a = 1.0
        args.R = 10.0
        args.t_end = 0.7
        args.mxstep=10000
        
        args.include_ejection_radius = True
        args.include_accretion_radius = True

        
        args.donor_normalized_spin = 1.0
        args.accretion_radius = 0.01

        args.model = "emt"
        
        args.ejection_radius_mode = 1
        data1 = integrator.integrate(args)


        args.ejection_radius_mode = 2
        data2 = integrator.integrate(args)
        
        cmd_args.plot1_loc = "upper right"
        cmd_args.plot2_loc = "upper center"        
        cmd_args.plot3_loc = "upper right"        

        cmd_args.label1 = r"$\mathrm{Low\,spin}$"
        cmd_args.label2 = r"$\mathrm{High\,}q$"
        cmd_args.filename = "binary_eccentric_nonzero_mode.pdf"
        plot_function_compare_size(cmd_args,data1,data2,ylim3=2.0*np.pi)

    elif application_id==4: ### eccentric binary -- zero vs nonzero tau
        args.M_d = 8.0
        args.M_d_dot_av = -1.0e-8

        args.rp_min = 10.0*RSun_in_AU
        args.M_a = 1.4
        args.e = 0.92
        args.a = 1.0
        args.R = 10.0
        args.t_end = 0.5
        args.mxstep=40000

        args.model = "emt"
        data1 = integrator.integrate(args)

        P_orb = 2.0*np.pi*np.sqrt(args.a**3/(G*(args.M_d+args.M_a)))
        t_hyd = np.sqrt( ((args.R*RSun_in_AU)**3)/(G*args.M_d**3))
        tau = 10*t_hyd
        print 't_hyd/yr',t_hyd,'t_hyd/d',t_hyd*365.25,'tau/yr',tau

        ### Modify lines below to determine whether tau should be constant, or n*tau should be constant
        args.tau = tau
        args.use_constant_tau = True
        #args.use_constant_tau = False
        #args.tau_mean_anomaly = 0.1
        data2 = integrator.integrate(args)

        cmd_args.plot1_loc = "lower right"
        cmd_args.plot2_loc = "center right"        
        cmd_args.plot3_loc = "upper right"        

        cmd_args.label1 = r"$\tau=0$"
        cmd_args.label2 = r"$\tau\neq0$"
        cmd_args.filename = "binary_eccentric_tau.pdf"
        plot_function_compare_size(cmd_args,data1,data2,plot_omega=True)

    if application_id==5: ### triple (fast mass transfer)

        args.R = 1.0
        args.rp_min = 1.0*RSun_in_AU
        args.mxstep=40000

        args.e = 0.001
        args.e_out = 0.6
        args.a = 1.0
        args.a_out = 150.0
        args.M_d = 1.0
        args.M_a = 0.1
        args.M_t = 1.0
        args.i_rel = 85.0*np.pi/180.0
        
        args.omega = 45.0*np.pi/180.0
        args.include_tertiary = True
        args.include_quadrupole_terms = True
        args.include_octupole_terms = True
        args.include_1PN_terms = True

        args.scaled_t_end = False
        args.t_end = 2.0e7
        args.M_d_dot_av = -1.0e-20

        args.model = "emt"
        
        import pickle
        if cmd_args.calc == True:
            args.M_d_dot_av = -1.0e-30
            data1 = integrator.integrate(args)
            
            args.M_d_dot_av = -1.0e-8
            data2 = integrator.integrate(args)

            data = data1,data2
            pickle.dump( data, open( "applications_data_id_5.pickle", "wb" ) )
    
        try:
            data = pickle.load( open( "applications_data_id_5.pickle", "rb" ) )
            data1,data2 = data
        except:
            print 'Error loading pickled data -- did you first run the script with --calc?'
            exit(-1)


        cmd_args.plot1_loc = "center left"
        cmd_args.plot2_loc = "center left"        
        cmd_args.plot3_loc = "center left"        
        cmd_args.label1 = r"$\mathrm{No \,mass\,transfer}$"
        cmd_args.label2 = r"$\mathrm{Mass\,transfer}$"
        cmd_args.filename = "triple_fast.pdf"

        plot_function_compare_size(cmd_args,data1,data2)


    if application_id==6: ### triple (slow mass transfer)

        args.R = 1.0
        args.rp_min = 1.0*RSun_in_AU
        args.mxstep=40000

        args.e = 0.001
        args.e_out = 0.6
        args.a = 1.0
        args.a_out = 150.0
        args.M_d = 1.0
        args.M_a = 0.1
        args.M_t = 1.0
        args.i_rel = 85.0*np.pi/180.0
        
        args.omega = 45.0*np.pi/180.0
        args.include_tertiary = True
        args.include_quadrupole_terms = True
        args.include_octupole_terms = True
        args.include_1PN_terms = True

        args.scaled_t_end = False
        args.N_steps = 10000
        args.t_end = 1.0e8
        args.M_d_dot_av = -1.0e-20

        args.model = "emt"
        
        import pickle
        if cmd_args.calc == True:
            args.M_d_dot_av = -1.0e-9
            data1 = integrator.integrate(args)
        
            args.M_d_dot_av = -1.0e-30
            data2 = integrator.integrate(args)
            
            
            data = data1,data2
            pickle.dump( data, open( "applications_data_id_6.pickle", "wb" ) )
    
        try:
            data = pickle.load( open( "applications_data_id_6.pickle", "rb" ) )
            data1,data2 = data
        except:
            print 'Error loading pickled data -- did you first run the script with --calc?'
            exit(-1)

        cmd_args.plot1_loc = "center left"
        cmd_args.plot2_loc = "center left"        
        cmd_args.plot3_loc = "center left"        
        cmd_args.label2 = r"$\mathrm{No \,mass\,transfer}$"
        cmd_args.label1 = r"$\mathrm{Mass\,transfer}$"
        cmd_args.filename = "triple_slow.pdf"

        plot_function_triple_slow(cmd_args,data1,data2)

import ctypes

import os
path = os.path.dirname(os.path.realpath(__file__))
emtlibrary = ctypes.CDLL(path+'/emtlibrary.so')

functions = [emtlibrary.fm,emtlibrary.fa,emtlibrary.fe,emtlibrary.fomega]

for i,f in enumerate(functions):

    f.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)
    f.restype = ctypes.c_double

functions = [emtlibrary.ga,emtlibrary.ge,emtlibrary.ha,emtlibrary.he]

for i,f in enumerate(functions):

    f.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double)
    f.restype = ctypes.c_double


f = emtlibrary.R_Lc_div_a
f.argtypes = (ctypes.c_double,)
f.restype = ctypes.c_double

f = emtlibrary.compute_eccentric_anomaly_from_mean_anomaly
f.argtypes = (ctypes.c_double,ctypes.c_double)
f.restype = ctypes.c_double

f = emtlibrary.XL0_q
f.argtypes = (ctypes.c_double,)
f.restype = ctypes.c_double

def temp_(x):
    y = ctypes.c_double(0.0)
    z = ctypes.c_double(0.0)
    emtlibrary.temp_.argtypes = (ctypes.c_double,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
    emtlibrary.temp_.restype = ctypes.c_void_p
    emtlibrary.temp_(x, ctypes.byref(y), ctypes.byref(z))
    
    return (y.value, z.value)

emtlibrary.temp = temp_


def triple_EOM_(CONST_G,CONST_C,e_in,e_out,g_in,g_out,cositot, \
    a_in,a_out,m1,m2,m3, \
    include_quadrupole_terms,include_octupole_terms,include_1PN_terms):

    de_in_dt = ctypes.c_double(0.0)
    de_out_dt = ctypes.c_double(0.0)
    dg_in_dt = ctypes.c_double(0.0)
    dg_out_dt = ctypes.c_double(0.0)
    dcositot_dt = ctypes.c_double(0.0)
    
    emtlibrary.triple_EOM_.argtypes = (ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,\
        ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,\
        ctypes.c_int,ctypes.c_int,ctypes.c_int,\
        ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),\
        ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
    
    emtlibrary.triple_EOM_.restype = ctypes.c_void_p
    emtlibrary.triple_EOM_(CONST_G,CONST_C,e_in,e_out,g_in,g_out,cositot, \
        a_in,a_out,m1,m2,m3, \
        include_quadrupole_terms,include_octupole_terms,include_1PN_terms,\
        ctypes.byref(de_in_dt), ctypes.byref(de_out_dt),\
        ctypes.byref(dg_in_dt), ctypes.byref(dg_out_dt),
        ctypes.byref(dcositot_dt))
    
    return (de_in_dt.value, de_out_dt.value, dg_in_dt.value, dg_out_dt.value, dcositot_dt.value)

emtlibrary.triple_EOM = triple_EOM_

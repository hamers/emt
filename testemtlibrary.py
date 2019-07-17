from wrapperemtlibrary import emtlibrary

def test_f_functions():
    assert round(emtlibrary.fm(0.8,1.1,0.8,0.1),7) == 0.0223401
    assert round(emtlibrary.fa(0.8,1.1,0.8,0.1),6) == 0.14934
    assert round(emtlibrary.fe(0.8,1.1,0.8,0.1),7) == 0.0285751
    assert round(emtlibrary.fomega(0.8,1.1,0.8,0.1),7) == 0.0036584
    
def test_g_functions():
    assert round(emtlibrary.ga(0.8,1.1,0.8),7) == 0.0866617
    assert round(emtlibrary.ge(0.8,1.1,0.8),7) == 0.0144378

def test_h_functions():
    assert round(emtlibrary.ha(0.8,1.1,0.8),6) == 0.406088
    assert round(emtlibrary.he(0.8,1.1,0.8),6) == 0.071878

def test_R_L_c():
    assert round(emtlibrary.R_Lc_div_a(1.0),11) == 0.37892051838

def test_compute_eccentric_anomaly_from_mean_anomaly():
    assert round(emtlibrary.compute_eccentric_anomaly_from_mean_anomaly(0.1,0.8),12) == 0.442716567419

def test_XL0_q():
    assert round(emtlibrary.XL0_q(0.1),12) == 0.305361080053


def test_triple_EOM():
    return_triple = emtlibrary.triple_EOM(1.0,1.0,0.1,0.1,0.1,0.1,0.1,1.0,100.0,1.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1,1,1,1)
    
    assert round(return_triple[1],15) == 2.6340081e-08
    assert round(return_triple[2],15) == 0.0
    assert round(return_triple[3],15) == 8.57099233534009
    assert round(return_triple[4],15) == 3.8817792e-08

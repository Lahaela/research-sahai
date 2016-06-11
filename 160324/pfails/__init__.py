from __future__ import division
from matplotlib import *
from numpy import *
import scipy as sp
from scipy import stats
from decimal import *
import pickle
from scipy.interpolate import interp1d
import operator as op

def nCr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

# For scipy Linear Interpolation
def hs_rs_table(rate, blocklength):
    k = (1-rate)*blocklength
    msg = rate*blocklength
    op_SNR = arange(0, 5, 0.01)
    pbitdrop = array([0.5*math.erfc(sqrt(10**opsnr/2)) for opsnr in op_SNR])
    hcerr = 1 - ((1-pbitdrop)**7 + 7*pbitdrop*(1-pbitdrop)**6)
    hcflst = 1 - (1-hcerr)**2
    reeddrop = [sum([nCr(blocklength, d)* hcf**d *(1-hcf)**(blocklength-d) for d in range(int(k/2)+1, blocklength)]) \
                for hcf in hcflst]
#     func = interp1d(op_SNR, reeddrop, kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    return [op_SNR, reeddrop]

def save_table(table, filename):
    f = open(filename, 'w')
    pickle.dump(table, f)
    f.close()
    
def load_table(filename):
    f = open(filename, 'r')
    table = pickle.load(f)
    f.close()
    return table

# Linear interpolation function
def p_single(codetable, op_SNR, endpoint, dfade, tabledSNR):
    fadexp = sp.stats.expon()
    psingle = []
    fade = arange(0, endpoint[0], dfade[0])
    snrlookup = codetable(op_SNR+log10(fade))
    fadepr = fadexp.pdf(fade)
    psingle.append(dfade[0]*dot(snrlookup, fadepr))
    for idx in range(1, len(endpoint)):
        fade = arange(endpoint[idx-1], endpoint[idx], dfade[idx])
        snrlookup = codetable(op_SNR+log10(fade))
        fadepr = fadexp.pdf(fade)
        psingle.append(dfade[idx]*dot(snrlookup, fadepr))
    return sum(psingle)

def p_combo(codetable, a, op_SNR, endpoint, dfade, tabledSNR):
    if a == 0:
        return 1.0
    fadexp = sp.stats.erlang(a)
    pcombo = []
    fade = arange(0, endpoint[0], dfade[0])
    snrlookup = codetable(op_SNR+log10(fade))
    fadepr = fadexp.pdf(fade)
    pcombo.append(dfade[0]*dot(snrlookup, fadepr))
    for idx in range(1, len(endpoint)):
        fade = arange(endpoint[idx-1], endpoint[idx], dfade[idx])
        snrlookup = codetable(op_SNR+log10(fade))
        fadepr = fadexp.pdf(fade)
        pcombo.append(dfade[idx]*dot(snrlookup, fadepr))
    return sum(pcombo)

# N = num_nodes
def p_protocol(codetable, N, op_SNR, endpoint, dfade, tabledSNR):
    psingle = p_single(codetable, op_SNR, endpoint, dfade, tabledSNR)
    return sum([nCr(N, a) * (1-psingle)**a * psingle**(N-a) * 
        (1-(1-p_combo(codetable, a, op_SNR, endpoint, dfade, tabledSNR))**(N-a)) for a in range(N)])

def shannon_combo(N, a, rate1, p1, SNR):
    rate2 = (N-a)/N*rate1 + 2*N/10000
    p2 = 1 - exp(-(2**rate2 -1)/(10**(SNR/10)))
    return p2**a * min(p2/p1, 1)

from __future__ import division
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import argrelextrema
from scipy.interpolate import interp1d
import operator as op
from decimal import *
import pickle
import math

def nCr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

def Q(x):
    return 0.5*sp.special.erfc(x/np.sqrt(2))

def Qinv(x):
    return (np.sqrt(2)*sp.special.erfinv(1-2*x))**2

# Table Generation
# rate = num_nodes * 160 / 10,000
# blocklength is for RS (not bits)
def hs_rs_table(op_SNR, rate, blocklength):
    rate = rate * 7/4
    k = (1-rate)*blocklength
    # msg = rate*blocklength
    pbitdrop = Q(np.sqrt(2*10**(op_SNR/10)))
    hcerr = 1 - ((1-pbitdrop)**7 + 7*pbitdrop*(1-pbitdrop)**6)
    hcflst = 1 - (1-hcerr)**3
    reeddrop = [sum([nCr(blocklength, d)* hcf**d *(1-hcf)**(blocklength-d) for d in range(int(k/2)+1, blocklength)]) \
                for hcf in hcflst]
    return [op_SNR, np.array(reeddrop)]

def shannon_table(rate):
    op_SNR = np.arange(-20, 5, 0.01)
    shannon = []
    for SNR in op_SNR:
        C = np.log2(1 + 10**(SNR/10))
        if rate > C:
            shannon.append(1.0) # pr failure
        else:
            shannon.append(0.0)
    return [op_SNR, np.array(shannon)]

def hard_bin_shannon_table(rate):
    op_SNR = np.arange(-20, 5, 0.01)
    shannon = []
    SNR = 2*10**(op_SNR/10)
    pbit = Q(np.sqrt(SNR))
    for p in pbit:
        C = 1 + p*np.log2(p) + (1-p)*np.log2(1-p)
        if rate > C:
            shannon.append(1.0) # pr failure
        else:
            shannon.append(0.0)
    return [op_SNR, np.array(shannon)]

# Doesn't make sense for error_prob > 0.5 -- ignore that part of waterfall plot
def polyanskiy_table(rate, blocklength=10000):
    op_SNR = np.arange(-25, 5, 0.01)
    pdrop = []
    for SNR in op_SNR:
        C = np.log2(1 + 10**(SNR/10))
        V = (np.log2(math.e))**2 * (1 - (1 / (1 + 10**(SNR/10)))**2)
        error_prob = Q(sqrt(blocklength/V) * (C - rate))
        pdrop.append(error_prob)
#         error_prob = Q(sqrt(2*10**(SNR/10)))
#         Cprime = C - np.sqrt(V/blocklength) * sqrt(2*10**(SNR/10)) # Qinv(error_prob)
    return [op_SNR, np.array(pdrop)]

def save_table(table, filename):
    f = open(filename, 'w')
    pickle.dump(table, f)
    f.close()
    
def load_table(filename):
    f = open(filename, 'r')
    table = pickle.load(f)
    f.close()
    return table

# Occupy CoW
def shannon_combo(N, a, rate1, p1, SNR):
#     rate2 = (N-a)/N*rate1 + 2*N/10000 # Rate adaptation
    rate2 = rate1
    p2 = 1 - np.exp(-(2**rate2 -1)/(10**(SNR/10)))
    return p2**a * min(p2/p1, 1)

def adaptive_endpoints(fx, fade):
    grad = np.gradient(fx, 10**(-3))*10**(-3) # Magic
#     grad2 = np.gradient(grad, 10**(-4))*10**(-4)
    
    ind = argrelextrema(grad, np.less)[0]
    ind = ind[np.argsort(grad[ind])]
    if len(ind) > 0:
        mid = fade[ind[0]]
        end1, end2 = mid - 5000*10**(-9), mid + 5000*10**(-9) # Magic
    else:
        mid = 0
        end1, end2 = fade[0], fade[0]+10**(-3)
    
    endpts = [end1, end2, 2]
    return endpts

# Adaptive, Energy Combining
def p_single(codetable, op_SNR, endpoint, dfade):
    fadexp = sp.stats.expon()
    fade = arange(0, endpoint, 10**(-3)) # Magic
    psingle = []
#     snrlookup = codetable(op_SNR+log10(fade))
#     fadepr = fadexp.pdf(fade)
    optimize = codetable(op_SNR+10*np.log10(fade)) * fadexp.pdf(fade)
    endpts = adaptive_endpoints(optimize, fade)
    fade = arange(0, endpts[0], dfade[0])
    fadepr = fadexp.pdf(fade)
    snrlookup = codetable(op_SNR+10*np.log10(fade))
    psingle.append(dfade[0]*np.dot(snrlookup, fadepr))
    for idx in range(1, len(endpts)):
        fade = arange(endpts[idx-1], endpts[idx], dfade[idx])
        snrlookup = codetable(op_SNR+10*np.log10(fade))
        fadepr = fadexp.pdf(fade)
        psingle.append(dfade[idx]*np.dot(snrlookup, fadepr))
    return sum(psingle)

def p_combo(codetable, a, op_SNR, endpoint, dfade):
    if a == 0:
        return 1.0
    fadexp = sp.stats.erlang(a)
    pcombo = []
    fade = arange(0, endpoint, 10**(-3))
    snrlookup = codetable(op_SNR+10*np.log10(fade))
    fadepr = fadexp.pdf(fade)
    optimize = snrlookup * fadepr
    endpts = adaptive_endpoints(optimize, fade)
    fade = arange(0, endpts[0], dfade[0])
    snrlookup = codetable(op_SNR+10*np.log10(fade))
    fadepr = fadexp.pdf(fade)
    pcombo.append(dfade[0]*np.dot(snrlookup, fadepr))
    for idx in range(1, len(endpts)):
        fade = arange(endpts[idx-1], endpts[idx], dfade[idx])
        snrlookup = codetable(op_SNR+10*np.log10(fade))
        fadepr = fadexp.pdf(fade)
        pcombo.append(dfade[idx]*np.dot(snrlookup, fadepr))
    return sum(pcombo)

def p_protocol(codetable, N, op_SNR, endpoint, dfade):
    psingle = p_single(codetable, op_SNR, endpoint, dfade)
    return sum([nCr(N, a) * (1-psingle)**a * psingle**(N-a) * 
                (1-(1-p_combo(codetable, a, op_SNR, endpoint, dfade))**(N-a)) for a in range(N)])

def energy_combining(codingscheme, dSNR, dfade, endpoint, threshold, start_SNR, start_nodes, end_nodes):
    nominal_needed = []
    for N in range(start_nodes, end_nodes):
        filename = codingscheme + '/n' + str(N) + '.in'
        codetable = load_table(filename)
        func = interp1d(codetable[0], codetable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
        SNR = start_SNR-dSNR # dB scale
        pprotocol = 1.0
        while pprotocol > threshold:
            SNR += dSNR
            pprotocol = p_protocol(func, N, SNR, end, dfade)
        print('P(protocol)', N, SNR, pprotocol)
        nominal_needed.append(SNR)

    plt.plot(range(start_nodes, end_nodes), np.array(nominal_needed), lw=2.0, label=codingscheme)
    plt.xlabel('Number of Nodes', fontsize=18)
    plt.ylabel('Nominal SNR Needed (dB)', fontsize=18, labelpad=20)
    plt.title('Hockey {0} ECC'.format(threshold), fontsize=24)
    
    return np.array(nominal_needed)

# Loudest Talker Functions
def padd(pcrit, k, rsblock, numwords):
    hcerr = 1 - ((1-pcrit)**7 + 7*pcrit*(1-pcrit)**6)
    hccrit = 1 - (1-hcerr)**numwords
    paddnoise = sum([nCr(rsblock, d)* hccrit**d *(1-hccrit)**(rsblock-d) for d in range(int(k/2)+1, rsblock)])

def loudest_talker(codingscheme, dSNR, target, paddratio, start_SNR, start_nodes, end_nodes):
    nomSNR = []

    for N in range(start_nodes, end_nodes):
        filename = codingscheme + '/n' + str(N) + '.in'
        codetable = load_table(filename)

        padd = paddratio*target
        actualSNR = codetable[0][np.where(np.array(codetable[1])<=padd)[0][0]]

        SNR = start_SNR
        pprotocol = 1.0
        while pprotocol > target:
            SNR += dSNR
            hcrit = 10**((actualSNR - SNR)/10) # linear fade
            pbadfade = 1 - np.exp(-hcrit)
            psingle = pbadfade + (1-pbadfade)*padd

            pprotocol = sum([nCr(N, a) * (1-psingle)**a * psingle**(N-a) * 
                             (1-(1-(pbadfade**a + (1-pbadfade**a) * padd))**(N-a)) for a in range(N)])
        nomSNR.append(SNR)
        # print('Loudest Speaker', N, SNR, actualSNR)

    # plt.plot(range(start_nodes, end_nodes), np.array(nomSNR), lw=2.0, label=codingscheme)
    # plt.xlabel('Number of Nodes', fontsize=18)
    # plt.ylabel('Nominal SNR Needed (dB)', fontsize=18, labelpad=20)
    # plt.title('Hockey 10^{-9} Loudest Talker', fontsize=24)
    
    return np.array(nomSNR)

# Uplink Functions
def uplink(codingscheme, dSNR, target, paddratio, start_SNR, start_nodes, end_nodes):
    nomSNR = []

    for N in range(start_nodes, end_nodes):
        filename = codingscheme + '/n' + str(N) + '.in'
        codetable = load_table(filename)

        padd = paddratio*target
        actualSNR = codetable[0][np.where(np.array(codetable[1])<=padd)[0][0]]

        SNR = start_SNR
        pprotocol = 1.0
        while pprotocol > target:
            SNR += dSNR
            hcrit = 10**((actualSNR - SNR)/10) # linear fade
            pbadfade = 1 - np.exp(-hcrit)
            psingle = pbadfade + (1-pbadfade)*padd

            psuccess = sum([nCr(N, a) * (1-psingle)**a * psingle**(N-a) *
                ((1-(pbadfade**a + (1-pbadfade**a)*padd)**a)*(1-padd))**(N-a) for a in range(1,N+1)])
            pprotocol = 1 - psuccess
        nomSNR.append(SNR)
    #     print('Protocol', N, SNR, actualSNR)

    # plt.plot(range(start_nodes, end_nodes), np.array(nomSNR), lw=2.0, label=codingscheme)
    # plt.xlabel('Number of Nodes', fontsize=18)
    # plt.ylabel('Nominal SNR Needed (dB)', fontsize=18, labelpad=20)
    # plt.title('Hockey 10^{-9} Uplink', fontsize=24)
    
    return np.array(nomSNR)

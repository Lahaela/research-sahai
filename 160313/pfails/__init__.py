from __future__ import division
from matplotlib import *
from numpy import *
import scipy as sp
from scipy import stats
import csv
from decimal import *
import operator as op


def nCr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

def load_table(filename):
    codetable = {}
    with open(filename, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if row[0] == 'SNR (dB)':
                continue
            codetable[Decimal(row[0])] = float(row[1])
    return codetable

def lookup(codetable, kmin, kmax, snr):
    if snr < kmin: return 1.0
    elif snr > kmax: return 0.0
    return codetable[snr]

def snr_lookup(codetable, kmin, kmax, opSNR, dSNR):
    op_SNR = Decimal(str(opSNR))
    if op_SNR < kmin: return 1.0
    elif op_SNR > kmax: return 0.0
    if op_SNR in codetable:
        return codetable[op_SNR]
    lsnr = op_SNR - op_SNR % Decimal(str(dSNR))
    rsnr = op_SNR + Decimal(str(dSNR)) - op_SNR % Decimal(str(dSNR))
#     lsnr = op_SNR.quantize(Decimal(str(dSNR)), rounding=ROUND_DOWN)
#     rsnr = op_SNR.quantize(Decimal(str(dSNR)), rounding=ROUND_UP)
    lp, rp = lookup(codetable, kmin, kmax, lsnr), lookup(codetable, kmin, kmax, rsnr)
    return (rp-lp)*(opSNR-float(lsnr))/dSNR + lp

# Changing bin width
def p_single(codetable, kmin, kmax, op_SNR, endpoint, dfade, tabledSNR):
    fadexp = sp.stats.expon()
    psingle = [dfade[0]]
    fade = arange(dfade[0], endpoint[0], dfade[0])
    psingle.append(dfade[0]*sum([snr_lookup(codetable, kmin, kmax, op_SNR+log10(f), tabledSNR) * fadexp.pdf(f) for f in fade]))
    for idx in range(1, len(endpoint)):
        fade = arange(endpoint[idx-1], endpoint[idx], dfade[idx])
        psingle.append(dfade[idx]*sum([snr_lookup(codetable, kmin, kmax, op_SNR+log10(f), tabledSNR) * fadexp.pdf(f) for f in fade]))
    return sum(psingle)

# Changing bins
def p_combo(codetable, kmin, kmax, a, op_SNR, endpoint, dfade, tabledSNR):
    if a == 0:
        return 1.0
    fadexp = sp.stats.erlang(a)
    pcombo = [dfade[0]]
    fade = arange(dfade[0], endpoint[0], dfade[0])
    pcombo.append(sum([snr_lookup(codetable, kmin, kmax, op_SNR+log10(f), tabledSNR) * fadexp.pdf(f) * dfade[0] for f in fade]))
    for idx in range(1, len(endpoint)):
        fade = arange(endpoint[idx-1], endpoint[idx], dfade[idx])
        pcombo.append(sum([snr_lookup(codetable, kmin, kmax, op_SNR+log10(f), tabledSNR) * fadexp.pdf(f) * dfade[idx] for f in fade]))
    return sum(pcombo)

# N = num_nodes
def p_protocol(codetable, kmin, kmax, N, op_SNR, endpoint, dfade, tabledSNR):
    psingle = p_single(codetable, kmin, kmax, op_SNR, endpoint, dfade, tabledSNR)
#     print('psingle', psingle, 'pcombo, ')
    return sum([nCr(N, a) * (1-psingle)**a * psingle**(N-a) * (1-(1-p_combo(codetable, kmin, kmax, a, op_SNR, endpoint, dfade, tabledSNR))**(N-a)) for a in range(N)])
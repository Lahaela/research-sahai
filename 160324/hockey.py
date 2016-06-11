from __future__ import division
from matplotlib.pyplot import *
from numpy import *
import scipy as sp
from scipy import stats
from decimal import *
import pickle
from scipy.interpolate import interp1d
import operator as op
from pfails import *

# Hockey Stick Plot
nominal_needed = []
shannon = []
dSNR, dfade, end = 0.01, [10**(-4), 10**(-2)], [1, 2]
threshold = 10**(-3)
for N in range(1, 35):
    filename = 'lookup_0-01/n' + str(N) + '.in'
    codetable = load_table(filename)
    func = interp1d(codetable[0], codetable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    SNR = 0 # dB scale
    pprotocol = 1.0
    while pprotocol > threshold:
        pprotocol = p_protocol(func, N, SNR, end, dfade, dSNR)
        SNR += 0.01
#         print(SNR, pprotocol)
    nominal_needed.append(SNR-0.01)
    
    rate = N*160/10000
    SNR = 0 # dB scale
    pshannon = 1.0
    while pshannon > threshold:
        plink = 1 - exp(-(2**rate -1)/(10**(SNR/10)))
        pshannon = sum([nCr(N, a) * (1-plink)**a * plink**(N-a) * 
                         (1-(1-shannon_combo(N, a, rate, plink, SNR))**(N-a)) for a in range(N)])
        SNR += 0.01
    shannon.append(SNR-0.01)


plot(range(1, 35), 10*array(nominal_needed), lw=2.0, label='4-3 Hamming')
plot(range(1, 35), array(shannon), lw=2.0, label='Shannon')
# legend(loc=0)
xlabel('Number of Nodes', fontsize=18)
ylabel('Nominal SNR Needed (dB)', fontsize=18, labelpad=20)
title('Hockey 10^{-3}', fontsize=24)
savefig('hockey_script_3.pdf', bbox='tight')

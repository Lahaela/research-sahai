from __future__ import division
from matplotlib.pyplot import *
from numpy import *
import scipy as sp
from scipy import stats
import csv
from decimal import *
import operator as op
from pfails import *

dSNR, dfade, end = 0.01, [0.01, 0.1], [1, 2]
N = 5
codetable = load_table('lookup_tables/n5.csv') # table resolution 0.01dB
SNRdB = arange(0, log10(10**5/50), 0.01)
# SNR_lst = 10**x
y = [p_protocol(codetable, 0.0, 5.0, N, SNR, end, dfade, dSNR) for SNR in SNRdB]
plot(10*SNRdB, y, lw=2.0)
title('P(Protocol fail)', fontsize=24)
xlabel('Operating SNR (dB)', fontsize=20)
savefig('testprotocol.pdf', bbox='tight')
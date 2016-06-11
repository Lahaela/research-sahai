from __future__ import division
from matplotlib.pyplot import *
from numpy import *
import scipy as sp
from scipy import stats
import csv
from decimal import *
import operator as op
from pfails import *

# Hockey Stick Plot

nominal_needed = []
for N in range(1, 35):
    filename = 'lookup_tables/n' + str(N) + '.csv'
    codetable = load_table(filename)
    tablemin, tablemax = min(codetable), max(codetable)
    dSNR, dfade, end = 0.01, [10**(-2), 10**(-1)], [1, 2]
    SNR = 0.005
    pprotocol = 1.0
    while pprotocol > 10**(-3):
        pprotocol = p_protocol(codetable, tablemin, tablemax, N, SNR, end, dfade, dSNR)
    nominal_needed.append(SNR)

plot(range(1, 35), nominal_needed, lw=2.0)
xlabel('Number of Nodes', fontsize=20)
ylabel('Nominal SNR Needed (dB)', fontsize=20)
savefig('hockey3.pdf', bbox='tight')
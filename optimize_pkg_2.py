from __future__ import division
import numpy as np
import scipy as sp
from scipy.stats import binom
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
from cow import *

protocol_target = 10**(-9)
extension = ".in"

# 21 from Hamming Code
# Quantize by 21 (to get all possible RS blocklengths)
# toy: 1400 = min bits for one phase, 4200 = entire protocol
# hs_rs_table(op_SNR, rate, blocklength)
class node_table():
    op_SNR = np.arange(-2.5, 11, 0.01)

    # num_nodes = n
    def __init__(self, n, protocol=40000, down=True):
        self.num_node = n
        minb = n * 160 * 7 / 4
#         quantize = (20000-minb*2)/100
        bitrange = np.arange(minb, protocol-minb, 21)
        self.bitrange = bitrange
        num_tables = len(bitrange)
        self.tables = {}
        for b in self.bitrange:
#             b = int(bitrange[bidx])
            block = b/21 if down else b/21/n
            self.tables[b] = hs_rs_table(self.op_SNR, n*160/b, int(block))

# 21 is arbitrary quantization taken from Hamming quantization
class node_table_shannon():
    op_SNR = np.arange(-16, -3, 0.01)

    def __init__(self, n):
        self.num_node = n
        minb = n * 160 * 7 / 4
        bitrange = np.arange(minb, 3*1400-minb, 21)
        self.bitrange = bitrange
        self.tables = {}
        for b in self.bitrange:
            self.tables[b] = shannon_table(n*160/b, self.op_SNR)

def leah_log10(x):
    return np.log10(x) if x else 0

def xor_analysis_new(N, tSNR, rSNRdu, rSNR3, p_a1=10**(-9), p_a2=10**(-9), p_a3=10**(-9)):
    if rSNRdu > tSNR: return 0
    h_du = 10**((rSNRdu - tSNR)/10) # linear fade
    h_xor = 10**((rSNR3 - tSNR)/10)
    # Probability fade is bad
    p_f1 = p_f2 = 1 - np.exp(-h_du)
    p_f3 = 1 - np.exp(h_du-h_xor) if h_xor > h_du else 0

    p_link_1 = p_f1 + (1 - p_f1) * p_a1
    p_link_2 = p_f2 + (1 - p_f2) * p_a2

    result = 0
    rv_gc = binom(N, 1 - p_f1)
    for gc in range(0, N+1):
        rv_ad = binom(gc, 1 - p_a1)
        for ad in range(0, gc+1):
            rv_ad_tilde = binom(ad, 1 - p_a2)
            bu = gc - ad
            rv_bu = binom(ad, 1-p_link_2)
            kad = np.arange(1, ad+1)
            s_bu_tilde = (1 - p_f3) + p_f3 * np.dot(rv_bu.pmf(kad), 1-np.power(p_f3, kad)) if p_f3 else 1
            # s_bu_tilde = (1 - p_f3) + p_f3 * sum([rv_bu.pmf(k) * (1 - p_f3**k) for k in range(1, ad+1)]) if p_f3 else 1
            q_bu_tilde = (1 - s_bu_tilde) + s_bu_tilde * p_a3
            for ad_tilde in range(0, ad+1):
                rv_ad_tilde_s = binom(ad_tilde, 1 - p_f3)
                rv_ad_hat_s = binom(ad - ad_tilde, 1 - p_f3)
                # if p_f3 = 0 then ad_tilde_s should = ad_tilde because ad_tilde_i should be empty
                for ad_tilde_s in range(0 if p_f3 else ad_tilde, ad_tilde+1):
                    # ad_tilde already succeeded
                    # if p_f3 = 0 then ad_hat_s should = ad - ad_tilde because ad_hat_i should be empty
                    for ad_hat_s in range(0 if p_f3 else ad-ad_tilde, ad-ad_tilde+1):
                        ad_s = ad_tilde_s + ad_hat_s
                        ad_i = ad - ad_tilde_s - ad_hat_s # ad - ad_s
                        rv_ads = binom(ad_s, 1 - p_link_2)
                        rv_adi = binom(ad_i, 1 - p_link_2)
                        q_ad_hat_s = p_a3
                        q_ad_hat_i = p_link_2**ad_s + (1 - p_link_2)**ad_s * p_a3
                        # This is the problem zone that doesn't support vectorization (?)
                        ks, ki = np.arange(1, ad_s+1), np.arange(0, ad_i+1)
                        f_e = p_link_2**ad_s + (np.dot(rv_ads.pmf(ks), np.power(p_f3, ks)) * np.dot(rv_adi.pmf(ki), np.power(p_f3, ki)) if p_f3 else 0)
                        q_e = f_e + (1 - f_e) * (1 - (1 - p_a3)**2)
                        rv_bu_tilde = binom(bu, 1 - p_a2)
                        bu_tilde = np.arange(0, bu+1)
                        log_pstates = np.log10(rv_gc.pmf(gc)) + np.log10(rv_ad.pmf(ad)) + np.log10(rv_ad_tilde.pmf(ad_tilde)) + np.log10(rv_ad_tilde_s.pmf(ad_tilde_s)) + np.log10(rv_ad_hat_s.pmf(ad_hat_s))
                        bu_tilde_res = np.dot(rv_bu_tilde.pmf(bu_tilde), np.multiply(np.power(1 - q_e, N-ad-bu_tilde), np.power(1 - q_bu_tilde, bu_tilde)))
                        res = 10**log_pstates * (1 - q_ad_hat_s)**ad_hat_s * (1 - q_ad_hat_i)**(ad - ad_tilde - ad_hat_s) * bu_tilde_res
                        result += res
                        # print res, result
                        # print gc, ad, ad_tilde, ad_tilde_s, ad_hat_s, bu_tilde
                        # print "\n"
    return result

# def xor_analysis_new(N, tSNR, rSNRdu, rSNR3, p_a1=10**(-9), p_a2=10**(-9), p_a3=10**(-9)):
#     if rSNRdu > tSNR: return 0
#     h_du = 10**((rSNRdu - tSNR)/10) # linear fade
#     h_xor = 10**((rSNR3 - tSNR)/10)
#     # Probability fade is bad
#     p_f1 = p_f2 = 1 - np.exp(-h_du)
#     p_f3 = 1 - np.exp(h_du-h_xor) if h_xor > h_du else 0

#     p_link_1 = p_f1 + (1 - p_f1) * p_a1
#     p_link_2 = p_f2 + (1 - p_f2) * p_a2

#     result = 0
#     rv_gc = binom(N, 1 - p_f1)
#     for gc in range(0, N+1):
#         rv_ad = binom(gc, 1 - p_a1)
#         for ad in range(0, gc+1):
#             rv_ad_tilde = binom(ad, 1 - p_a2)
#             bu = gc - ad
#             rv_bu = binom(ad, 1-p_link_2)
#             kad = np.arange(1, ad+1)
#             s_bu_tilde = (1 - p_f3) + p_f3 * np.dot(rv_bu.pmf(kad), 1-np.power(p_f3, kad)) if p_f3 else 1
#             # s_bu_tilde = (1 - p_f3) + p_f3 * sum([rv_bu.pmf(k) * (1 - p_f3**k) for k in range(1, ad+1)]) if p_f3 else 1
#             q_bu_tilde = (1 - s_bu_tilde) + s_bu_tilde * p_a3
#             for ad_tilde in range(0, ad+1):
#                 rv_ad_tilde_s = binom(ad_tilde, 1 - p_f3)
#                 rv_ad_hat_s = binom(ad - ad_tilde, 1 - p_f3)
#                 # if p_f3 = 0 then ad_tilde_s should = ad_tilde because ad_tilde_i should be empty
#                 for ad_tilde_s in range(0 if p_f3 else ad_tilde, ad_tilde+1):
#                     # ad_tilde already succeeded
#                     # if p_f3 = 0 then ad_hat_s should = ad - ad_tilde because ad_hat_i should be empty
#                     for ad_hat_s in range(0 if p_f3 else ad-ad_tilde, ad-ad_tilde+1):
#                         q_ad_hat_s = p_a3
#                         ad_s = ad_tilde_s + ad_hat_s
#                         ad_i = (ad_tilde - ad_tilde_s) + (ad - ad_tilde - ad_hat_s)
#                         q_ad_hat_i = p_link_2**ad_s + (1 - p_link_2)**ad_s * p_a3
#                         rv_ads = binom(ad_s, 1 - p_link_2)
#                         rv_adi = binom(ad_i, 1 - p_link_2)
#                         ks, ki = np.arange(1, ad_s+1), np.arange(0, ad_i+1)
#                         f_e = p_link_2**ad_s + (np.dot(rv_ads.pmf(ks), np.power(p_f3, ks)) * np.dot(rv_adi.pmf(ki), np.power(p_f3, ki)) if p_f3 else 0)
#                         # f_e = p_link_2**ad_s + (sum([rv_ads.pmf(ks) * p_f3**ks for ks in range(1, ad_s+1)]) if p_f3 else 0) * (sum([rv_adi.pmf(ki) * p_f3**ki for ki in range(0, ad_i+1)]) if p_f3 else 0)
#                         q_e = f_e + (1 - f_e) * (1 - (1 - p_a3)**2)
#                         rv_bu_tilde = binom(bu, 1 - p_a2)
#                         for bu_tilde in range(0, bu+1):
#                             log_pstates = np.log10(rv_gc.pmf(gc)) + np.log10(rv_ad.pmf(ad)) + np.log10(rv_ad_tilde.pmf(ad_tilde)) + np.log10(rv_bu_tilde.pmf(bu_tilde)) + np.log10(rv_ad_tilde_s.pmf(ad_tilde_s)) + np.log10(rv_ad_hat_s.pmf(ad_hat_s))
#                             # print [np.log10(rv_gc.pmf(gc)) , np.log10(rv_ad.pmf(ad)) , np.log10(rv_ad_tilde.pmf(ad_tilde)) , np.log10(rv_bu_tilde.pmf(bu_tilde)) , np.log10(rv_ad_tilde_s.pmf(ad_tilde_s)) , np.log10(rv_ad_hat_s.pmf(ad_hat_s))]
#                             # print 10**log_pstates, (1 - q_bu_tilde)**bu_tilde
#                             res = 10**log_pstates * (1 - q_ad_hat_s)**ad_hat_s * (1 - q_ad_hat_i)**(ad - ad_tilde - ad_hat_s) * (1 - q_bu_tilde)**bu_tilde * (1-q_e)**(N - ad - bu_tilde)
#                             # log_res = log_pstates + ad_hat_s * np.log10(1 - q_ad_hat_s) + (ad - ad_tilde - ad_hat_s) * np.log10(1 - q_ad_hat_i) + bu_tilde * np.log10(1 - q_bu_tilde) + (N - ad - bu_tilde) * np.log10(1 - q_e)
#                             # print log_pstates, log_res
#                             result += res
#                             # print res, result
#                             # print gc, ad, ad_tilde, ad_tilde_s, ad_hat_s, bu_tilde
#                             # print "\n"
#     return result

xor_counts = [0, 0, 21, 56, 126, 252, 462, 792, 1287, 2002, 3003, 4368, 6188, 8568, 11628, 15504, 20349, 26334, 33649, 42504, 53130, 65780, 80730, 98280, 118755, 142506, 169911, 201376, 237336, 278256, 324632, 376992, 435897, 501942, 575757, 658008]

def xor_analysis_pruned(N, tSNR, rSNRdu, rSNR3, p_a1=10**(-9), p_a2=10**(-9), p_a3=10**(-9)):
    if rSNRdu > tSNR: return 0
    h_du = 10**((rSNRdu - tSNR)/10) # linear fade
    h_xor = 10**((rSNR3 - tSNR)/10)
    # Probability fade is bad
    p_f1 = p_f2 = 1 - np.exp(-h_du)
    p_f3 = 1 - np.exp(h_du-h_xor) if h_xor > h_du else 0

    p_link_1 = p_f1 + (1 - p_f1) * p_a1
    p_link_2 = p_f2 + (1 - p_f2) * p_a2
    
    # test = []
#     test2 = []
    pstates_thresh = np.log10(10**(-10)/xor_counts[N])
    
    log_pstates = np.zeros((5,))
    
    result = 0
    rv_gc = binom(N, 1 - p_f1)
    for gc in range(0, N+1):
        log_pstates[0] = np.log10(rv_gc.pmf(gc))
        rv_ad = binom(gc, 1 - p_a1)
        for ad in range(0, gc+1):
            log_pstates[1] = np.log10(rv_ad.pmf(ad))
            if sum(log_pstates) < pstates_thresh:
                continue
            rv_ad_tilde = binom(ad, 1 - p_a2)
            bu = gc - ad
            rv_bu = binom(ad, 1-p_link_2)
            kad = np.arange(1, ad+1)
            s_bu_tilde = (1 - p_f3) + p_f3 * np.dot(rv_bu.pmf(kad), 1-np.power(p_f3, kad)) if p_f3 else 1
            # s_bu_tilde = (1 - p_f3) + p_f3 * sum([rv_bu.pmf(k) * (1 - p_f3**k) for k in range(1, ad+1)]) if p_f3 else 1
            q_bu_tilde = (1 - s_bu_tilde) + s_bu_tilde * p_a3
            for ad_tilde in range(0, ad+1):
                log_pstates[2] = np.log10(rv_ad_tilde.pmf(ad_tilde))
                if sum(log_pstates) < pstates_thresh:
                    continue
                rv_ad_tilde_s = binom(ad_tilde, 1 - p_f3)
                rv_ad_hat_s = binom(ad - ad_tilde, 1 - p_f3)
                # if p_f3 = 0 then ad_tilde_s should = ad_tilde because ad_tilde_i should be empty
                for ad_tilde_s in range(0 if p_f3 else ad_tilde, ad_tilde+1):
                    log_pstates[3] = np.log10(rv_ad_tilde_s.pmf(ad_tilde_s))
                    if sum(log_pstates) < pstates_thresh:
                        continue
                    # ad_tilde already succeeded
                    # if p_f3 = 0 then ad_hat_s should = ad - ad_tilde because ad_hat_i should be empty
                    for ad_hat_s in range(0 if p_f3 else ad-ad_tilde, ad-ad_tilde+1):
                        log_pstates[4] = np.log10(rv_ad_hat_s.pmf(ad_hat_s))
#                         if [gc, ad, ad_tilde, ad_tilde_s, ad_hat_s] == [3, 3, 2, 2, 1]:
#                             print log_pstates
                        if sum(log_pstates) < pstates_thresh:
                            continue
                        ad_s = ad_tilde_s + ad_hat_s
                        ad_i = ad - ad_tilde_s - ad_hat_s # ad - ad_s
                        rv_ads = binom(ad_s, 1 - p_link_2)
                        rv_adi = binom(ad_i, 1 - p_link_2)
                        q_ad_hat_s = p_a3
                        q_ad_hat_i = p_link_2**ad_s + (1 - p_link_2)**ad_s * p_a3
                        # This is the problem zone that doesn't support vectorization (?)
                        ks, ki = np.arange(1, ad_s+1), np.arange(0, ad_i+1)
                        f_e = p_link_2**ad_s + (np.dot(rv_ads.pmf(ks), np.power(p_f3, ks)) * np.dot(rv_adi.pmf(ki), np.power(p_f3, ki)) if p_f3 else 0)
                        q_e = f_e + (1 - f_e) * (1 - (1 - p_a3)**2)
                        rv_bu_tilde = binom(bu, 1 - p_a2)
                        bu_tilde = np.arange(0, bu+1)
#                         log_pstates = np.log10(rv_gc.pmf(gc)) + np.log10(rv_ad.pmf(ad)) + np.log10(rv_ad_tilde.pmf(ad_tilde)) + np.log10(rv_ad_tilde_s.pmf(ad_tilde_s)) + np.log10(rv_ad_hat_s.pmf(ad_hat_s))
                        bu_tilde_res = np.dot(rv_bu_tilde.pmf(bu_tilde), np.multiply(np.power(1 - q_e, N-ad-bu_tilde), np.power(1 - q_bu_tilde, bu_tilde)))
                        res = 10**sum(log_pstates) * (1 - q_ad_hat_s)**ad_hat_s * (1 - q_ad_hat_i)**(ad - ad_tilde - ad_hat_s) * bu_tilde_res
                        result += res
                        # print res, result
                        # print gc, ad, ad_tilde, ad_tilde_s, ad_hat_s, bu_tilde
                        # print "\n"
                        # test.append([gc, ad, ad_tilde, ad_tilde_s, ad_hat_s, res, sum(log_pstates)])
#                         test.append([res, log_pstates, (1 - q_ad_hat_s)**ad_hat_s, (1 - q_ad_hat_i)**(ad - ad_tilde - ad_hat_s), bu_tilde_res])
    return result

def xor_analysis_opt(N, tSNR, rSNR1, rSNR2, p_add_3, p_add_1=10**(-9), p_add_2=10**(-9)):
    if rSNR1 > tSNR or rSNR2 > tSNR:
        return 0
    hcrit_1 = 10**((rSNR1 - tSNR)/10) # linear fade
    p_fade_1 = 1 - np.exp(-hcrit_1)
    hcrit_2 = 10**((rSNR2 - tSNR)/10) # linear fade
    p_fade_2 = 1 - np.exp(-hcrit_2)

    p_link_1 = p_fade_1 + (1 - p_fade_1) * p_add_1
    p_link_2 = p_fade_2 + (1 - p_fade_2) * p_add_2
    q = p_add_1 / p_link_1

    result = 0
    for ad in range(1, N+1):
        p_ad = nCr(N, ad) * (1 - p_link_1)**ad * p_link_1**(N-ad)
        adn_result = 0
        for adu in range(1, ad+1):
            p_adu = nCr(ad, adu) * (1 - p_add_2)**adu * p_add_2**(ad-adu)
            p_adn_success = (1 - p_add_3)**(ad - adu)
            adn_result += p_adu * p_adn_success

        b_result = 0
        for b in range(0, N-ad+1):
            bu_result = 0
            for bu in range(0, b+1):
                p_bu = nCr(b, bu) * (1 - p_add_2)**bu * p_add_2**(b-bu)
                p_bu_success = (1 - p_add_3)**bu
                p_rest_success = ((1 - p_link_2**ad) * (1 - p_add_3)**2)**(N-ad-bu)
                bu_result += p_bu * p_bu_success * p_rest_success

            p_b = nCr(N-ad, b) * q**b * (1-q)**(N-ad-b)
            b_result += p_b*bu_result

        result += p_ad * adn_result * b_result
    return result

def optimize_1(N, tSNR_start, filepath_down, filepath_up, precision, protocol=4*10**4, downtarget=10**(-10), uptarget=10**(-10)):
    """The most naive form of optimization, where we assume the blocklength is divided evenly between all three phases.

    We also set the minimum SNR experienced at the receiver to be such that all three phases meet their additive noise targets. Additive noise is defined to be the noise other than fade that would cause the error correcting code to be corrupted.

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      tx_SNR_range {np.arange} -- [description]
      filepath_down {string} -- [description]
      filepath_up {string} -- [description]

    Keyword Arguments:
      protocol {int} -- The total blocklength in bits for all phases combined (default: {40,000})
      downtarget {float: fraction} -- [description] (default: {10**(-10)})
      uptarget {float: fraction} -- [description] (default: {10**(-10)})

    Returns:
      np.array of length 7 --
      0. The SNR sent at the transmitter (pre-fade)
      1. The SNR threshold at the downlink receiver (post-fade)
      2. The SNR threshold at the uplink receiver (post-fade)
      3. The SNR needed at the XOR receiver
      4. The additive noise in downlink
      5. The additive noise in uplink
      6. The additive noise in the XOR phase
    """
    p_add_3 = 10**(-10)

    downbit = upbit = xorbit = protocol/3
    downTable = load_table(filepath_down+str(N)+extension)
    downfunc = interp1d(downTable[0], downTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    upTable = load_table(filepath_up+str(N)+extension)
    upfunc = interp1d(upTable[0], upTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    xor_table = load_table(filepath_up+str(N)+extension)
    xor_func = interp1d(xor_table[0], xor_table[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    down_SNR = downTable[0][np.where(np.array(downTable[1])<=downtarget)[0][0]]
    up_SNR = upTable[0][np.where(np.array(upTable[1])<=uptarget)[0][0]]
    xor_SNR = xor_table[0][np.where(np.array(xor_table[1])<=p_add_3)[0][0]]
    rSNRdu = max(down_SNR, up_SNR)
    # actual_SNR = max(down_SNR, up_SNR, xor_SNR)
    
    for dSNR in 10**np.arange(0.0, precision-1, -1):
        for tSNR in np.arange(tSNR_start, 90, dSNR):
            xor_opt = xor_analysis_new(N, tSNR, rSNRdu, xor_SNR, downfunc(down_SNR), upfunc(up_SNR), xor_func(xor_SNR))
            if 1-xor_opt <= protocol_target:
                tSNR_start = tSNR - dSNR
                break
    return (tSNR, down_SNR, up_SNR, xor_SNR, downfunc(down_SNR), upfunc(up_SNR), xor_func(xor_SNR))

def optimize_3(N, rSNRdu_start, rSNR3_start, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, d_rSNR=0.1, precision=-2):
    """We assume the blocklength is evenly divided among all three phases.
    We don't enforce any Downlink or Uplink additive noise targets. Instead, the protocol finds the "right" transmitter SNR and receiver SNR pair such that the overall protocol meets its reliability target.

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      rx_SNR_start {int} -- The number to start the Receiver SNR range (ends at Transmitter SNR)
      tx_SNR_range {np.arange} -- [description]
      filepath_down {string} -- [description]
      filepath_up {string} -- [description]

    Keyword Arguments:
      protocol {int} -- The length of the entire protocol (all phases combined) in bits (default: {40,000})

    Returns:
      np.array of length 5 --
      0. The SNR sent at the transmitter (pre-fade) -- this is what we care about minimizing
      1. The SNR experienced at the receiver (post-fade)
      2. The blocklength of the Downlink Phase in bits
      3. The blocklength of the Uplink Phase in bits
      4. The blocklength of the XOR Phase in bits
    """
    downbit = upbit = xorbit = protocol/3
    downTable = load_table(filepath_down+str(N)+extension)
    downfunc = interp1d(downTable[0], downTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    upTable = load_table(filepath_up+str(N)+extension)
    upfunc = interp1d(upTable[0], upTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    xor_table = load_table(filepath_up+str(N)+extension)
    xor_func = interp1d(xor_table[0], xor_table[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    # tSNR_start = tx_SNR_range[0]
    # for dSNR in 10**np.arange(0.0, precision-1, -1):
    for tSNR in tx_SNR_range:
        rSNRdu_test = [float("inf")]
        for rSNR3 in np.arange(rSNR3_start, min(15, tSNR), d_rSNR):
            p_add_3 = xor_func(rSNR3)
            xoropt_last = 0.0
            for rSNRdu in np.arange(rSNRdu_start, rSNR3+d_rSNR, d_rSNR):
                xor_opt = xor_analysis_pruned(N, tSNR, rSNRdu, rSNR3, downfunc(rSNRdu), upfunc(rSNRdu), xor_func(rSNR3))
                if rSNRdu == rSNRdu_start:
                    if 1 - xor_opt > rSNRdu_test[-1]:
                        break
                    rSNRdu_test.append(1 - xor_opt)
                # xor_opt = xor_analysis_opt(N, tSNR, rSNR1, rSNR2, p_add_3, downfunc(rSNR1), upfunc(rSNR2))
                if 1-xor_opt <= protocol_target:
                    print 1 - xor_opt
                    return (tSNR, rSNRdu, rSNR3, downfunc(rSNRdu), upfunc(rSNRdu), p_add_3)
                if xor_opt < xoropt_last:
                    break
                else:
                    xoropt_last = xor_opt
        # print rSNRdu_test
    # return (float("inf"), float("inf"), downbit, upbit, xorbit)

def optimize_3_graphical(N, rSNRdu_range, rSNR3_range, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, inner=True):
    downbit = upbit = xorbit = protocol/3
    downTable = load_table(filepath_down+str(N)+extension)
    downfunc = interp1d(downTable[0], downTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    upTable = load_table(filepath_up+str(N)+extension)
    upfunc = interp1d(upTable[0], upTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    xor_table = load_table(filepath_up+str(N)+extension)
    xor_func = interp1d(xor_table[0], xor_table[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    if inner:
        xn, yn = len(rSNR3_range), len(rSNRdu_range)
    else:
        xn, yn = len(tx_SNR_range), len(rSNR3_range)

    result = np.zeros((xn, yn))

    print "xn, yn", result.shape, "size of matrix", xn * yn

    for tSNR in tx_SNR_range:
        for rSNR3 in rSNR3_range:
            p_add_3 = xor_func(rSNR3)
            for rSNRdu in rSNRdu_range:
                xor_opt = xor_analysis_pruned(N, tSNR, rSNRdu, rSNR3, downfunc(rSNRdu), upfunc(rSNRdu), xor_func(rSNR3))
                if inner:
                    xk, yk = np.where(rSNR3_range == rSNR3)[0][0], np.where(rSNRdu_range == rSNRdu)[0][0]
                else:
                    xk, yk = np.where(tx_SNR_range == tSNR)[0][0], np.where(rSNR3_range == rSNR3)[0][0]
                result[xk, yk] = 1 - xor_opt
                # print result[xk, yk]
    return result

def optimize_4(N, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, downtarget=10**(-9), uptarget=10**(-9), d_rSNR=0.1):
    """We enforce a Downlink additive noise target and an Uplink additive noise target. We do NOT assume the blocklength is evenly divided among all three phases. Instead, we allocate the minimum blocklength so that Downlink meets its additive noise target. Then we allocate the minimum blocklength so that Uplink meets its additive noise target. The remaining blocklength is allocated to the XOR phase, which determines the XOR additive noise. The optimization module then finds the transmitter SNR and receiver SNR pair so that the combination of parameters will allow the protocol to meet its overall reliability target.

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      tx_SNR_range {np.arange} -- [description]
      filepath_down {string} -- [description]
      filepath_up {string} -- [description]

    Keyword Arguments:
      protocol {int} -- The length of the entire protocol (all phases combined) in bits (default: {40,000})
      downtarget {float: fraction} -- [description] (default: {10**(-9)})
      uptarget {float: fraction} -- [description] (default: {10**(-9)})

    Returns:
      np.array of length 5 --
      0. The SNR experienced at the receiver (post-fade)
      1. The SNR sent at the transmitter (pre-fade) -- this is what we care about minimizing
      2. The blocklength of the Downlink Phase in bits
      3. The blocklength of the Uplink Phase in bits
      4. The blocklength of the XOR Phase in bits
    """
    downFile = filepath_down + str(N) + '.in'
    downNode = load_table(downFile)
    upFile = filepath_up + str(N) + '.in'
    upNode = load_table(upFile)

    for tSNR in tx_SNR_range:
        for rSNR3 in np.arange(-1, tSNR, d_rSNR):
            for rSNR1 in np.arange(-1, tSNR, d_rSNR):
                for rSNR2 in np.arange(-1, tSNR, d_rSNR):
                    downbit, upbit = float("inf"), float("inf")
                    for bit in sorted(downNode.tables.keys()):
                        bittable = downNode.tables[bit]
                        downfunc = interp1d(bittable[0], bittable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
                        if downfunc(rSNR1) <= downtarget:
                            downbit = bit
                            break
                    for bit in sorted(upNode.tables.keys()):
                        bittable = upNode.tables[bit]
                        upfunc = interp1d(bittable[0], bittable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
                        if upfunc(rSNR2) <= uptarget:
                            upbit = bit
                            break
                    xorbit = protocolsla - downbit - upbit
                    if xorbit <= 0: continue
                    # We calculate p_add_3 each time because the rate changes every time (new table)
                    blocklength = int(xorbit / 21 / N)
                    rate = N * 160 / xorbit * 7 / 4 if xorbit else float("inf")
                    if rate > 1: continue
                    k = (1-rate)*blocklength
                    pbitdrop = Q(np.sqrt(2*10**(rSNR3/10)))
                    hcerr = 1 - ((1-pbitdrop)**7 + 7*pbitdrop*(1-pbitdrop)**6)
                    hcf = 1 - (1-hcerr)**3
                    p_add_3 = 1-binom.cdf(int(k/2), blocklength, hcf)
                    # {p_add_3} reeddrop = sum([nCr(blocklength, d)* hcf**d *(1-hcf)**(blocklength-d) for d in range(int(k/2)+1, blocklength)])
                    xor_opt = xor_analysis_opt(N, tSNR, rSNR1, rSNR2, p_add_3, downfunc(rSNR1), upfunc(rSNR2))
                    if 1-xor_opt <= protocol_target:
                        return (tSNR, rSNR1, rSNR2, rSNR3, downbit, upbit, xorbit, downfunc(rSNR1), upfunc(rSNR2), p_add_3)

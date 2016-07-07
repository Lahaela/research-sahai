from __future__ import division
import numpy as np
import scipy as sp
from scipy.stats import binom
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

def xor_analysis_opt(N, p_add_3, nom_SNR, actual_SNR, p_add_1=10**(-9), p_add_2=10**(-9)):
    hcrit_1 = hcrit_2 = 10**((actual_SNR - nom_SNR)/10) # linear fade
    p_fade_1 = 1 - np.exp(-hcrit_1)
#     hcrit_2 = 10**((actual_SNR - nom_SNR)/10) # linear fade
    p_fade_2 = 1 - np.exp(-hcrit_2)
#     p_fade_1 = p_fade_2 = 1 - np.exp(-1)
    p_link_1 = p_fade_1 + (1 - p_fade_1) * p_add_1
    p_link_2 = p_fade_2 + (1 - p_fade_2) * p_add_2
    q = p_add_1 / p_link_1
#     p_add_3 = xorfunc(SNR)

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

def optimize_1(N, rx_SNR_range, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, downtarget=10**(-10), uptarget=10**(-10)):
    """The most naive form of optimization, where we assume the blocklength is divided evenly between all three phases.

    We also set the minimum SNR experienced at the receiver to be such that all three phases meet their additive noise targets. Additive noise is defined to be the noise other than fade that would cause the error correcting code to be corrupted.

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      rx_SNR_range {np.arange} -- [description]
      tx_SNR_range {np.arange} -- [description]
      filepath_down {string} -- [description]
      filepath_up {string} -- [description]

    Keyword Arguments:
      protocol {int} -- The total blocklength in bits for all phases combined (default: {40,000})
      downtarget {float: fraction} -- [description] (default: {10**(-10)})
      uptarget {float: fraction} -- [description] (default: {10**(-10)})

    Returns:
      np.array of length 5 --
      0. The SNR experienced at the receiver (post-fade)
      1. The SNR sent at the transmitter (pre-fade) -- this is what we care about minimizing
      2. The blocklength of the Downlink Phase in bits
      3. The blocklength of the Uplink Phase in bits
      4. The blocklength of the XOR Phase in bits
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
    actual_SNR = max(down_SNR, up_SNR, xor_SNR)

    for nominal_SNR in tx_SNR_range:
        xor_opt = xor_analysis_opt(N, xor_func(actual_SNR), nominal_SNR, actual_SNR, downfunc(actual_SNR), upfunc(actual_SNR))
        if 1-xor_opt <= protocol_target:
            return (actual_SNR, nominal_SNR, downbit, upbit, xorbit)
    return (float("inf"), float("inf"), downbit, upbit, xorbit)

def optimize_2(N, rx_SNR_start, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, d_actual_SNR=0.1, downtarget=10**(-9), uptarget=10**(-9)):
    """We assume the blocklength is divided evenly among all three phases.
    The receiver SNR is set to be the minimum such that both Downlink and Uplink meet their additive noise targets, but we do not set a XOR additive noise target. Instead, the receiver SNR is also the minimum such that given the XOR blocklength (protocol / 3), the transmitter SNR (which is also optimized over), and receiver SNR, the protocol finds the "right" XOR additive noise target so that the protocol meets its overall protocol reliability target.

    [description]

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      rx_SNR_range {np.arange} -- [description]
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
    downbit = upbit = xorbit = protocol/3
    downTable = load_table(filepath_down+str(N)+extension)
    downfunc = interp1d(downTable[0], downTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
    upTable = load_table(filepath_up+str(N)+extension)
    upfunc = interp1d(upTable[0], upTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    xor_table = load_table(filepath_up+str(N)+extension)
    xor_func = interp1d(xor_table[0], xor_table[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(rx_SNR_start, nominal_SNR, d_actual_SNR):
            if downfunc(actual_SNR) > downtarget: continue
            if upfunc(actual_SNR) > uptarget: continue

            p_add_3 = xor_func(actual_SNR)
            xor_opt = xor_analysis_opt(N, p_add_3, nominal_SNR, actual_SNR, downfunc(actual_SNR), upfunc(actual_SNR))
            if 1-xor_opt <= protocol_target:
                return (actual_SNR, nominal_SNR, downbit, upbit, xorbit)
    return (float("inf"), float("inf"), downbit, upbit, xorbit)

def optimize_3(N, rx_SNR_start, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, d_actual_SNR=0.1):
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
      0. The SNR experienced at the receiver (post-fade)
      1. The SNR sent at the transmitter (pre-fade) -- this is what we care about minimizing
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

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(rx_SNR_start, nominal_SNR, d_actual_SNR):
            p_add_3 = xor_func(actual_SNR)
            xor_opt = xor_analysis_opt(N, p_add_3, nominal_SNR, actual_SNR, downfunc(actual_SNR), upfunc(actual_SNR))
            if 1-xor_opt <= protocol_target:
                return (actual_SNR, nominal_SNR, downfunc(actual_SNR), upfunc(actual_SNR), p_add_3)
    return (float("inf"), float("inf"), downbit, upbit, xorbit)


def optimize_4(N, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4, downtarget=10**(-9), uptarget=10**(-9)):
    """We enforce a Downlink additive noise target and an Uplink additive noise target. We do NOT assume the blocklength is evenly divided among all three phases. Instead, we allocate the minimum blocklength so that Downlink meets its additive noise target. Then we allocate the minimum blocklength so that Uplink meets its additive noise target. The remaining blocklength is allocated to the XOR phase, which determines the XOR additive noise. The optimization module then finds the transmitter SNR and receiver SNR pair so that the combination of parameters will allow the protocl to meet its overall reliability target.

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

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(max(-1, nominal_SNR-90), nominal_SNR, 0.1):
            downbit, upbit = float("inf"), float("inf")
            for bit in sorted(downNode.tables.keys()):
                bittable = downNode.tables[bit]
                downfunc = interp1d(bittable[0], bittable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
                if downfunc(actual_SNR) <= downtarget:
                    downbit = bit
                    break
            for bit in sorted(upNode.tables.keys()):
                bittable = upNode.tables[bit]
                upfunc = interp1d(bittable[0], bittable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
                if upfunc(actual_SNR) <= uptarget:
                    upbit = bit
                    break
            xorbit = protocol - downbit - upbit
            if xorbit <= 0: continue
#             xorbit = max(0, 4200 - downbit - upbit)
            # We calculate reeddrop each time because the rate changes every time (new table)
            blocklength = int(xorbit / 21 / N)
            rate = N * 160 / xorbit * 7 / 4 if xorbit else float("inf")
            if rate > 1: continue
            k = (1-rate)*blocklength
            pbitdrop = Q(np.sqrt(2*10**(actual_SNR/10)))
            hcerr = 1 - ((1-pbitdrop)**7 + 7*pbitdrop*(1-pbitdrop)**6)
            hcf = 1 - (1-hcerr)**3
            reeddrop = 1-binom.cdf(int(k/2), blocklength, hcf)
            # reeddrop = sum([nCr(blocklength, d)* hcf**d *(1-hcf)**(blocklength-d) for d in range(int(k/2)+1, blocklength)])
            xor_opt = xor_analysis_opt(N, reeddrop, nominal_SNR, actual_SNR, downfunc(actual_SNR), upfunc(actual_SNR))
            if 1-xor_opt <= protocol_target:
                return (actual_SNR, nominal_SNR, downbit, upbit, xorbit, downfunc(actual_SNR), upfunc(actual_SNR), reeddrop)
    return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan) # default behavior when nothing is returned

# At the moment, not using this
def optimize_5(N, tx_SNR_range, filepath_down, filepath_up, protocol=4*10**4):
    """We do NOT assume blocklength is evenly divided among the three phases. We don't enforce additive noise constraints on Downlink/Uplink/XOR, instead letting the optimization find the right Downlink/Uplink additive noise. As of this implementation we assume Downlink/Uplink has the same additive noise constraint but this can change in the future. We find the right transmitter SNR and receiver SNR pair and blocklength allocation for the overall protocol to meet its reliability target.

    Arguments:
      N {int} -- The number of nodes/users total in the control system
      tx_SNR_range {np.arange} -- [description]
      filepath_down {string} -- [description]
      filepath_up {string} -- [description]

    Keyword Arguments:
      protocol {int} -- The length of the entire protocol (all phases combined) in bits (default: {40,000})

    Returns:
      [type] -- [description]
    """
    linktarget = [10**i for i in range(-10, -3)]
    optimize_target = zeros((7, 6))
    for idx in range(7):
        optimize_target[idx] = concatenate((optimize_4(N, tx_SNR_range, filepath_down, filepath_up, protocol, linktarget[idx], linktarget[idx]), array([linktarget[idx]])))
        print optimize_target[idx]
    return optimize_target[np.argmin(optimize_target[:,1])]

def optimize_3_downlink(N, rx_SNR_start, tx_SNR_range, filepath_down, protocol=4*10**4, d_actual_SNR=0.1):
    # downbit = protocol/3
    downTable = load_table(filepath_down+str(N)+extension)
    downfunc = interp1d(downTable[0], downTable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(rx_SNR_start, nominal_SNR, d_actual_SNR):
            padd = downfunc(actual_SNR)
            hcrit = 10**((actual_SNR - nominal_SNR)/10) # linear fade
            pbadfade = 1 - np.exp(-hcrit)
            psingle = pbadfade + (1-pbadfade)*padd
            rv = binom(N, psingle)
            pprotocol = sum([rv.pmf(a) * (1-(1-(pbadfade**a + (1-pbadfade**a) * padd))**(N-a)) for a in range(N)])
            if pprotocol <= protocol_target: return (actual_SNR, nominal_SNR, downfunc(actual_SNR))
    return (float("inf"), float("inf"), 0)

def optimize_3_2(N, rx_SNR_start, tx_SNR_range, filepath_down, filepath_up, p_add_1, p_add_2, protocol=4*10**4, d_actual_SNR=0.1):
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
      0. The SNR experienced at the receiver (post-fade)
      1. The SNR sent at the transmitter (pre-fade) -- this is what we care about minimizing
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

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(rx_SNR_start, nominal_SNR, d_actual_SNR):
            p_add_3 = xor_func(actual_SNR)
            xor_opt = xor_analysis_opt(N, p_add_3, nominal_SNR, actual_SNR, p_add_1, p_add_2)
            if 1-xor_opt <= protocol_target:
                return (actual_SNR, nominal_SNR, p_add_1, p_add_2, p_add_3)
    return (float("inf"), float("inf"), downbit, upbit, xorbit)

def optimize_4_shannon(N, rx_SNR_start, tx_SNR_range, filepath, protocol=4*10**4, phasetarget=10**(-9)):
    """We enforce a Downlink additive noise target and an Uplink additive noise target. We do NOT assume the blocklength is evenly divided among all three phases. Instead, we allocate the minimum blocklength so that Downlink meets its additive noise target. Then we allocate the minimum blocklength so that Uplink meets its additive noise target. The remaining blocklength is allocated to the XOR phase, which determines the XOR additive noise. The optimization module then finds the transmitter SNR and receiver SNR pair so that the combination of parameters will allow the protocl to meet its overall reliability target.

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
    shannonFile = filepath + str(N) + extension
    shannonNode = load_table(shannonFile)

    for nominal_SNR in tx_SNR_range:
        for actual_SNR in np.arange(rx_SNR_start, nominal_SNR, 0.1):
            downbit, upbit = float("inf"), float("inf")
            for bit in shannonNode.bitrange:
                bittable = shannonNode.tables[bit]
                func = interp1d(bittable[0], bittable[1], kind='linear', bounds_error=False, fill_value=(1.0, 0.0))
                if func(actual_SNR) <= phasetarget:
                    downbit = upbit = bit
                    break
            xorbit = protocol - downbit - upbit
            if xorbit <= 0: continue
#             xorbit = max(0, 40,000 - downbit - upbit)
            # We calculate reeddrop each time because the rate changes every time (new table)
            rate = N * 160 / xorbit if xorbit else float("inf")
            if rate > 1: continue
            C = np.log2(1 + 10**(actual_SNR/10))
            reeddrop = 1.0 if rate > C else 0.0
            xor_opt = xor_analysis_opt(N, reeddrop, nominal_SNR, actual_SNR, func(actual_SNR), func(actual_SNR))
            if 1-xor_opt <= protocol_target:
                return (actual_SNR, nominal_SNR, downbit, upbit, xorbit, func(actual_SNR), func(actual_SNR), reeddrop)
    return (float("inf"), float("inf"), float("inf"), float("inf"), float("inf"), 0, 0, 0) # default behavior when nothing is returned

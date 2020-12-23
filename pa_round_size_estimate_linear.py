"""
Quick script to get a first round estimate for a Minerva audit
of the PA presidential race in 2020.

Oliver Broadrick 2020
"""

import time
import numpy as np
import scipy as sp
import scipy.stats
import scipy.optimize
from ballot_comparison import ballot_comparison_pvalue
from sprt import ballot_polling_sprt
import matplotlib.pyplot as plt
import numpy.testing
from contest import ContestType
from contest import Contest
from minerva_s import Minerva_S
from fishers_combination import create_modulus, maximize_fisher_combined_pvalue, calculate_lambda_range, maximize_stouffers_combined_pvalue
from scipy.stats import binom
import math
import matplotlib.pyplot as plt
from simulations import minerva_pvalue_direct_count, r2bravo_pvalue_direct_count


# pa nums here (2020 presidential race)
Nw = 3459923 #biden
Nl = 3378263 #trump
N = Nw + Nl
margin = Nw - Nl
stopping_probability = .9 #prob stop under alt hyp, 90% in first round
alpha = .1

# binary search (jagged graph so not perfect, but quick for estimate)
# binary search bounds 
n = 50000
while(1):
    n = n - 1

    # compute the 1 - stopping_probability quantile of the alt dist
    # kmax where pr[k >= kmax | alt] = stopping_probability
    # floor because we need to ensure at least a stopping_probability prob of stopping
    kmax = math.floor(binom.ppf(1 - stopping_probability, n, Nw / N))

    pvalue = minerva_pvalue_direct_count(winner_votes=kmax, n=n, popsize=N, alpha=alpha, \
                            Vw=Nw, Vl=Nl, null_margin=0)

    # if left = right then the initial right bound was too small

    # when and right converge, right is the minimum round size that achieves stopping_probability
    if (pvalue > alpha):
        #print(combination_results['refined'])
        print("pvalue(one below):",pvalue)
        print("round_size:",n)
        break



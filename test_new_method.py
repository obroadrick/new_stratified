"""
Script to test the new stratified method.

Oliver Broadrick 2020
"""
from new_method import maximize_joint_pvalue, find_minimum_round_size
import time
from scipy.stats import binom
from round_sizes import find_sample_size_for_stopping_prob_efficiently, \
                find_sample_size_for_stopping_prob_efficiently_r2bravo

# set up a contest
N_relevant = 100
fractional_margin = .2
polling_proportion = .5
N_w = round(N_relevant * (1 + fractional_margin) / 2)
N_l = N_relevant - N_w
N_w_fraction = N_w / N_relevant
N_2 = round(polling_proportion * N_relevant)
N_1 = N_relevant - N_2
N_w1 = round(N_w_fraction * N_1)
N_l1 = N_1 - N_w1
N_w2 = N_w - N_w1
N_l2 = N_2 - N_w2

"""
###############################################################################
## COMPUTE PVALUE
###############################################################################
# define some samples
n1 = 10
n2 = 10
k_c = n1
k_p = int(binom.ppf(.9, n2, N_w2 / (N_w2 + N_l2)))
print("kmax:",k_p)
 
# compute pvalue (tracking time)
start = time.time()
results = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None, plot=True, print_data=False)
print("time:",(time.time()-start)/60,"minutes")
print(results)

"""

#"""
###############################################################################
## COMPUTE MINIMUM POLLING ROUND SIZE
###############################################################################
# comparison stratum sample size
n1 = 10
# desired stopping probability
stop_prob = .9
# risk limit
alpha = .1

# Print beginning message
print("Finding minimum polling round sizes that give .9 stopping probablility. . .")

# Find minimum round size with new method
start = time.time()
results = find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha)
#print("time:",(time.time()-start)/60,"minutes")
print("New method")
print(results['round_size'])
print(results['pvalue'])

# Also get suite round sizes (with minerva)
suite_minerva = find_sample_size_for_stopping_prob_efficiently(stop_prob, \
        N_w1, N_l1, N_w2, N_l2, n1, alpha)
print("\nSUITE with Minerva")
print("round size:", suite_minerva['round_size'])
print("pvalue:",suite_minerva['combined_pvalue'])

#"""





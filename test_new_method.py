"""
Script to test the new stratified method.

Oliver Broadrick 2020
"""
from new_method import maximize_joint_pvalue, find_minimum_round_size
from new_method_na import maximize_joint_pvalue_na, find_minimum_round_size_na
import time
from scipy.stats import binom
from round_sizes import find_sample_size_for_stopping_prob_efficiently, \
                find_sample_size_for_stopping_prob_efficiently_r2bravo

# set up a contest
N_relevant = 400
fractional_margin = .3
polling_proportion = .8
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
n1 = 20
# desired stopping probability
stop_prob = .9
# risk limit
alpha = .1

# Print beginning message
print("Consider an audit of a contest with the following values:")
print(N_relevant, "relevant ballots,")
print(N_w1, "winner ballots in the comparison stratum,")
print(N_l1, "loser ballots in the comparison stratum,")
print(N_w2, "winner ballots in the polling stratum,")
print(N_l2, "loser ballots in the polling stratum,")
print("a risk limit of", alpha)
print("and a comparison sample of", n1, "matches.")
print()
print("Finding polling stratum minimum round sizes that give .9 stopping probablility. . .")

# New method
start = time.time()
results = find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha)
#print("time:",(time.time()-start)/60,"minutes")
print("New method")
print(results['round_size'], "minimum round size")
print(results['pvalue'], "pvalue for kmax", results['kmax'])

results_na = find_minimum_round_size_na(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha)
#print("time:",(time.time()-start)/60,"minutes")
print("\nNew method not athenized")
print(results_na['round_size'], "minimum round size")
print(results_na['pvalue'], "pvalue for kmax", results['kmax'])


# Minerva SUITE
suite_minerva = find_sample_size_for_stopping_prob_efficiently(stop_prob, \
        N_w1, N_l1, N_w2, N_l2, n1, alpha, right=1000)
print("\nSUITE with Minerva")
print(suite_minerva['round_size'], "minimum round size")
print(suite_minerva['combined_pvalue'], "pvalue for kmax")

# Minerva SUITE
suite_r2bravo = find_sample_size_for_stopping_prob_efficiently_r2bravo(stop_prob, \
        N_w1, N_l1, N_w2, N_l2, n1, alpha, right = 1000)
print("\nSUITE with R2 Bravo")
print(suite_r2bravo['round_size'], "minimum round size")
print(suite_r2bravo['combined_pvalue'], "pvalue for kmax")

#"""





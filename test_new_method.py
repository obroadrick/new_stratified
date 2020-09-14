"""
Script to test the new stratified method.

Oliver Broadrick 2020
"""
from new_method import maximize_joint_pvalue, find_minimum_round_size
import time
from scipy.stats import binom

# set up a contest
N_relevant = 104000
fractional_margin = .0192
percent_polling = .05
N_w = round(N_relevant * (1 + fractional_margin) / 2)
N_l = N_relevant - N_w
N_w_fraction = N_w / N_relevant
N_2 = round(percent_polling * N_relevant)
N_1 = N_relevant - N_2
N_w1 = round(N_w_fraction * N_1)
N_l1 = N_1 - N_w1
N_w2 = N_w - N_w1
N_l2 = N_2 - N_w2

"""
# TEST THE MAXIMIZE PVALUE FUNCTION
# samples
n1 = 750
n2 = 710
k_c = n1
k_p = int(binom.ppf(.9, n2, N_w2 / (N_w2 + N_l2)))
print("kmax:",k_p)
 
# compute pvalue (tracking time)
start = time.time()
results = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None, plot=True, print_data=False)

# print results
print("time:",(time.time()-start)/60,"minutes")
print(results)
"""

# TEST THE FIND MINIMUM ROUND SIZE FUNCTION
# comparison stratum sample size
n1 = 750
# desired stopping probability
stop_prob = .9
# risk limit
alpha = .1

# find minimum round size (tracking time)
start = time.time()
results = find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha)

# print results
print("time:",(time.time()-start)/60,"minutes")
print(results)








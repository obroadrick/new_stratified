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
import numpy as np
import json

# store data in json
data_file_name = "round_sizes.json"
data = {}
data['audits'] = []

# set up a contest
N_relevant = 10000
n1 = 200
# desired stopping probability
stop_prob = .9
# risk limit
alpha = .1

print("\nConsider auditing a contest with", N_relevant, "relevant ballots,")
print("a risk limit of", alpha)
print("and a comparison sample of", n1, "matches.")
print()
print("Now we find the first round minimum sample sizes for the polling stratum that give .9 stopping probablility. . .")


for fractional_margin in np.linspace(.02, .2, 19):
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

    # Print beginning messages
    print("\nfor a margin of",fractional_margin,". . .")

    # New method
    start = time.time()
    results = find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha)
    #print("time:",(time.time()-start)/60,"minutes")
    print("New method")
    print(results['round_size'], "minimum round size")
    print(results['pvalue'], "pvalue for kmax", results['kmax'])

    # SUITE with R2 Bravo
    suite_r2bravo = find_sample_size_for_stopping_prob_efficiently_r2bravo(stop_prob, \
            N_w1, N_l1, N_w2, N_l2, n1, alpha, right = 15000)
    print("\nSUITE with R2 Bravo")
    print(suite_r2bravo['round_size'], "minimum round size")
    print(suite_r2bravo['combined_pvalue'], "pvalue for kmax")

    # update data struct
    data['audits'].append({
        'fractional_margin':fractional_margin,
        'prop_polling':polling_proportion,
        'N_relevant':N_relevant,
        'N_w':N_w,
        'N_l':N_l,
        'N_2':N_2,
        'N_1':N_1,
        'N_w1':N_w1,
        'N_l1':N_l1,
        'N_w2':N_w2,
        'N_l2':N_l2,
        'comparisons':n1,
        'round_size_suite_r2bravo': suite_r2bravo['round_size'],
        'pvalue_suite_r2bravo': suite_r2bravo['combined_pvalue'],
        'round_size_new_method': results['round_size'],
        'pvalue_new_method': results['pvalue'],
        'kmax_new_method': results['kmax']
    })

    # add to json file
    with open(data_file_name, 'w') as outfile:
            json.dump(data, outfile, indent=2)

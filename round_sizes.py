"""
The functions in this file are used to find the minimum round sizes that
produce a desired probability of stopping. Contest-wide Minerva and R2 Bravo 
audits as well as 2-strata audits with either Minerva or R2 Bravo are 
included. 

Additionally there are functions for computing the probablility distribution
over possible pvalues for given round sizes in a 2-strata audit.

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

def compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None, combine_func=None):
    """
    Computes and returns lists of k values, their associated combined pvalue,
    and their probability under the null hypothesis for a 2-strata audit using
    Minerva in the ballot polling stratum. Assumes no errors in the comparisons.

    Args:
        N_w1 (int): reported number of votes for the winner in the comparison stratum
        N_l1 (int): reported number of votes for the loser in the comparison stratum
        N_w2 (int): reported number of votes for the winner in the polling stratum
        N_l2 (int): reported number of votes for the loser in the polling stratum
        n1 (int): number of comparisons
        n2 (int): first round size in the polling stratum
        alpha (float): risk limit
        underlying (dict): feature not yet implemented (coming soon to a repo near you!)
        combine_func (function): Optional function for combining pvalues (default to Fisher's method)

    Return {}:
        possible_winner_votes ([int]): possible number of winner votes in the polling sample
        dist_over_winner_votes ([float]): probability of each possible number of winner votes
        pvalues ([float]): combined pvalue resulting from each possible number of winner votes

    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2
    
    feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    possible_winner_votes = range(0, n2 + 1)
    dist_over_winner_votes = binom.pmf(possible_winner_votes, n2, N_w2 / N_2)

    pvalues = []
    lambdas = []

    for k, pr_k in zip(possible_winner_votes, dist_over_winner_votes):

        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, k, n2 - k, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            minerva_pvalue_direct_count(winner_votes=k, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        results = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N_1, N_w2, N_l2, N_2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                               modulus=mod, alpha=alpha, \
                               feasible_lambda_range=feasible_lambda_range, combine_func=combine_func)
        pvalue = results['max_pvalue']
        refined = results['refined']
        alloc_lambda = results['allocation lambda']
        
        pvalues.append(pvalue)
        lambdas.append(alloc_lambda)
        print("for k="+str(k)+"   pval="+str(round(pvalue,4))+"   refined:"+str(refined))

    return {
        "possible_winner_votes":possible_winner_votes,
        "dist_over_winner_votes":dist_over_winner_votes,
        "pvalues":pvalues,
        "alloc_lambda":lambdas
    }
        

def compute_stopping_probability_whole_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):
    """
    Computes the stopping probability for the given polling stratum first
    round sizes in a 2-strata audit with Minerva.
    Computes the full probability distribution over pvalues to do so. (AKA real slow)

    Note/Plan: Come back and search for kmin then find pr[k >= kmin | alt]
                Should work and be faster...

    Args:
        N_w1 (int): reported number of votes for the winner in the comparison stratum
        N_l1 (int): reported number of votes for the loser in the comparison stratum
        N_w2 (int): reported number of votes for the winner in the polling stratum
        N_l2 (int): reported number of votes for the loser in the polling stratum
        n1 (int): number of comparisons
        n2 (int): first round size in the polling stratum
        alpha (float): risk limit
        underlying (dict): feature not yet implemented (coming soon to a repo near you!)

    Return (float):
        the probability of stopping for the given round sizes
    """

    results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

    possible_winner_votes = results["possible_winner_votes"]
    dist_over_winner_votes = results["dist_over_winner_votes"]
    pvalues = results["pvalues"]

    # find the index of the first pvalue that passes the stopping condition
    index = None
    for i,pvalue in zip(range(0, n2 + 1), pvalues):
        if (pvalue <= alpha):
            index = i
            break

    # if there is not such index then the probability of stopping is 0
    if (index is None):
        return 0

    prob_stop = sum(dist_over_winner_votes[index:])

    return prob_stop

def compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):
    """
    Computes the stopping probability for the given strata sample sizes for 
    a 2-strata audit with Minerva for the polling stratum by finding
    kmin, then finding Pr[k >= kmin | alt].

    NOT YET DONE

    Args:
        N_w1 (int): reported number of votes for the winner in the comparison stratum
        N_l1 (int): reported number of votes for the loser in the comparison stratum
        N_w2 (int): reported number of votes for the winner in the polling stratum
        N_l2 (int): reported number of votes for the loser in the polling stratum
        n1 (int): number of comparisons
        n2 (int): first round size in the polling stratum
        alpha (float): risk limit
        underlying (dict): feature not yet implemented (coming soon to a repo near you!)

    Return (float):
        the probability of stopping for the given round sizes
    """

    results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

    possible_winner_votes = results["possible_winner_votes"]
    dist_over_winner_votes = results["dist_over_winner_votes"]
    pvalues = results["pvalues"]

    # find the index of the first pvalue that passes the stopping condition
    index = None
    for i,pvalue in zip(range(0, n2 + 1), pvalues):
        if (pvalue <= alpha):
            index = i
            break

    # if there is not such index then the probability of stopping is 0
    if (index is None):
        return 0

    prob_stop = sum(dist_over_winner_votes[index:])

    return prob_stop

def find_sample_size_for_stopping_prob_efficiently(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=None, combine_func=None, stouffers=False):
    """
    This function will also compute minimum round size for the 
    passed stopping probability, but it will do so much more 
    efficiently. At each point in the search only one pvalue 
    will be computed. Should have done it this way to begin with.

    Uses Minerva for the ballot polling stratum.
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2

    feasible_lambda_range = calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    left = 1
    if (right is None):
        right = N_2

    print("right: "+str(right))
     
    while(1):
        n2 = math.ceil((left + right) / 2 )

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n2, N_w2 / N_2))

        # compute pvalue for this kmax
        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, kmax, n2 - kmax, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            minerva_pvalue_direct_count(winner_votes=kmax, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        if stouffers is False:
            combination_results = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                                   N_1, N_w2, N_l2, N_2, \
                                   pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                                   alpha=alpha, \
                                   feasible_lambda_range=feasible_lambda_range)

        else:
             combination_results = maximize_stouffers_combined_pvalue(N_w1, N_l1, \
                                   N_1, N_w2, N_l2, N_2, \
                                   pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                                   alpha=alpha, \
                                   stouffers=combine_func)
        pvalue = combination_results['max_pvalue']
        pvalue_comparison = combination_results['pvalue1']
        pvalue_polling = combination_results['pvalue2']
        alloc_lambda = combination_results['allocation lambda']
       
        # update binary search bounds
        if (pvalue > alpha):
            left = n2
        elif (pvalue <= alpha):
            right = n2

        # if left = right then the initial right bound was too small
        if (right == left):
            print("required round size is too larger")
            return None # not sure if returning None is proper thing to do here...

        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1 and n2 == right):
            #print(combination_results['refined'])
            return  {
                        "round_size":right,
                        "combined_pvalue":pvalue,
                        "comparison_pvalue":pvalue_comparison,
                        "polling_pvalue":pvalue_polling,
                        "alloc_lambda":alloc_lambda
                    }


def find_sample_size_for_stopping_prob_minerva(stopping_probability, N_w, N_l, alpha, underlying=None, right=None):
    """
    Finds the first round size that achieves the passed stopping_probability
    for a Minerva audit (with no stratification). 
    """

    N = N_w + N_l 

    left = 1
    if (right is None):
        right = N
     
    while(1):
        n = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n, N_w / N))

        # compute pvalue for this kmax
        pvalue = minerva_pvalue_direct_count(winner_votes=kmax, n=n, popsize=N, alpha=alpha, Vw=N_w, Vl=N_l, null_margin=0)

        # update binary search bounds
        if (pvalue > alpha):
            left = n
        elif (pvalue <= alpha):
            right = n
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1):
            if (right == N):
                print("required round size is greater than stratum size")
            return right

def find_sample_size_for_stopping_prob_r2bravo(stopping_probability, N_w, N_l, alpha, underlying=None, right=None):
    """
    Finds the first round size that achieves the passed stopping_probability
    for an R2 Bravo audit (with no stratification). 
    """

    N = N_w + N_l 

    left = 1
    right = N
     
    while(1):
        n = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n, N_w / N))

        # compute pvalue for this kmax
        pvalue = r2bravo_pvalue_direct_count(winner_votes=kmax, n=n, popsize=N, alpha=alpha, Vw=N_w, Vl=N_l, null_margin=0)

        # update binary search bounds
        if (pvalue > alpha):
            left = n
        elif (pvalue <= alpha):
            right = n
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1):
            if (right == N):
                print("required round size is greater than stratum size")
            return right


def find_sample_size_for_stopping_prob_efficiently_r2bravo(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None, right=None, combine_func=None, stouffers=False):
    """
    This function will also compute minimum round size for the 
    passed stopping probability, but it will do so much more 
    efficiently. At each point in the search only one pvalue 
    will be computed. Should have done it this way to begin with.
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2

    feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    left = 1
    right = N_2
     
    while(1):
        n2 = math.ceil((left + right) / 2)

        # compute the 1 - stopping_probability quantile of the alt dist
        # kmax where pr[k >= kmax | alt] = stopping_probability
        # floor because we need to ensure at least a stopping_probability prob of stopping
        kmax = math.floor(binom.ppf(1 - stopping_probability, n2, N_w2 / N_2))

        # compute pvalue for this kmax
        cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                       o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N_1,
                                   null_lambda=alloc)

        mod = create_modulus(n1, n2, kmax, n2 - kmax, N_1, margin, 1.03905)

        nocvr_pvalue = lambda alloc: \
            r2bravo_pvalue_direct_count(winner_votes=kmax, n=n2, popsize=N_2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)

        if stouffers is False:
            combination_results = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                                   N_1, N_w2, N_l2, N_2, \
                                   pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                                   modulus=mod, alpha=alpha, \
                                   feasible_lambda_range=feasible_lambda_range, combine_func=combine_func)
        else:
            combination_results = maximize_stouffers_combined_pvalue(N_w1, N_l1, \
                                   N_1, N_w2, N_l2, N_2, \
                                   pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                                   alpha=alpha, \
                                   stouffers=combine_func)

        pvalue = combination_results['max_pvalue']
        pvalue_comparison = combination_results['pvalue1']
        print(pvalue_comparison)
        pvalue_polling = combination_results['pvalue2']
        alloc_lambda = combination_results['allocation lambda']

        
        # update binary search bounds
        if (pvalue > alpha):
            left = n2
        elif (pvalue <= alpha):
            right = n2
 
        # when and right converge, right is the minimum round size that achieves stopping_probability
        if (left == right - 1 and n2 == right):
            if (right == N_2):
                print("requried round size is greater than stratum size")
            return  {
                        "round_size":right,
                        "combined_pvalue":pvalue,
                        "comparison_pvalue":pvalue_comparison,
                        "polling_pvalue":pvalue_polling,
                        "alloc_lambda":alloc_lambda
                    }




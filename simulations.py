"""
Functions in this file perform simulations of 2-strata audits. Since my
work has mostly transitioned to analytical methods at the moment, rather 
than simulations, these functions are not as capable as I plan to make them.

Future work will build these functions out to handle various true underlying
tallies and will allow them to be used to confirm analytical methods as well
as observe how closely the claimed risk of 2-strata audits predict the true
risk.

Additionally, this file contains functions for computing first round pvalues
for Minerva and R2 Bravo. Multiple versions of these functions exist that 
either use r2b2 code (thanks Grant!) or my own code, and each can handle 
different input types (full, ordered sample of 1's and 0's or raw winner
tally).

Oliver Broadrick 2020
"""


import time
import numpy as np
import scipy as sp
import scipy.stats
import scipy.optimize
from ballot_comparison import ballot_comparison_pvalue
from hypergeometric import trihypergeometric_optim
from sprt import ballot_polling_sprt
import matplotlib.pyplot as plt
import numpy.testing
from contest import ContestType
from contest import Contest
from minerva_s import Minerva_S
from fishers_combination import create_modulus, maximize_fisher_combined_pvalue, calculate_lambda_range
from scipy.stats import binom
import math
import matplotlib.pyplot as plt

def simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha,
    reps=10000, verbose=False, feasible_lambda_range=None, underlying=None):
    """
    Simulate the Fisher method of combining a ballot comparison audit
    and ballot polling minerva audit, assuming the true results contain
    underlying winner votes. 
    Return the fraction of simulations where the the audit successfully
    confirmed the election results for each of several audits.
    
    Parameters
    ----------
    N_w1 : int
        votes for the reported winner in the ballot comparison stratum
    N_l1 : int
        votes for the reported loser in the ballot comparison stratum
    N_w2 : int
        votes for the reported winner in the ballot polling stratum
    N_l2 : int
        votes for the reported loser in the ballot polling stratum
    n1 : int
        sample size in the ballot comparison stratum
    n2 : int
        sample size in the ballot polling stratum
    alpha : float
        risk limit
    reps : int
        number of times to simulate the audit. Default 10,000
    verbose : bool
        Optional, print iteration number if True
    feasible_lambda_range : array-like
        lower and upper limits to search over lambda. Optional, but will speed up the search
    underlying : int
        true count of votes for winner overall (default assumes alt)
    
    Returns
    -------
    dict : fractions of simulations where the the audits successfully
    confirmed the election results
    """
    if underlying is None:
        underlying = N_w2

    N1 = N_w1 + N_l1
    N2 = N_w2 + N_l2
    margin = (N_w1 + N_w2) - (N_l1 + N_l2)

    # Population generated based on 'underlying' (assumed winner count)
    pop2 = [1]*underlying + [0]*(N2 - underlying)
    
    cvr_pvalue = lambda alloc: ballot_comparison_pvalue(n=n1, gamma=1.03905, \
                                   o1=0, u1=0, o2=0, u2=0,
                                   reported_margin=margin, N=N1,
                                   null_lambda=alloc)

    fisher_pvalues_r2_bravo = np.zeros(reps)
    fisher_pvalues_r2_bravo_direct = np.zeros(reps)
    fisher_pvalues_minerva = np.zeros(reps)
    fisher_pvalues_minerva_direct = np.zeros(reps)
   
    # Generate samples
    samples = []
    for i in range(reps):
        sam = np.random.choice(pop2, n2, replace=True)
        samples.append(sam)
    """
    # R2 BRAVO
    start = time.time()
    for i, sam in zip(range(len(samples)),samples):
        nw2 = np.sum(sam == 1)
        nl2 = np.sum(sam == 0)
        mod = create_modulus(n1, n2, nw2, nl2, N1, margin, 1.03905)
        nocvr_pvalue_r2_bravo = lambda alloc: \
            ballot_polling_sprt(sample=sam, popsize=N2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)['pvalue']
        fisher_pvalues_r2_bravo[i] = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N1, N_w2, N_l2, N2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue_r2_bravo], \
                               modulus=mod, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']
    r2_bravo_time = time.time() - start
    """

    # R2 BRAVO (direct)
    start = time.time()
    for i, sam in zip(range(len(samples)),samples):
        nw2 = np.sum(sam == 1)
        nl2 = np.sum(sam == 0)
        mod = create_modulus(n1, n2, nw2, nl2, N1, margin, 1.03905)
        nocvr_pvalue_r2_bravo_direct = lambda alloc: \
            r2_bravo_pvalue_direct(sample=sam, popsize=N2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)
        fisher_pvalues_r2_bravo_direct[i] = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N1, N_w2, N_l2, N2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue_r2_bravo_direct], \
                               modulus=mod, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']
    r2_bravo_direct_time = time.time() - start
    
    """
    # Minerva (Grant/r2b2)
    start = time.time()
    for i, sam in zip(range(len(samples)),samples):
        nw2 = np.sum(sam == 1)
        nl2 = np.sum(sam == 0)
        mod = create_modulus(n1, n2, nw2, nl2, N1, margin, 1.03905)
        nocvr_pvalue_minerva = lambda alloc: \
            minerva_pvalue(sample=sam, popsize=N2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)
        fisher_pvalues_minerva[i] = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N1, N_w2, N_l2, N2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue_minerva], \
                               modulus=mod, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']
    minerva_time = time.time() - start
    """ 
    """
    # Minerva (direct)
    start = time.time()
    for i, sam in zip(range(len(samples)),samples):
        nw2 = np.sum(sam == 1)
        nl2 = np.sum(sam == 0)
        mod = create_modulus(n1, n2, nw2, nl2, N1, margin, 1.03905)
        nocvr_pvalue_minerva_direct = lambda alloc: \
            minerva_pvalue_direct(sample=sam, popsize=N2, alpha=alpha, \
                                Vw=N_w2, Vl=N_l2, \
                                null_margin=(N_w2-N_l2) - alloc*margin)
        fisher_pvalues_minerva_direct[i] = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N1, N_w2, N_l2, N2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue_minerva_direct], \
                               modulus=mod, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']
    minerva_direct_time = time.time() - start
    """
    """
    "r2_bravo" : np.mean(fisher_pvalues_r2_bravo <= alpha),
    "r2_bravo_time" : r2_bravo_time,
    "r2_bravo_avg_pval" : np.mean(fisher_pvalues_r2_bravo),

    "r2_bravo_direct" : np.mean(fisher_pvalues_r2_bravo_direct <= alpha),
    "r2_bravo_direct_time" : r2_bravo_direct_time,
    "r2_bravo_direct_avg_pval" : np.mean(fisher_pvalues_r2_bravo_direct),

    "minerva" : np.mean(fisher_pvalues_minerva <= alpha),
    "minerva_time" : minerva_time,
    "minerva_avg_pval" : np.mean(fisher_pvalues_minerva),

    "minerva_direct" : np.mean(fisher_pvalues_minerva_direct <= alpha),
    "minerva_direct_time" : minerva_direct_time,
    "minerva_direct_avg_pval" : np.mean(fisher_pvalues_minerva_direct)
    """
    return  {
    "r2_bravo_direct" : np.mean(fisher_pvalues_r2_bravo_direct <= alpha),
    "r2_bravo_direct_time" : r2_bravo_direct_time,
    "r2_bravo_direct_avg_pval" : np.mean(fisher_pvalues_r2_bravo_direct),
            }

def minerva_pvalue(sample, popsize, alpha, Vw, Vl, null_margin):
    """Computes the pvalue for a one-round minerva audit with the passed values.
    Uses an adapted version of Grant's Minerva code in r2b2 (adapted for null margins).

    Parameters:
        sample : list of 1's (vote for winner) and 0's (vote for loser)
        popsize : total ballots in stratum
        alpha : risk limit
        Vw : reported votes for winner in stratum
        Vl : reported votes for loser in stratum
        null_margin : the margin in votes assumed under the null

    Returns:
        float : the minerva pvalue
    """
    contest = Contest(popsize,{'A':Vw,'B':Vl},1,['A'],ContestType.PLURALITY)
    audit = Minerva_S(alpha, 1.0, contest, null_margin)

    n = len(sample)

    audit.rounds.append(n)
    audit.current_dist_reported()
    audit.current_dist_null()

    k = np.sum(sample == 1)
    x = (popsize + null_margin) / 2

    if x < k or popsize - x < np.sum(sample == 0):
        return 0

    pvalue = audit.stopping_condition(k)['pvalue']

    return min(pvalue)

def minerva_pvalue_direct(sample, popsize, alpha, Vw, Vl, null_margin):
    """Computes the pvalue for a one-round minerva audit with the passed values.
    Makes computations directly (rather than with Grant's r2b2 Minerva code).

    Parameters:
        sample : list of 1's (vote for winner) and 0's (vote for loser)
        popsize : total ballots in stratum
        alpha : risk limit
        Vw : reported votes for winner in stratum
        Vl : reported votes for loser in stratum
        null_margin : the margin in votes assumed under the null

    Returns:
        float : the minerva pvalue
    """
    n = len(sample)
    k = np.sum(sample == 1)
    x = (popsize + null_margin) / 2

    if x < k or popsize - x < np.sum(sample == 0):
        return 0

    alt_dist = binom.pmf(range(0, n + 1), n, Vw / popsize)
    null_dist = binom.pmf(range(0, n + 1), n, x / popsize)

    alt_tail = sum(alt_dist[k:])
    null_tail = sum(null_dist[k:])

    pvalue = null_tail / alt_tail

    return pvalue

def minerva_pvalue_direct_count(winner_votes, n, popsize, alpha, Vw, Vl, null_margin):
    """Computes the pvalue for a one-round minerva audit with the passed values.
    Makes computations directly (rather than with Grant's r2b2 Minerva code).
    Uses the count of winner votes rather than the sample structure that SUITE uses.

    Parameters:
        winner_votes : number of votes for the winner in the sample (rest for loser)
        popsize : total ballots in stratum
        alpha : risk limit
        Vw : reported votes for winner in stratum
        Vl : reported votes for loser in stratum
        null_margin : the margin in votes assumed under the null

    Returns:
        float : the minerva pvalue
    """
    k = winner_votes
    x = (popsize + null_margin) / 2

    if x < k or popsize - x < n - k:
        return 0

    alt_dist = binom.pmf(range(0, n + 1), n, Vw / popsize)
    null_dist = binom.pmf(range(0, n + 1), n, x / popsize)

    alt_tail = sum(alt_dist[k:])
    null_tail = sum(null_dist[k:])

    pvalue = null_tail / alt_tail

    return pvalue


def r2_bravo_pvalue_direct(sample, popsize, alpha, Vw, Vl, null_margin):
    """Computes the pvalue for a one-round R2 BRAVO audit with the passed values.

    Parameters:
        sample : list of 1's (vote for winner) and 0's (vote for loser)
        popsize : total ballots in stratum
        alpha : risk limit
        Vw : reported votes for winner in stratum
        Vl : reported votes for loser in stratum
        null_margin : the margin in votes assumed under the null

    Returns:
        float : the R2 BRAVO pvalue
    """
    n = len(sample)
    k = np.sum(sample == 1)
    x = (popsize + null_margin) / 2

    if x < k or popsize - x < np.sum(sample == 0):
        return 0

    alt = binom.pmf(k, n, Vw / popsize)
    null = binom.pmf(k, n, x / popsize)

    pvalue = null / alt

    return pvalue

def estimate_round_size_for_stopping_prob(prob, N_w1, N_l1, N_w2, N_l2, alpha, underlying=None):
    """
    Estimate the round size required to produce the desired stopping probability.

    Parameters
    ----------
    N_w1 : int
        votes for the reported winner in the ballot comparison stratum
    N_l1 : int
        votes for the reported loser in the ballot comparison stratum
    N_w2 : int
        votes for the reported winner in the ballot polling stratum
    N_l2 : int
        votes for the reported loser in the ballot polling stratum
    alpha : float
        risk limit
    underlying : int
        true count of votes for winner overall (default assumes alt)
    
    Returns
    -------
    dict : fractions of simulations where the the audits successfully
    confirmed the election results
    """
    N1 = N_w1 + N_l1
    N2 = N_w2 + N_l2
    N = N1 + N2

    # Binary search for sample size
    left = 1
    right = N / 8 # more than third seems excessive

    tol = .01

    while (1):
        mid = (right + left) / 2
        
        # For now, just assume a 50-50 allocation between the strata (later can also search to maximize stopping prob for this allocation)
        n1 = math.ceil(mid / 2)
        n2 = math.floor(mid / 2)

        # Simulate for this round size (and allocation) to obtain stopping probality (later can obtain and compare analytically)
        reps = 50
        results = simulate_fisher_combined_audits(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, reps=reps, feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2))
        pr_stop = results['minerva']

        #print("round size: "+str(mid)+"   pr_stop: "+str(pr_stop))

        if pr_stop < prob:
            left = mid
        elif pr_stop < prob + tol:
            return mid
        else:
            right = mid

def compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):
    """
    Hopefully compute a distribution over possible pvalues
    """

    N_1 = N_w1 + N_l1
    N_2 = N_w2 + N_l2
    margin = N_w1 + N_w2 - N_l1 - N_l2
    
    feasible_lambda_range=calculate_lambda_range(N_w1, N_l1, N_1, N_w2, N_l2, N_2)

    possible_winner_votes = range(0, n2 + 1)
    dist_over_winner_votes = binom.pmf(possible_winner_votes, n2, N_w2 / N_2)

    pvalues = []

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

        pvalue = maximize_fisher_combined_pvalue(N_w1, N_l1, \
                               N_1, N_w2, N_l2, N_2, \
                               pvalue_funs=[cvr_pvalue, nocvr_pvalue], \
                               modulus=mod, \
                               feasible_lambda_range=feasible_lambda_range)['max_pvalue']

        pvalues.append(pvalue)
        #print("for k="+str(k)+"   pval="+str(pvalue))

    return {
        "possible_winner_votes":possible_winner_votes,
        "dist_over_winner_votes":dist_over_winner_votes,
        "pvalues":pvalues
    }
        

def compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None):

    results = compute_dist_over_pvalues(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

    possible_winner_votes = results["possible_winner_votes"]
    dist_over_winner_votes = results["dist_over_winner_votes"]
    pvalues = results["pvalues"]


    index = -1
    # find the index of the first pvalue less than the risk limit
    for i,pvalue in zip(range(0, n2 + 1), pvalues):
        if (pvalue <= alpha):
            index = i
            break
    prob_stop = sum(dist_over_winner_votes[index:])

    return prob_stop


def find_sample_size_for_stopping_prob(stopping_probability, N_w1, N_l1, N_w2, N_l2, n1, alpha, underlying=None):
    """
    Hopefully will find the minimum sample size for the ballot polling stratum
    which will achieve the passed stopping_probability.
    """
    start = time.time()

    N_2 = N_w2 + N_l2
     
    left = 1
    right = N_2 / 8
    while(1):
        # current value to test
        n2 = math.ceil((left + right) / 2)

        # analytically compute the stopping probability
        stop = compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, n2, alpha, underlying=None)

        # print results for fun
        #print("n2: "+str(n2)+"  pr_stop: "+str(stop))

        # update binary search bounds
        if (stop < stopping_probability):
            left = n2
        elif (stop > stopping_probability):
            right = n2
        
        # the left and right bounds should slowly converge,
        # eventually leaving the left bound one less than the right
        # in this case the right bound is the desired value
        if (left == right - 1):
            right_pr_stop = compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, right, alpha, underlying=None)
            print("n2: "+str(right)+"  pr_stop: "+str(right_pr_stop)+"  took: "+str((time.time()-start)/60)+" minutes")
            left_pr_stop = compute_stopping_probability(N_w1, N_l1, N_w2, N_l2, n1, left, alpha, underlying=None)
            print("one lower for confirmation: n2: "+str(left)+"  pr_stop: "+str(left_pr_stop))
            return {
                'round_size':right,
                'stopping_prob':right_pr_stop,
                'one_lower':left,
                'one_lower_prob':left_pr_stop
            }

def r2bravo_pvalue_direct_count(winner_votes, n, popsize, alpha, Vw, Vl, null_margin):
    """Computes the pvalue for a one-round r2bravo audit with the passed values.
    Uses the count of winner votes rather than the sample structure that SUITE uses.

    Parameters:
        sample : list of 1's (vote for winner) and 0's (vote for loser)
        popsize : total ballots in stratum
        alpha : risk limit
        Vw : reported votes for winner in stratum
        Vl : reported votes for loser in stratum
        null_margin : the margin in votes assumed under the null

    Returns:
        float : the minerva pvalue
    """
    k = winner_votes
    x = (popsize + null_margin) / 2

    if x < k or popsize - x < n - k:
        return 0

    alt_pr = binom.pmf(k, n, Vw / popsize)
    null_pr = binom.pmf(k, n, x / popsize)

    pvalue = null_pr / alt_pr

    return min([pvalue,1])














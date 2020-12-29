"""
In this file I'll write some functions for using the stratified method
specifically for 2-strata audits in which both strata are polling.
Hopefully, I'll quickly get to the point where I can find first round
sizes for both strata that are the same and achieve a given 
stopping probability.

Oliver Broadrick 2020
"""

import math
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from matplotlib import cm
import time

def generate_polling_dists(n, Vw, Vl, null_margin=0, plot=False):
    """
    Generates first round distributions over winner votes for 
    ballot polling audit.

    Parameters:
        k : number of votes for the winner in the sample (rest for loser)
        N : total ballots
        Vw : reported votes for winner
        Vl : reported votes for loser
        null_margin : the margin in votes assumed under the null
                Optional: default to 0 (for stratified application)
        plot : Optional: if true, will plot the distributions (default false)

    Returns: dict with
        dist_range : values of k to which the distributions correspond
        alt_dist : probability distribution over winner votes under alt
        null_dist : probability distribution over winner votes under null
    """
    # relevant ballots
    N = Vw + Vl

    # winner votes under the null
    x = (N + null_margin) / 2

    # don't bother with bayesian prior other than .5 and .5 at x and Vw
    # could easily extend to accept other priors, same as in comparison stratum
    dist_range = range(0, n + 1)
    alt_dist = binom.pmf(dist_range, n, Vw / N)
    null_dist = binom.pmf(dist_range, n, x / N)

    if plot:
        # plot for viewing pleasure
        plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
        plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
        plt.legend()
        plt.title("Polling Stratum Distributions over Winner Votes")
        plt.show()

    return {
        'dist_range': dist_range,
        'alt_dist': alt_dist,
        'null_dist': null_dist
    }

def generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, null_margin1=0, \
                        plot=False, print_data=True):
    """
    Generates the joint probability distribution for the first round
    of a 2-strata (comparison and polling) audit.

    Parameters:
        N_w1 : winner votes in comparison stratum
        N_l1 : loser votes in comparison stratum
        N_w2 : winner votes in polling stratum
        N_lw : loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        null_margin1 : the margin in votes assumed under the null in the comparison stratum
                        (defaults to zero)
        plot : Optional: if true, will plot the distributions (default false)

    Returns:
        alt_joint : joint probability distribution under the alternative hypothesis
        null_joint : joint probability distribution under the null hypothesis
    """

    # generate comparison dists
    polling_dists1 = generate_polling_dists(n1, N_w1, N_l1, null_margin1)
    dist_range_polling1 = polling_dists1['dist_range']
    alt_dist_polling1 = polling_dists1['alt_dist']
    null_dist_polling1 = polling_dists1['null_dist']

    # generate polling dists
    polling_dists2 = generate_polling_dists(n2, N_w2, N_l2, -1 * null_margin1)
    dist_range_polling2 = polling_dists2['dist_range']
    alt_dist_polling2 = polling_dists2['alt_dist']
    null_dist_polling2 = polling_dists2['null_dist']

    # compute joint dists
    alt_joint = np.empty((len(dist_range_polling1), len(dist_range_polling2)))
    null_joint = np.empty((len(dist_range_polling1), len(dist_range_polling2)))
    for k_p1 in dist_range_polling1:
        for k_p2 in dist_range_polling2:
            alt_joint[k_p1][k_p2] = alt_dist_polling1[k_p1] * alt_dist_polling[k_p2]
            null_joint[k_p1][k_p2] = null_dist_polling1[k_p1] * null_dist_polling[k_p2]

    return {
        'alt_joint': alt_joint,
        'null_joint': null_joint
    }

def compute_pvalue_from_joint_dist(k_c, k_p, alt_joint, null_joint):
    """
    Compute the pvalue for the passed joint distribution.

    Parameters:
        k_c: matches in comparison sample
        k_p: winner votes in polling sample
        alt_joint: joint
        alt_joint: the joint probability distribution under the alternative hypothesis
        null_joint: the joint probability distribution under the null hypothesis

    Returns:
        float: pvalue (ratio of corners)
    """
    # sum the 'corners' of the joint dists
    alt_corner = 0
    for k in range(k_c, len(alt_joint)):
        alt_corner += alt_joint[k][k_p:].sum()
    null_corner = 0
    for k in range(k_c, len(null_joint)):
        null_corner += null_joint[k][k_p:].sum()

    # pvalue is ratio of the 'corners'
    return null_corner / alt_corner
    
def compute_pvalue(k_c, k_p, alt_c, alt_p, null_c, null_p):
    """
    Compute the pvalue for the passed sample and distributions.
    Directly calculate each joint probability from the 
    independent distributions passed.

    Parameters:
        k_c: matches in comparison sample
        k_p: winner votes in polling sample
        alt_c : alternative distribution for comparison stratum
        alt_p : alternative distribution for comparison stratum
        null_c : null distribution for comparison stratum
        null_p : null distribution for comparison stratum

    Returns:
        float: pvalue (ratio of corners)
    """
    alt_corner = 0
    null_corner = 0

    # compute sum of the 'corner' of the joint dist
    for i in range(k_c, len(alt_c)):
        for j in range (k_p, len(alt_p)):
            alt_corner += alt_c[i] * alt_p[j]
            null_corner += null_c[i] * null_p[j]

    # pvalue is ratio of the 'corners'
    return null_corner / alt_corner

def find_kmin_pairs(alpha, alt_joint, null_joint):
    """
    Finds all valid kmin pairs for the passed joint distributions 
    and risk limit.

    Parameters:
        alpha: risk limit
        alt_joint: the joint probability distribution under the alternative hypothesis
        null_joint: the joint probability distribution under the null hypothesis

    Returns: dict with
        k_mins_c: values of kmin for comparison stratum
        k_mins_p: corresponding values of kmin for polling stratum
    """

    k_mins_c = []
    k_mins_p = []
    for k_c in range(len(alt_joint)):
        for k_p in range(len(alt_joint[0])):

            # test if this pair of k's meet the stopping condition
            pvalue = compute_pvalue_from_joint_dist(k_c, k_p, alt_joint, null_joint)

            if not math.isnan(pvalue) and pvalue <= alpha:
                # add these k mins to the lists
                k_mins_c.append(k_c)
                k_mins_p.append(k_p)

                """
                print("k_c:",k_c," k_p:",k_p)
                print("pvalue:",pvalue)
                """

                # break after k_p min is found for this value of k_c
                break

    return {
        'k_mins_c': k_mins_c,
        'k_mins_p': k_mins_p
    }

def plot_joint_dist(alt_joint, null_joint, k_mins_c=None, k_mins_p=None):
    """
    For viewing pleasure, plots the joint probability distribution as 
    well as the kmin line if kmins are provided.

    Parameters:
        alt_joint : joint probability distribution under the alternative hypothesis
        null_joint : joint probability distribution under the null hypothesis
        k_mins_c: values of kmin for comparison stratum
        k_mins_p: corresponding values of kmin for polling stratum
    """
    fig = plt.figure()
    axes = fig.gca(projection='3d')

    dist_range_polling = range(len(alt_joint[0] + 1))
    dist_range_comparison = range(len(alt_joint + 1))
    X, Y = np.meshgrid(dist_range_polling, dist_range_comparison)

    surf1 = axes.plot_surface(X, Y, alt_joint, label='joint dist under alt', cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # matplotlib error fixed with these two lines
    surf1._facecolors2d=surf1._facecolors3d
    surf1._edgecolors2d=surf1._edgecolors3d

    surf2 =axes.plot_surface(X, Y, null_joint, label='joint dist under null', cmap=cm.magma, linewidth=0, antialiased=False)
    # matplotlib error fixed with these two lines
    surf2._facecolors2d=surf2._facecolors3d
    surf2._edgecolors2d=surf2._edgecolors3d

    if k_min_c is not None and k_min_p is not None:
        for k_min_c, k_min_p in zip(k_mins_c, k_mins_p):
            z = [0,.06]
            axes.plot([k_min_c,k_min_c], [k_min_p,k_min_p], z, linewidth=3, color='g')

    axes.legend(fontsize=20)

    axes.set_xlabel('Polling Stratum: Winner Ballots', fontsize=20, labelpad=24)
    axes.set_ylabel('Comparison Stratum: Matching Ballots', fontsize=20, labelpad=24)
    axes.set_zlabel('Probability', fontsize=20, labelpad=24)

    plt.setp(axes.get_xticklabels(), fontsize=18)
    plt.setp(axes.get_yticklabels(), fontsize=18)
    plt.setp(axes.get_zticklabels(), fontsize=18)

    plt.show()

def compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2, k_p=0):
    """
    Computes upper and lower bounds for the number of winner votes in 
    the comparison stratum, x1. 

    Parameters:
        N_w1 : winner ballots in comparison stratum
        N_l1 : loser ballots in comparison stratum
        N_w2 : winner ballots in polling stratum
        N_lw : loser ballots in polling stratum
        k_p : winner ballots already drawn
            Optional: defaults to zero

    Return: dict with
        x1_l : lower bound on x1
        x1_u : upper bound on x1
    """
    
    N1 = N_w1 + N_l1
    N2 = N_w2 + N_l2
    N = N1 + N2

    winner_votes_null = math.floor(N / 2)

    x1_l = max(0, winner_votes_null - N2)
    x1_u = min(N1 - k_p, winner_votes_null)

    return {
        'x1_l' : x1_l,
        'x1_u' : x1_u
    }

def maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, \
        lower_bound=None, upper_bound=None, plot=False, print_data=False):
    """
    Maximizes the joint pvalue for the given sample by searching
    over all possible allocations of winner votes under the null
    hypothesis. If maximum pvalue is greater than 1, 1 is returned.

    Parameters:
        k_c : matches in comparison sample
        k_p : winner votes in polling sample
        N_w1 : reported winner votes in comparison stratum
        N_l1 : reported loser votes in comparison stratum
        N_w2 : reported winner votes in polling stratum
        N_lw : reported loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        lower_bound : lower bound for winner votes in comparison stratum under the null
        upper_bound : upper bound for winner votes in comparison stratum under the null
        plot : optional: set true for plot of pvalues vs. winner vote allocations

    Return: dict with
        pvalue : maximum pvalue found
    """
    
    if lower_bound is None or upper_bound is None:
        # compute bounds for search
        bounds = compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2, k_p)
        lower_bound = bounds['x1_l']
        upper_bound = bounds['x1_u']

    # define a step size and generate lists for testing
    # each test x is a number of winner votes in the comparison stratum
    divisions = 10 # could tweak this for effiency
    step_size = math.ceil((upper_bound - lower_bound) / divisions)
    test_xs = np.arange(lower_bound, upper_bound, step_size)
    pvalues = np.empty_like(test_xs, dtype=float)

    for i in range(len(test_xs)):
        # compute null margin based on winner votes
        x1 = test_xs[i]
        null_margin1 = x1 - (N_w1 + N_l1 - x1)

        """
        # generate joint distributions
        dists = generate_joint_dist(N_w1, N_l1, N_w2, N_l2, n1, n2, null_margin1)
        alt_joint = dists['alt_joint']
        null_joint = dists['null_joint']
        """
        # generate comparison dists
        polling1_dists = generate_polling_dists(n1, N_w1, N_l1, null_margin1, plot=False)
        alt_c = polling1_dists['alt_dist']
        null_c = polling1_dists['null_dist']

        # generate polling dists
        polling2_dists = generate_polling_dists(n2, N_w2, N_l2, -1 * null_margin1, plot=False)
        alt_p = polling2_dists['alt_dist']
        null_p = polling2_dists['null_dist']

        # compute pvalue
        pvalue = compute_pvalue(k_c, k_p, alt_c, alt_p, null_c, null_p)
        pvalues[i] = pvalue
        
        # if pvalue greater than 1, just return 1
        if pvalue > 1:
            return {
                        'pvalue' : 1,
                        'x1' : x1,
                        'search_iterations' : 1
                    }

    # get the maximum pvalue found
    max_index = np.argmax(pvalues)
    max_pvalue = pvalues[max_index]
    max_x = test_xs[max_index]

    if plot:
        # plot for viewing pleasure
        plt.plot(test_xs, pvalues, label='pvals for testxs', marker='o', color='b', linestyle='solid')
        plt.show()

    """
    print("max_index:",max_index)
    print("max_pvalue:",max_pvalue)
    print("max_x:",max_x)
    """
    # when step_size has reached 1, search is over
    if step_size == 1:
        return {
            'pvalue' : max_pvalue,
            'x1' : max_x,
            'search_iterations' : 1
        }
    else:
        # set bounds for refined search
        lower = max_x - step_size
        upper = max_x + step_size

        # perform refined search
        refined = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=lower, upper_bound=upper)

        # increase iterations by 1 and return results
        refined['search_iterations'] += 1
        return refined

def find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, stop_prob, alpha):
    """
    Finds the minimum polling first round size that achieves the 
    desired probability of stopping, stop_prob, assuming no errors
    in the comparison sample, by linear search starting at 1.

    Parameters:
        N_w1 : reported winner votes in comparison stratum
        N_l1 : reported loser votes in comparison stratum
        N_w2 : reported winner votes in polling stratum
        N_lw : reported loser votes in polling stratum
        n1 : comparison stratum sample size
        stop_prob : desired stopping probability
        alpha : risk limit

    Returns:
        int : minimum first polling round size
    """

    n = 1
    while(1):
        # find kmax such that Pr[k >= kmax] = stop_prob
        kmax = math.floor(binom.ppf(1 - math.sqrt(stop_prob), n, N_w2 / (N_w2 + N_l2)))

        # get pvalue for kmax
        pvalue = maximize_joint_pvalue(kmax, kmax, N_w1, N_l1, N_w2, N_l2, n, n, plot=False)['pvalue']
        #print(pvalue)

        # print n and pvalue for viewing pleasure
        #print("n:",n,"pvalue:",pvalue,"kmax:",kmax)

        # return n when pvalue for kmax meets stopping condition
        if (pvalue < alpha):
            return {
                "round_size": n,
                "pvalue": pvalue,
                "kmax" : kmax
            }

        # increment round size
        n += 1


### FOR NOW I'LL JUST START THE SCRIPT HERE...

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
stop_prob = .9
alpha = .1

print("new method both polling:")
results = find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, stop_prob, alpha)
print(results)
print("total round size (both strata):",results['round_size']*2)

### OK SO THAT WORKS, BUT ALSO WHAT WOULD BE NORMAL R2 BRAVO OF THE WHOLE THING...?

def minerva_pval(k, n, N_w, N_l):
    num = 0
    denom = 0
    for i in range(k,n+1):
        num += binom.pmf(i, n, .5)
        denom += binom.pmf(i, n, N_w/(N_w+N_l))
    return num / denom

def find_minimum_round_size_r2bravo(N_w, N_l, stop_prob, alpha):
    n = 1
    while(1):
        # find kmax such that Pr[k >= kmax] = stop_prob
        kmax = math.floor(binom.ppf(1 - stop_prob, n, N_w / (N_w+N_l)))

        # get pvalue for kmax
        pvalue = minerva_pval(kmax, n, N_w, N_l)

        # print n and pvalue for viewing pleasure
        #print("n:",n,"pvalue:",pvalue,"kmax:",kmax)

        # return n when pvalue for kmax meets stopping condition
        if (pvalue < alpha):
            return {
                "round_size": n,
                "pvalue": pvalue,
                "kmax" : kmax
            }

        # increment round size
        n += 1

print("minerva of whole contest:")
print(find_minimum_round_size_r2bravo(N_w1+N_w2, N_l1+N_l2, stop_prob, alpha))


"""
OK FINALLY GETTING THE EXPECTED RESULT THAT THIS NEW METHOD IS NOT IN FACT
BETTER THAN JUST DOING A MINERVA AUDIT OF THE WHOLE CONTEST IN ONE GO...
OF COURSE AS I WAS WHIPPING THIS SCRIPT TOGETHER AT FIRST I WAS GETTING
THAT RESULT AND FINALLY FIXING THE BUGS AND GETTING THE EXPECTED RESULT
PRINTED OUT WAS THE BIGGEST RELEASE OF STRESS EVER... BIG HAPPY SIGH.. PHEW
"""

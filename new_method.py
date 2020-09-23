"""
Functions for a new method of stratified audit with 2 strata (ballot
comparison and ballot polling) that finds a sharper p-value by 
analytically finding the overall pvalue, rather than combining
independent p-values.

Oliver Broadrick 2020
"""

import math
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt
from matplotlib import cm
import time

def generate_comparison_dists(n, Vw, Vl, null_margin=0, \
                                plot=False, print_data=False):
    """
    Generates the first round probability distributions over number 
    of matches for both the null and alternative hypothesis of 
    a comparison audit.

    For simplicity (this is just proof of concept) I'm considering 
    all overstatement errors 2-vote overstatements

    Parameters:
        n : sample size
        N : total ballots
        Vw : reported votes for winner
        Vl : reported votes for loser
        null_margin : the margin in votes under the null
                    Optional: defaults to 0 (used for stratified application)
        plot : Optional: if true, will plot the distributions (default false)

    Returns: dict with
        dist_range : values of k to which the distributions correspond
        alt_dist : probability distribution over matches under the alternative hypothesis
        null_dist : probability distribution over matches under the null hypothesis
    """
    # useful bits
    N = Vw + Vl
    reported_margin = Vw - Vl
    
    # compute minimum mismatches under the null
    mismatches_null = math.ceil((reported_margin - null_margin) / 2)
    if mismatches_null < 0:
        # null_margin > reported margin, no overstatements needed
        mismatches_null = 0
    #print("mismatches_null:",mismatches_null)

    # maximum matches under the null
    matches_null = N - mismatches_null
    #print("matches_null:", matches_null)

    # if print_data is true, print xnull and bnull
    if print_data is True:
        # get winner votes under null from null margin
        x1 = (N + null_margin) / 2
        # if in between rougly 270 - 80
        if x1 > 272 and x1 < 283:
            print("x1:",x1)
            print("matches_null:",matches_null)
            print("mismatches_null:",mismatches_null)

    # risk limiting prior (0's, .5 at matches_null, then zeros, then uniform)
    if matches_null >= N:
        prior = np.concatenate([np.zeros(matches_null), np.array([1])])
    else:
        prior = np.concatenate([np.zeros(matches_null), np.array([.5]), np.full(N - matches_null, .5 / (N - matches_null))])

    """
    # plot for viewing pleasure
    plt.scatter(range(N+1), prior, label='prior (comparison)', marker='o', color='b')
    plt.show()
    """

    # generate distributions
    dist_range = range(0, n + 1)

    # null dist is only affected by one point on prior
    null_dist = binom.pmf(dist_range, n, matches_null / N)

    if matches_null == N:
        # in case of matches > N, alt should just be binom
        alt_dist = binom.pmf(dist_range, n, matches_null / N)
    else:
        # otherwise alt dist affected by uniform prior
        alt_dist = np.empty(n + 1)
        for k in dist_range:
            for x in range(matches_null + 1, N + 1):
                alt_dist[k] += binom.pmf(k, n, x / N) * prior[x]
        # normalize
        alt_dist = alt_dist / sum(alt_dist)

    if (plot):
        # plot for viewing pleasure
        plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
        plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
        plt.show()

    return {
        'dist_range': dist_range,
        'alt_dist': alt_dist,
        'null_dist': null_dist
    }
 
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

    if (plot):
        # plot for viewing pleasure
        plt.scatter(dist_range, alt_dist, label='alt', marker='o', color='b')
        plt.scatter(dist_range, null_dist, label='null', marker='x', color='r')
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
    comparison_dists = generate_comparison_dists(n1, N_w1, N_l1, null_margin1, print_data=True)
    dist_range_comparison = comparison_dists['dist_range']
    alt_dist_comparison = comparison_dists['alt_dist']
    null_dist_comparison = comparison_dists['null_dist']

    # generate polling dists
    polling_dists = generate_polling_dists(n2, N_w2, N_l2, -1 * null_margin1)
    dist_range_polling = polling_dists['dist_range']
    alt_dist_polling = polling_dists['alt_dist']
    null_dist_polling = polling_dists['null_dist']

    # compute joint dists
    alt_joint = np.empty((len(dist_range_comparison), len(dist_range_polling)))
    null_joint = np.empty((len(dist_range_comparison), len(dist_range_polling)))
    for k_c in dist_range_comparison:
        for k_p in dist_range_polling:
            alt_joint[k_c][k_p] = alt_dist_comparison[k_c] * alt_dist_polling[k_p]
            null_joint[k_c][k_p] = null_dist_comparison[k_c] * null_dist_polling[k_p]

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

def maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=None, upper_bound=None, plot=False, print_data=False):
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
        comparison_dists = generate_comparison_dists(n1, N_w1, N_l1, null_margin1, print_data=True)
        alt_c = comparison_dists['alt_dist']
        null_c = comparison_dists['null_dist']

        # generate polling dists
        polling_dists = generate_polling_dists(n2, N_w2, N_l2, -1 * null_margin1)
        alt_p = polling_dists['alt_dist']
        null_p = polling_dists['null_dist']

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
        lower_bound = max_x - step_size
        upper_bound = max_x + step_size

        # perform refined search
        refined = maximize_joint_pvalue(k_c, k_p, N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound, upper_bound)

        # increase iterations by 1 and return results
        refined['search_iterations'] += 1
        return refined

def find_minimum_round_size(N_w1, N_l1, N_w2, N_l2, n1, stop_prob, alpha):
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

    n2 = 1
    while(1):
        # find kmax such that Pr[k >= kmax] = stop_prob
        kmax = math.floor(binom.ppf(stop_prob, n2, N_w2 / (N_w2 + N_l2)))

        # get pvalue for kmax
        pvalue = maximize_joint_pvalue(n1, kmax, N_w1, N_l1, N_w2, N_l2, n1, n2)['pvalue']

        # print n2 and pvalue for viewing pleasure
        #print("n2:",n2,"pvalue:",pvalue)

        # return n2 when pvalue for kmax meets stopping condition
        if (pvalue < alpha):
            return {
                "round_size": n2,
                "pvalue": pvalue
            }

        # increment round size
        n2 += 1





"""
From SUITE work (Ottobani et al)

Modified to allow for combination methods other than Fisher's

"""

import numpy as np
import math
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

def fisher_combined_pvalue(pvalues):
    """
    Find the p-value for Fisher's combined test statistic

    Parameters
    ----------
    pvalues : array_like
        Array of p-values to combine

    Returns
    -------
    float
        p-value for Fisher's combined test statistic
    """
    if np.any(np.array(pvalues)==0):
        return 0
    obs = -2*np.sum(np.log(pvalues))
    return 1-scipy.stats.chi2.cdf(obs, df=2*len(pvalues))


def create_modulus_old(n1, n2, n_w2, n_l2, N1, V_wl, gamma):
    """
    The modulus of continuity for the Fisher's combined p-value.
    This function returns the modulus of continuity, as a function of
    the distance between two lambda values.
    
    n1 : int
        sample size in the ballot comparison stratum
    n2 : int
        sample size in the ballot polling stratum
    n_w2 : int
        votes for the reported winner in the ballot polling sample
    n_l2 : int
        votes for the reported loser in the ballot polling sample
    N1 : int
        total number of votes in the ballot comparison stratum
    V_wl : int
        margin (in votes) between w and l in the whole contest
    gamma : float
        gamma from the ballot comparison audit
    """
    Wn = n_w2; Ln = n_l2; Un = n2-n_w2-n_l2
    assert Wn>=0 and Ln>=0 and Un>=0
    
    if N1 == 0:
        T2 = lambda delta: 0
    else:
        T2 = lambda delta: 2*n1*np.log(1 + V_wl*delta/(2*N1*gamma))
    
    return lambda delta: 2*Wn*np.log(1 + V_wl*delta) + 2*Ln*np.log(1 + 2*V_wl*delta) + \
            2*Un*np.log(1 + 3*V_wl*delta) + T2(delta)
            
            
def create_modulus(n1, n2, n_w2, n_l2, N1, V_wl, gamma):
    """
    The modulus of continuity for the Fisher's combined p-value.
    This function returns the modulus of continuity, as a function of
    the distance between two lambda values.
    
    n1 : int
        sample size in the ballot comparison stratum
    n2 : int
        sample size in the ballot polling stratum
    n_w2 : int
        votes for the reported winner in the ballot polling sample
    n_l2 : int
        votes for the reported loser in the ballot polling sample
    N1 : int
        total number of votes in the ballot comparison stratum
    V_wl : int
        margin (in votes) between w and l in the whole contest
    gamma : float
        gamma from the ballot comparison audit
    """
    Wn = n_w2; Ln = n_l2; Un = n2-n_w2-n_l2
    assert Wn>=0 and Ln>=0 and Un>=0
    
    if N1 == 0:
        T2 = lambda delta: 0
    else:
        T2 = lambda delta: 2*n1*np.log(1 + V_wl*delta/(2*N1*gamma))
    
    return lambda delta: 2*Wn*np.log(1 + V_wl*delta) + 2*Ln*np.log(1 + V_wl*delta) + \
            2*Un*np.log(1 + 2*V_wl*delta) + T2(delta)


def maximize_fisher_combined_pvalue(N_w1, N_l1, N1, N_w2, N_l2, N2,
    pvalue_funs, stepsize=0.05, modulus=None, alpha=0.05, feasible_lambda_range=None, combine_func=None):
    """
    Grid search to find the maximum P-value.

    Find the smallest Fisher's combined statistic for P-values obtained 
    by testing two null hypotheses at level alpha using data X=(X1, X2).

    Parameters
    ----------
    N_w1 : int
        votes for the reported winner in the ballot comparison stratum
    N_l1 : int
        votes for the reported loser in the ballot comparison stratum
    N1 : int
        total number of votes in the ballot comparison stratum
    N_w2 : int
        votes for the reported winner in the ballot polling stratum
    N_l2 : int
        votes for the reported loser in the ballot polling stratum
    N2 : int
        total number of votes in the ballot polling stratum
    pvalue_funs : array_like
        functions for computing p-values. The observed statistics/sample and known parameters should be
        plugged in already. The function should take the lambda allocation AS INPUT and output a p-value.
    stepsize : float
        size of the grid for searching over lambda. Default is 0.05
    modulus : function
        the modulus of continuity of the Fisher's combination function.
        This should be created using `create_modulus`.
        Optional (Default is None), but increases the precision of the grid search.
    alpha : float
        Risk limit. Default is 0.05.
    feasible_lambda_range : array-like
        lower and upper limits to search over lambda. 
        Optional, but a smaller interval will speed up the search.
    combine_func : function
        function used to combine strata pvalues
        Optional (Defaults to Fisher's method), but allows for trying others

    Returns
    -------
    dict with 

    max_pvalue: float
        maximum combined p-value
    min_chisq: float
        minimum value of Fisher's combined test statistic
    allocation lambda : float
        the parameter that minimizes the Fisher's combined statistic/maximizes the combined p-value
    refined : bool
        was the grid search refined after the first pass?
    stepsize : float
        the final grid step size used
    tol : float
        if refined is True, this is an upper bound on potential approximation error of min_chisq
    """	

    if (combine_func is None):
        combine_func = fisher_combined_pvalue

    assert len(pvalue_funs)==2
    
    # find range of possible lambda
    if feasible_lambda_range is None:
        feasible_lambda_range = calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2)
    (lambda_lower, lambda_upper) = feasible_lambda_range

    test_lambdas = np.arange(lambda_lower, lambda_upper+stepsize, stepsize)
    if len(test_lambdas) < 5:
        stepsize = (lambda_upper + 1 - lambda_lower)/5
        test_lambdas = np.arange(lambda_lower, lambda_upper+stepsize, stepsize)
    fisher_pvalues = np.empty_like(test_lambdas)
    pvalue1s = np.empty_like(test_lambdas)
    pvalue2s = np.empty_like(test_lambdas)
    for i in range(len(test_lambdas)):
        pvalue1 = float(np.min([1, pvalue_funs[0](test_lambdas[i])]))
        pvalue1s[i] = pvalue1
        pvalue2 = float(np.min([1, pvalue_funs[1](1-test_lambdas[i])]))
        pvalue2s[i] = pvalue2
        fisher_pvalues[i] = combine_func([pvalue1, pvalue2])
        
    pvalue = np.max(fisher_pvalues)
    max_index = np.argmax(fisher_pvalues)
    alloc_lambda = test_lambdas[max_index]
    pvalue1 = pvalue1s[max_index]
    pvalue2 = pvalue2s[max_index]

    
    """
    # REMOVED BECAUSE I WANT ACCURATE PVALUES EVEN WHEN GREATER THAN RISK LIMIT
    # If p-value is over the risk limit, then there's no need to refine the
    # maximization. We have a lower bound on the maximum.
    if pvalue > alpha or modulus is None:
        return {'max_pvalue' : pvalue,
                'pvalue1' : pvalue1,
                'pvalue2' : pvalue2,
                'min_chisq' : sp.stats.chi2.ppf(1 - pvalue, df=4),
                'allocation lambda' : alloc_lambda,
                'tol' : None,
                'stepsize' : stepsize,
                'refined' : False
                }
    """
    
    # Use modulus of continuity for the Fisher combination function to check
    # how close this is to the true max
    fisher_fun_obs = scipy.stats.chi2.ppf(1-pvalue, df=4)
    fisher_fun_alpha = scipy.stats.chi2.ppf(1-alpha, df=4)
    dist = np.abs(fisher_fun_obs - fisher_fun_alpha)
    mod = modulus(stepsize)

    if mod <= dist:
        return {'max_pvalue' : pvalue,
                'pvalue1' : pvalue1,
                'pvalue2' : pvalue2,
                'min_chisq' : fisher_fun_obs,
                'allocation lambda' : alloc_lambda,
                'stepsize' : stepsize,
                'tol' : mod,
                'refined' : False,
                'num_refined' : 0
                }
    else:
        lambda_lower = alloc_lambda - 2*stepsize
        lambda_upper = alloc_lambda + 2*stepsize
        refined = maximize_fisher_combined_pvalue(N_w1, N_l1, N1, N_w2, N_l2, N2,
            pvalue_funs, stepsize=stepsize/10, modulus=modulus, alpha=alpha, 
            feasible_lambda_range=(lambda_lower, lambda_upper))
        refined['refined'] = True
        refined['num_refined'] += 1

        return refined

def maximize_stouffers_combined_pvalue(N_w1, N_l1, N_w2, N_l2, n1, n2, pvalue_funs, stouffers, lower_bound=None, upper_bound=None, alpha=0.05):
    """
    Maximizes the stouffer combined pvalue for the given sample by searching
    over all possible allocations of winner votes under the null
    hypothesis.

    Parameters:
        N_w1 : reported winner votes in comparison stratum
        N_l1 : reported loser votes in comparison stratum
        N_w2 : reported winner votes in polling stratum
        N_lw : reported loser votes in polling stratum
        n1 : comparison sample size (first round)
        n2 : polling sample size (first round)
        lower_bound : lower bound for winner votes in comparison stratum under the null
        upper_bound : upper bound for winner votes in comparison stratum under the null
        stouffers : combination function with desired weights (pvalues as input)
        pvalue_funs : functions for computing the independent pvalues

    Return: dict with
        pvalue : maximum pvalue found
    """
    
    if lower_bound is None or upper_bound is None:
        # compute bounds for search
        bounds = compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2)
        lower_bound = bounds['x1_l']
        upper_bound = bounds['x1_u']

    """
    print("lower_bound:",lower_bound)
    print("upper_bound:",upper_bound)
    """

    # define a step size and generate lists for testing
    divisions = 10 # could tweak this for effiency
    step_size = math.ceil((upper_bound - lower_bound) / divisions)
    test_xs = np.arange(lower_bound, upper_bound + step_size, step_size)
    stouffers_pvalues = np.empty_like(test_xs, dtype=float)
    pvalue1s = np.empty_like(stouffers_pvalues)
    pvalue2s = np.empty_like(stouffers_pvalues)
    lambda1s = np.empty_like(stouffers_pvalues)

    # convert test_xs into lambdas for comparison stratum
    for i in range(len(test_xs)):
        # compute null margin based on winner votes
        x1 = test_xs[i]
        null_margin1 = x1 - (N_w1 + N_l1 - x1)

        # convert null_margin1 to lambda
        reported_margin_overall = (N_w1 + N_w2) - (N_l1 + N_l2)
        reported_margin1 = N_w1 - N_l1
        lambda1 = (null_margin1 - reported_margin1) / (reported_margin_overall)
        lambda1s[i] = lambda1
        #print("lambda1:",lambda1)

    for i in range(len(test_xs)):
        # get lambda
        lambda1 = lambda1s[i]

        # compute independent pvalues
        pvalue1 = np.min([1, pvalue_funs[0](lambda1)])
        pvalue1s[i] = pvalue1
        pvalue2 = np.min([1, pvalue_funs[1](1-lambda1)])
        pvalue2s[i] = pvalue2
        stouffers_pvalues[i] = stouffers([pvalue1, pvalue2])
  
        """
        # print for viewing pleasure
        print("N_w1:",N_w1)
        print("N_l1:",N_l1)
        print("N_w2:",N_w2)
        print("N_l2:",N_l2)
        print("x1:",x1)
        print("null_margin1:",null_margin1)
        print("comparison pvalue:",pvalue1)
        print("polling pvalue:",pvalue2)
        print("stouffers_combined:",stouffers_pvalues[i])
        print("i+1:",i+1,"of",len(test_xs))
        """

    # get the maximum pvalue found
    max_index = np.argmax(stouffers_pvalues)
    max_pvalue = stouffers_pvalues[max_index]
    max_x = test_xs[max_index]
    pvalue1 = pvalue1s[max_index]
    pvalue2 = pvalue2s[max_index]

    print("comparison pvalue:",pvalue1)
    print("polling pvalue:",pvalue2)
    print("stouffers_combined:",stouffers_pvalues[i])
 
    """
    print("max_index:",max_index)
    print("max_pvalue:",max_pvalue)
    print("max_x:",max_x)
    """

    # if step_size has reached 1, search is over
    if step_size == 1:
        return {
            'max_pvalue' : max_pvalue,
            'pvalue1':pvalue1,
            'pvalue2':pvalue2,
            'x1' : max_x,
            'search_iterations' : 1,
            'allocation lambda' : None
        }
    else:
        # set bounds for refined search
        x1_l = max_x - step_size
        x1_u = max_x + step_size

        # perform refined search
        refined = maximize_stouffers_combined_pvalue(N_w1, N_l1, N_w2, N_l2, n1, n2, lower_bound=x1_l, upper_bound=x1_u, pvalue_funs=pvalue_funs, stouffers=stouffers)

        # increase iterations by 1 and return results
        refined['search_iterations'] += 1
        return refined

def compute_winner_vote_bounds(N_w1, N_l1, N_w2, N_l2, k_c=0, k_p=0):
    """
    Computes upper and lower bounds for the number of winner votes in 
    the comparison stratum, x1. 

    Parameters:
        N_w1 : winner votes in comparison stratum
        N_l1 : loser votes in comparison stratum
        N_w2 : winner votes in polling stratum
        N_lw : loser votes in polling stratum
        k_c : matches in comparison sample
            Optional: defaults to zero
        k_p : winner votes in polling sample
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



def plot_fisher_pvalues(N, overall_margin, pvalue_funs, alpha=None):
    """
    Plot the Fisher's combined p-value for varying error allocations 
    using data X=(X1, X2) 

    Parameters
    ----------
    N : array_like
        Array of stratum sizes
    overall_margin : int
        the difference in votes for reported winner and loser in the population
    pvalue_funs : array_like
        functions for computing p-values. The observed statistics/sample and known parameters should be plugged in already. The function should take the lambda allocation AS INPUT and output a p-value.
    alpha : float
        Optional, desired upper percentage point
    
    Returns
    -------
    dict with 
    
    float
        maximum combined p-value
    float
        minimum value of Fisher's combined test statistic
    float
        lambda, the parameter that minimizes the Fisher's combined statistic/maximizes the combined p-value
    """
    assert len(N)==2
    assert len(pvalue_funs)==2
        
    # find range of possible lambda
    (lambda_lower, lambda_upper) = calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2)

    fisher_pvalues = []
    cvr_pvalues = []
    nocvr_pvalues = []
    for lam in np.arange(lambda_lower, lambda_upper+1, 0.5):
        pvalue1 = np.min([1, pvalue_funs[0](lam)])
        pvalue2 = np.min([1, pvalue_funs[1](1-lam)])
        fisher_pvalues.append(fisher_combined_pvalue([pvalue1, pvalue2]))

    plt.scatter(np.arange(lambda_lower, lambda_upper+1, 0.5), fisher_pvalues, color='black')
    if alpha is not None:
        plt.axhline(y=alpha, linestyle='--', color='gray')
    plt.xlabel("Allocation of Allowable Error")
    plt.ylabel("Fisher Combined P-value")
    plt.show()
    
def calculate_lambda_range(N_w1, N_ell1, N_1, N_w2, N_ell2, N_2):
    '''
    Find the largest and smallest possible values of lambda.
    
    Input:
    ------
        N_w1 : int
            reported votes for overall reported winner w in stratum 1
        N_ell1 : int
            reported votes for overall reported loser ell in stratum 1
        N1 : int
            ballots cast in stratum 1
        N_w2 : int
            reported votes for overall reported winner w in stratum 2
        N_ell2 : int
            reported votes for overall reported loser ell in stratum 2
        N1 : int
            ballots cast in stratum 2
   
    Returns:
    --------
        (lb, ub): real ordered pair. lb is a sharp lower bound on lambda; ub is a sharp upper bound
    
    Derivation:
    -----------
    Let V denote the overall reported margin of w over ell across both strata, i.e., 
        V = (N_w1 + N_w2) - (N_ell1 + N_ell2)
        
    The overstatement error in votes in stratum s is *at most* the difference between the 
    reported margin in that stratum and the margin that would result if all ballots in 
    stratum s had votes for the loser; i.e.,
       margin_error_s <= (N_ws - N_ells)+N_s.
    Thus 
        lambda*V <= N_w1 - N_ell1 + N_1.
        (1-lambda)*V <= N_w2 - N_ell2 + N_2, i.e., lambda*V >= V - (N_w2 - N_ell2 + N_2)
    
    The overstatement error in votes in a stratum is at least the difference between the 
    reported margin in that stratum and the margin that would result if all ballots in the
    stratum had votes for the winner; i.e.,
       margin_error_s >= (N_ws - N_ells)-N_s.
    Thus
       lambda*V >=  N_w1 - N_ell1 - N_1 
       (1 - lambda)*V >=  N_w2 - N_ell2 - N_2, i.e., lambda*V <= V - (N_w2 - N_ell2 - N_2)
    
    Combining these yields
       lambda >= max( N_w1 - N_ell1 - N_1, V - (N_w2 - N_ell2 + N_2) )/V
       lambda <= min( N_w1 - N_ell1 + N_1, V - (N_w2 - N_ell2 - N_2) )/V.
    '''
    V = N_w1 + N_w2 - N_ell1 - N_ell2
    lb = np.amax([N_w1 - N_ell1 - N_1, V - (N_w2 - N_ell2 + N_2)] )/V
    ub = np.amin([ N_w1 - N_ell1 + N_1, V - (N_w2 - N_ell2 - N_2)] )/V       
    return (lb, ub)
    
    
def bound_fisher_fun(N_w1, N_l1, N1, N_w2, N_l2, N2,
                     pvalue_funs, feasible_lambda_range=None, stepsize=0.05):
        """
        DEPRECATED: Create piecewise constant upper and lower bounds for the Fisher's 
        combination function for varying error allocations

        Parameters
        ----------
        N_w1 : int
            votes for the reported winner in the ballot comparison stratum
        N_l1 : int
            votes for the reported loser in the ballot comparison stratum
        N1 : int
            total number of votes in the ballot comparison stratum
        N_w2 : int
            votes for the reported winner in the ballot polling stratum
        N_l2 : int
            votes for the reported loser in the ballot polling stratum
        N2 : int
            total number of votes in the ballot polling stratum
        pvalue_funs : array_like
            functions for computing p-values. The observed statistics/sample and known parameters
            should be plugged in already. 
            The function should take the lambda allocation AS INPUT and output a p-value.
        feasible_lambda_range : array-like
            lower and upper limits to search over lambda. Optional, but will speed up the search
        stepsize : float
            size of the mesh to calculate bounds; default 0.5
    
        Returns
        -------
        dict with 
    
        array-like
           sample_points : Fisher's combining function evaluated at the grid points
        array-like
           lower_bounds : piecewise constant lower bound on Fisher's combining function between the grid points
        array-like
           upper_bounds : piecewise constant upper bound on Fisher's combining function between the grid points
        array-like
           grid : grid of lambdas
        """
        if feasible_lambda_range is None:
            feasible_lambda_range = calculate_lambda_range(N_w1, N_l1, N1, N_w2, N_l2, N2)
        (lambda_lower, lambda_upper) = feasible_lambda_range
        
        cvr_pvalue = pvalue_funs[0]
        nocvr_pvalue = pvalue_funs[1]
        cvr_pvalues = []
        nocvr_pvalues = []
        
        for lam in np.arange(lambda_lower, lambda_upper+1, stepsize):
            pvalue1 = np.min([1, cvr_pvalue(lam)])
            pvalue2 = np.min([1, nocvr_pvalue(1-lam)])
            cvr_pvalues.append(pvalue1)
            nocvr_pvalues.append(pvalue2)
            
        lower_bounds = [fisher_combined_pvalue([cvr_pvalues[i+1], nocvr_pvalues[i]]) for i in range(len(cvr_pvalues)-1)]
        upper_bounds = [fisher_combined_pvalue([cvr_pvalues[i], nocvr_pvalues[i+1]]) for i in range(len(cvr_pvalues)-1)]
        sample_points = [fisher_combined_pvalue([cvr_pvalues[i], nocvr_pvalues[i]]) for i in range(len(cvr_pvalues))]
        
        return {'sample_points' : sample_points,
                'upper_bounds' : upper_bounds,
                'lower_bounds' : lower_bounds,
                'grid' : np.arange(lambda_lower, lambda_upper+1, stepsize)
                }


################################################################################
############################## Unit tests ######################################
################################################################################

def test_modulus1():
    N1 = 1000
    N2 = 100
    Vw1 = 550
    Vl1 = 450
    Vw2 = 60
    Vl2 = 40
    Vu2 = N2-Vw2-Vl2
    Vwl = (Vw1 + Vw2) - (Vl1 + Vl2)

    Nw1_null = N1/2
    n1 = 100
    o1 = o2 = u1 = u2 = 0
    gamma = 1.03

    nw2 = 6; Wn = nw2
    nl2 = 4; Ln = nl2
    n2 = 10; Un = n2 - nl2 - nw2

    def c(lam):
        return Vw2 - Vl2 - lam*Vwl

    def Nw_null(lam):
        return (N2 - Un + c(lam))/2

    def fisher_fun(lam):
        Nw_null_val = Nw_null(lam)
        c_val = c(lam)
        fisher_fun = -2*np.sum(np.log(Nw_null_val - np.arange(Wn))) +\
        -2*np.sum(np.log(Nw_null_val - c_val - np.arange(Ln))) +\
        -2*np.sum(np.log(N2 - 2*Nw_null_val + c_val - np.arange(Un))) +\
        2*np.sum(np.log(Vw2 - np.arange(Wn))) +\
        2*np.sum(np.log(Vl2 - np.arange(Ln))) +\
        2*np.sum(np.log(Vu2 - np.arange(Un))) +\
        -2*n1*np.log(1 - (1-lam)*Vwl/(2*N1*gamma)) + \
        2*o1*np.log(1 - 1/(2*gamma)) + \
        2*o2*np.log(1 - 1/(gamma)) + \
        2*u1*np.log(1 + 1/(2*gamma)) + \
        2*u2*np.log(1 + 1/(gamma))
        return fisher_fun
        
    mod = create_modulus(n1, n2, nw2, nl2, N1, Vwl, gamma)
    
    v1 = np.abs(fisher_fun(0.6 + 0.1) - fisher_fun(0.6))
    v2 = mod(0.1)
    np.testing.assert_array_less(v1, v2)
    v1 = np.abs(fisher_fun(0.2 + 0.01) - fisher_fun(0.2))
    v2 = mod(0.01)
    np.testing.assert_array_less(v1, v2)
    v1 = np.abs(fisher_fun(0.8 + 0.001) - fisher_fun(0.8))
    v2 = mod(0.001)
    np.testing.assert_array_less(v1, v2)
    
    
    mod_old = create_modulus_old(n1, n2, nw2, nl2, N1, Vwl, gamma)
    np.testing.assert_array_less(mod(0.1), mod_old(0.1))
    np.testing.assert_array_less(mod(0.01), mod_old(0.01))
    np.testing.assert_array_less(mod(0.001), mod_old(0.001))


if __name__ == "__main__":
    test_modulus1()

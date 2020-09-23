from __future__ import division
import math
import numpy as np
import numpy.random
import scipy as sp
import scipy.stats
import scipy.optimize


def ballot_polling_sprt(sample, popsize, alpha, Vw, Vl, 
                        null_margin=0, number_invalid=None):
    """
    Wald's SPRT for the difference in true number of votes for the winner, Nw, and the loser, Nl:

    H_0: Nw = Nl + null_margin
    H_1: Nw = Vw, Nl = Vl

    The type II error rate, usually denoted beta, is set to 0%:
    if the data do not support rejecting the null, there is a full hand count.
    Because beta=0, the reciprocal of the likelihood ratio is a conservative p-value.

    Parameters
    ----------
    sample : array-like
        audit sample. Elements should equal 1 (ballots for w), 0 (ballots for l), or np.nan (the rest)
    popsize : int
        number of ballots cast in the stratum
    alpha : float
        desired type 1 error rate
    Vw : int
        number of votes for w in the stratum, under the alternative hypothesis
    Vl : int
        total number of votes for l in the stratum, under the alternative hypothesis
    null_margin : int
        vote margin between w and l under the null hypothesis; optional
        (default 0)
    number_invalid : int
        total number of invalid  ballots, undervoted ballots, or ballots for 
        other candidates in the stratum, if known; optional (default None)

    Returns
    -------
    dict
    """
    
    # Set parameters
    upper = 1/alpha
    n = len(sample)
    assert n <= popsize, "Sample size greater than population size"
    sample = np.array(sample)
    Wn = np.sum(sample == 1)
    Ln = np.sum(sample == 0)
    Un = n - Wn - Ln
    decision = "None"

    # Set up likelihood for null and alternative hypotheses
    Vw = int(Vw)
    Vl = int(Vl)
    Vu = int(popsize - Vw - Vl)
    assert Vw >= Wn and Vl >= Ln and Vu >= Un, "Alternative hypothesis isn't consistent with the sample"

    alt_logLR = np.sum(np.log(Vw - np.arange(Wn))) + \
                np.sum(np.log(Vl - np.arange(Ln))) + \
                np.sum(np.log(Vu - np.arange(Un)))
        
    null_logLR = lambda Nw: (Wn > 0)*np.sum(np.log(Nw - np.arange(Wn))) + \
                (Ln > 0)*np.sum(np.log(Nw - null_margin - np.arange(Ln))) + \
                (Un > 0)*np.sum(np.log(popsize - 2*Nw + null_margin - np.arange(Un)))
    
    # This is for testing purposes. In practice, number_invalid will be unknown.
    if number_invalid is not None:
        assert isinstance(number_invalid, int)
        assert number_invalid < popsize
        nuisance_param = (popsize - number_invalid + null_margin)/2
        if nuisance_param < Wn or (nuisance_param - null_margin) < Ln \
            or number_invalid < Un:
            return {'decision' : 'Number invalid is incompatible with the null and the data',
                    'upper_threshold' : upper,
                    'LR' : np.inf,
                    'pvalue' : 0,
                    'sample_proportion' : (Wn/n, Ln/n, Un/n),
                    'Nu_used' : number_invalid,
                    'Nw_used' : nuisance_param
                    }
        res = alt_logLR - null_logLR(nuisance_param)
        LR = np.exp(res)

    else:
        upper_Nw_limit = (popsize - Un + null_margin)/2
        lower_Nw_limit = np.max([Wn, Ln+null_margin])
        
        # For extremely small or large null_margins, the limits do not
        # make sense with the sample values.
        if upper_Nw_limit < Wn or (upper_Nw_limit - null_margin) < Ln:
            return {'decision' : 'Null is impossible, given the sample',
                    'upper_threshold' : upper,
                    'LR' : np.inf,
                    'pvalue' : 0,
                    'sample_proportion' : (Wn/n, Ln/n, Un/n),
                    'Nu_used' : None,
                    'Nw_used' : None
                    }
        
        if lower_Nw_limit > upper_Nw_limit:
            lower_Nw_limit, upper_Nw_limit = upper_Nw_limit, lower_Nw_limit
        
        LR_derivative = lambda Nw: np.sum([1/(Nw - i) for i in range(Wn)]) + \
                    np.sum([1/(Nw - null_margin - i) for i in range(Ln)]) - \
                    2*np.sum([1/(popsize - 2*Nw + null_margin - i) for i in range(Un)])

        # Sometimes the upper_Nw_limit is too extreme, causing illegal 0s.
        # Check and change the limit when that occurs.
        if np.isinf(null_logLR(upper_Nw_limit)) or np.isinf(LR_derivative(upper_Nw_limit)):
            upper_Nw_limit -= 1

        # Check if the maximum occurs at an endpoint: deriv has no sign change
        if LR_derivative(upper_Nw_limit)*LR_derivative(lower_Nw_limit) > 0: 
            nuisance_param = upper_Nw_limit if null_logLR(upper_Nw_limit)>=null_logLR(lower_Nw_limit) else lower_Nw_limit
            #print(nuisance_param)
        # Otherwise, find the (unique) root of the derivative of the log likelihood ratio
        else:
            nuisance_param = sp.optimize.brentq(LR_derivative, lower_Nw_limit, upper_Nw_limit)
        number_invalid = popsize - nuisance_param*2 + null_margin
        logLR = alt_logLR - null_logLR(nuisance_param)
        LR = np.exp(logLR)

    if LR <= 0:
        # accept the null and stop
        decision = 0
                
    if LR >= upper:
        # reject the null and stop
        decision = 1
            
    return {'decision' : decision,
            'upper_threshold' : upper,
            'LR' : LR,
            'pvalue' : min(1, 1/LR),
            'sample_proportion' : (Wn/n, Ln/n, Un/n),
            'Nu_used' : number_invalid,
            'Nw_used' : nuisance_param
            }


###################### Unit tests ############################

def test_sprt_functionality():
    trials = np.zeros(100)
    trials[0:50] = 1
    res = ballot_polling_sprt(trials, popsize=1000, alpha=0.05, Vw=500, Vl=450)
    assert res['decision']=='None'
    assert res['lower_threshold']==0.0
    assert res['upper_threshold']==20.0
    assert res['pvalue']>0.05
    assert res['LR']<1
    assert res['sample_proportion']==(0.5, 0.5, 0)
    
    trials[50:60] = 1
    res = ballot_polling_sprt(trials, popsize=1000, alpha=0.05, Vw=600, Vl=400)
    assert res['decision']=='None'
    assert res['lower_threshold']==0.0
    assert res['upper_threshold']==20.0
    assert res['pvalue']>0.05
    assert res['LR']>1
    assert res['sample_proportion']==(0.6, 0.4, 0)
    
    trials = np.zeros(250)
    trials[0:110] = 1
    trials[110:150] = np.nan
    res = ballot_polling_sprt(trials, popsize=1000, alpha=0.05, Vw=500, Vl=450)
    assert res['decision']=='None'
    assert res['lower_threshold']==0.0
    assert res['upper_threshold']==20.0
    assert res['pvalue']>0.05
    assert res['LR']<1
    assert res['sample_proportion']==(0.44, 0.4, 0.16)
            
    trials = np.zeros(100)
    trials[0:40] = 1        
    res = ballot_polling_sprt(trials, popsize=1000, alpha=0.05, Vw=500, Vl=450)
    assert res['decision']=='None'
    assert res['pvalue']==1
    assert res['LR']<1
    assert res['sample_proportion']==(0.4, 0.6, 0)


def test_sprt_analytic_example():
    sample = [0, 0, 1, 1]
    population = [0]*5 + [1]*5
    popsize = len(population)
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=5, Vl=4, number_invalid=0)
    np.testing.assert_almost_equal(res['LR'], 0.6)
    np.testing.assert_almost_equal(res['Nu_used'], 0)
    res2 = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=5, Vl=4)
    np.testing.assert_almost_equal(res2['LR'], 0.6, decimal=2)
    np.testing.assert_almost_equal(res2['Nu_used'], 0)
    
    sample = [0, 1, 1, 1]
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4, number_invalid=0)
    np.testing.assert_almost_equal(res['LR'], 1.6)
    res2 = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4)
    np.testing.assert_almost_equal(res2['LR'], 1.6, decimal=2)
    np.testing.assert_almost_equal(res2['Nu_used'], 0, decimal=2)

    sample = [0, 1, 1]
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4, number_invalid=0)
    np.testing.assert_almost_equal(res['LR'], 1.2)
    res2 = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4)
    np.testing.assert_almost_equal(res2['LR'], 1.2, decimal=2)
    np.testing.assert_almost_equal(res2['Nu_used'], 0, decimal=2)


def test_edge_cases():
    sample = [0, 0, 1, 1]
    population = [0]*5 + [1]*5
    popsize = len(population)

    # if nuisance_param < 0 or nuisance_param > popsize
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4, number_invalid=12)
    np.testing.assert_almost_equal(res['decision'], 'Number invalid is incompatible with the null and the data')

    # if nuisance_param < Wn or (nuisance_param - null_margin) < Ln \
    #    or number_invalid < Un:
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4, number_invalid=7)
    np.testing.assert_almost_equal(res['decision'], 'Number invalid is incompatible with the null and the data')

    # if upper_Nw_limit < Wn or (upper_Nw_limit - null_margin) < Ln:
    res = ballot_polling_sprt(sample, popsize, alpha=0.05, Vw=6, Vl=4, null_margin=7)
    np.testing.assert_almost_equal(res['decision'], 'Null is impossible, given the sample')


if __name__ == 'main':
    test_sprt_functionality()
    test_sprt_analytic_example()
    test_edge_cases()

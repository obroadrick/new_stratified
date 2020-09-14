# new_stratified
My work developing and testing a new stratified audit approach.

## Context (my previous work)
Summer 2020, I worked on combining [SUITE](https://github.com/pbstark/CORLA18/tree/master/code) and [Athena](https://github.com/gwexploratoryaudits/r2b2/tree/master) for 2-strata audits. That work is in a [separate repo](https://github.com/obroadrick/stratified_athena/tree/master).

## The new approach
The new idea is to avoid the conservative p-values we get from using combining functions to combine distinct strata p-values. We will do this by combining probability distributions from individual strata into one joint probability distribution. Then, a stopping condition, similar to that presented in [Athena](https://github.com/gwexploratoryaudits/r2b2/tree/master), will be used.



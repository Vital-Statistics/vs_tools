"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
# Auto-generated from vsUtilityFunctions/ package modules.

# ---- source: vsUtilityFunctions/delong_test.py ----
#!/usr/bin/env python3

######### THIS IS PRODUCING WRONG RESULTS

# def delong_test(y_true, y_pred1, y_pred2):
#     from scipy.stats import norm
#     import numpy as np
#     from sklearn.metrics import roc_auc_score
#     """ DeLong test for comparing the difference between two AUCs """
#     def compute_delong_auc_variance(auc, n_pos, n_neg):
#         q1 = auc / (2 - auc)
#         q2 = 2 * auc**2 / (1 + auc)
#         return (auc * (1 - auc) + (n_pos - 1) * (q1 - auc**2) + (n_neg - 1) * (q2 - auc**2)) / (n_pos * n_neg)

#     auc1 = roc_auc_score(y_true, y_pred1)
#     auc2 = roc_auc_score(y_true, y_pred2)
    
#     n_pos = np.sum(y_true)
#     n_neg = len(y_true) - n_pos

#     var1 = compute_delong_auc_variance(auc1, n_pos, n_neg)
#     var2 = compute_delong_auc_variance(auc2, n_pos, n_neg)
#     cov12 = np.sqrt(var1 * var2)  # Assuming independence

#     diff = auc1 - auc2
#     se = np.sqrt(var1 + var2 - 2 * cov12)
#     z = diff / se
#     p_value = 2 * (1 - norm.cdf(abs(z)))

#     print('AUC 1: '+str(round(auc1,4)))
#     print('AUC 2: '+str(round(auc2,4)))
#     print('z: '+str(round(z,4)))
#     print('p-value: '+str(round(p_value,4)))

#     return auc1, auc2, z, p_value

# import numpy as np
# from scipy.stats import norm
# from sklearn.metrics import roc_auc_score

# def compute_auc_variance(ground_truth, predictions):
#     """
#     Computes DeLong variance for AUC
#     """
#     order = np.argsort(-predictions)  # Sort in descending order
#     ranked_truth = ground_truth[order]
    
#     pos_scores = np.cumsum(ranked_truth == 1) / np.sum(ranked_truth == 1)
#     neg_scores = np.cumsum(ranked_truth == 0) / np.sum(ranked_truth == 0)
    
#     auc = roc_auc_score(ground_truth, predictions)
#     var = (auc * (1 - auc) + (len(pos_scores) - 1) * np.var(pos_scores) + (len(neg_scores) - 1) * np.var(neg_scores)) / (len(pos_scores) * len(neg_scores))
    
#     return auc, var

# def delong_test(y_true, y_pred1, y_pred2):
#     """ Performs DeLong test for two correlated ROC AUC scores """
#     auc1, var1 = compute_auc_variance(np.array(y_true), np.array(y_pred1))
#     auc2, var2 = compute_auc_variance(np.array(y_true), np.array(y_pred2))
    
#     cov12 = np.sqrt(var1 * var2)  # Covariance term, approximated conservatively
#     diff = auc1 - auc2
#     se = np.sqrt(var1 + var2 - 2 * cov12)  # Standard error
    
#     z = diff / se  # Z-score
#     p_value = 2 * (1 - norm.cdf(abs(z)))  # Two-tailed p-value
    
#     print('AUC 1: '+str(round(auc1,4)))
#     print('AUC 2: '+str(round(auc2,4)))
#     print('z: '+str(round(z,4)))
#     print('p-value: '+str(round(p_value,4)))

#     return auc1, auc2, z, p_value


# ---- source: vsUtilityFunctions/rhoToZ.py ----


def rhoToZ(rho):
    import math
    return(.5*math.log((1+rho)/(1-rho)))

# ---- source: vsUtilityFunctions/standardize.py ----


import pandas as pd

def standardize(x,axis=0):
    import numpy as np
    col=list(x)
    if axis==1:
        z=(np.array(x)-np.array(x).mean(axis=1).reshape(-1,1))/np.array(x).std(axis=1).reshape(-1,1) 
    else:
        z=(np.array(x)-np.array(x).mean(axis=0))/np.array(x).std(axis=0)
    if type(x)==pd.DataFrame:
        z=pd.DataFrame(z,columns=list(x),index=x.index)
    return(z)

# ---- source: vsUtilityFunctions/twoCorrTest.py ----


############ Reference:
###########  Hinkle DE, Wiersma W, Jurs SG (1988) Applied statistics for the behavioral sciences. 2nd ed. Boston: Houghton Mifflin Company
############ Fisher r to z plus t-test


def twoCorrTest(rho0,n0,rho1,n1):
    import math
    from scipy.stats import norm
    z=(rhoToZ(rho0)-rhoToZ(rho1))/math.sqrt(1/(n0-3)+1/(n1-3))
    p=(1-norm.cdf(abs(z)))*2
    # p=norm.logcdf(abs(z))
    return((z,p))

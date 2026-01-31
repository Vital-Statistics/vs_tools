"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
# Auto-generated from vsPosteriorDensities/ package modules.

# ---- source: vsPosteriorDensities/ll_MVNormal.py ----


import numpy as np

def ll_MVNormal(X,prior=1e-10):
    import scipy.stats as stats
    mvn=stats.multivariate_normal
    
    n,c=X.shape
    sig=np.array(X.transpose())@np.array(X)/n+prior*np.identity(c)

    return(sum(np.log(mvn.pdf(X,cov=sig))))

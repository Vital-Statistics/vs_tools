# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 13:41:55 2021

@author: joest
"""

import numpy as np

def ll_MVNormal(X,prior=1e-10):
    import scipy.stats as stats
    mvn=stats.multivariate_normal
    
    n,c=X.shape
    sig=np.array(X.transpose())@np.array(X)/n+prior*np.identity(c)

    return(sum(np.log(mvn.pdf(X,cov=sig))))

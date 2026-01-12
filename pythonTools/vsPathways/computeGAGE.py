# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 09:11:43 2021

@author: joest
"""

import math
import pandas as pd

def computeGAGE(pv):
    ### pv should be a list of p-values for analytes in the pathway
    import scipy.stats as stats
    pv=[max(1e-30,v) for v in pv if not pd.isna(v)]
    gam=stats.gamma(len(pv))
    return(max(1e-30,1-gam.cdf(sum([-math.log(v) for v in pv]))))


# ############## Testing for accuracy
# nSamp=10000
# sampSize=10
# gam=stats.gamma(sampSize)
# res=pd.Series(0,index=list(range(nSamp)))
# for i in res.index.values:
#     rr=np.random.uniform(size=sampSize)
#     # rr=-np.log(np.random.uniform(size=sampSize))
#     # rr=sum(rr)
#     # rr=np.mean(gam.rvs(100))
#     # res.loc[i]= 1-gam.cdf(rr)
#     res.loc[i]=computeGAGE(rr)
# res.hist(bins=np.arange(0,1.01,.01))


# np.mean(gam.rvs(10000))

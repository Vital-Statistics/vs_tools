# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:29:12 2017

@author: jel2
"""

def v_buildSparse(x,y,delta=None):
    import scipy
    import pandas as pd
    import numpy as np
    if delta is None:
        delta=[1]*len(x)
    i=x.codes
    j=y.codes
    rLbl=pd.Series(x.categories)
    cLbl=pd.Series(y.categories)
    M=scipy.sparse.coo_matrix((np.array(delta),(i,j)),shape=(len(rLbl),len(cLbl)))
    M=M.tocsr()
    
    return(pd.DataFrame.sparse.from_spmatrix(M,columns=cLbl,index=rLbl))

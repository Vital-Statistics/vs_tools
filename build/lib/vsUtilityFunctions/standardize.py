# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:33:58 2022

@author: joest
"""

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

# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:05:36 2017

@author: jel2
"""

import pandas as pd

#### this version resets the index of X in the returned table
def getSubtable(X,dictColumn='DETAILS'):
    X=X.reset_index()
    t=list(X)
    t.remove(dictColumn)
    return(X[t].join(pd.DataFrame.from_records(X.DETAILS.tolist(),index=X.index)))
    
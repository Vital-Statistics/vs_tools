#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 12:46:51 2023

@author: rudy
"""

import numpy as np
import pandas as pd

def predict_proba(mod,Z):
    ll=(Z.isna().sum(axis=1)==0)
    yhat=pd.Series(np.nan,index=Z.index)
    yhat[ll]=mod.predict_proba(Z.loc[ll])[:,1]
    return(yhat)

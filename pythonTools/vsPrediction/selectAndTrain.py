#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 09:07:54 2023

@author: rudy
"""

import numpy as np
import pandas as pd

from vsPrediction.xval_train import xval_train
from vsVisualizations.loopProgress import loopProgress

def selectAndTrain(X,y,lbl='Setting regression penalty',returnYhat=False,verbose=False,leaveOneOut=False):
    from sklearn.metrics import roc_auc_score
    from scipy import stats as stat
    from scipy.special import logit
    from statsmodels.sandbox.stats.multicomp import multipletests
    from sklearn.linear_model import LogisticRegression

    ### Set Regression Penalty with cross validation
    nNan=pd.concat([X,y],axis=1).isna().sum(axis=1)
    nDrop=(nNan>0).sum()
    if nDrop>0:
        print('Dropping '+str(nDrop)+' observations due to missingness.')
        X=X.loc[nNan==0]
        y=y.loc[nNan==0]
    xvp=pd.DataFrame({'tt':np.arange(.025,2.0,.025),'cstat':.5,'P-Value':1.0,'FDR':1.0})
    for i,tt in enumerate(xvp.tt):
        loopProgress(i,len(xvp),lbl)
        mod=LogisticRegression(max_iter=1000,penalty='l1',solver='liblinear',C=tt)
        yhat=xval_train(mod,X,y,leaveOneOut=leaveOneOut)
        
        xvp.loc[xvp.tt==tt,'P-Value']=stat.ranksums(yhat.loc[y==0],yhat.loc[y==1],alternative='less')[1]  ## One sided test since the model is designed to assign higher values when outcome is 1
        
        xvp.loc[xvp.tt==tt,'cstat']=roc_auc_score(y,yhat)
        if xvp.loc[i,'cstat']==xvp.cstat.max():
            bestYHat=yhat
    xvp.sort_values('cstat',inplace=True,ascending=False)
    
    ### Retrain the model using the identified regression penalty
    xvp['FDR']=multipletests(xvp['P-Value'],method='fdr_bh')[1]
    tt=xvp.iloc[0].tt
    
    mod=LogisticRegression(max_iter=1000,penalty='l1',solver='liblinear',C=tt)
    mod.fit(X,y)
    
    ## identify the variables that were used
    aUsd=pd.DataFrame({'Analyte':list(X),'beta':mod.coef_.squeeze()})
    kpZ=list(aUsd.loc[aUsd.beta!=0,'Analyte'])
    
    if verbose:
        print(xvp.head(10))
    
    if returnYhat:
        res=(mod,kpZ,bestYHat)
    else:
        res=(mod,kpZ)

    ########## Accuracy on validation data
    return(res)

# We used logistic regression with a lasso penalty to simultaneously perform variable selection
# and model fitting.  We tested a range of penalties, in each case testing accuracy on hold outs
# using leave-one-out cross-validation.  All performance statistics are on the highest performing
# model, but are computed on the held out samples only.  Variables identified as included in the
# model are those that had non-zero regression coefficients when using the optimal penalty on the
# full data set.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 08:48:39 2023

@author: rudy
"""

def regressAll(vIndep,Y,studyVars=None,lbl=None,showIntercept=False):
    import numpy as np
    import statsmodels.api as sm
    from statsmodels.sandbox.stats.multicomp import multipletests
    
    # if type(dem.Group)==pd.core.series.Series:
    #     vIndep=pd.DataFrame(vIndep)
    
    drp=(Y.std()<1e-10)
    if drp.sum()>0:
        print('Dropping '+str(drp.sum())+' columns from dependent variable list because they are almost constant.')
        Y=Y.drop(columns=drp.loc[drp].index.values)
    if lbl is None:
        lbl='Regressing'
    if studyVars is None:
        studyVars=list(vIndep)
    if type(studyVars)!=list:
        print('studyVars must be a list of column headers')
        return
    drp=vIndep.isna().sum(axis=1)>0
    if drp.sum()>0:
        print('Dropping '+str(drp.sum())+' sample(s) due to missingness in independent variables.')
        print(vIndep.isna().sum())
    vIndep=vIndep.loc[~drp]
    ### Create a blank dataframe to contain results
    blank={'# Missing':np.nan}
    for v in studyVars:
        blank['beta, '+v]=np.nan
        blank['P-Value, '+v]=np.nan
        blank['95 %tile, '+v]=''
        blank['FDR, '+v]=np.nan
    for v in [b for b in list(vIndep) if b not in studyVars]:
        blank['beta, '+v]=np.nan
        blank['P-Value, '+v]=np.nan
    res=pd.DataFrame(blank,index=list(Y))
    N=len(res)
    Z=vIndep.join(Y)
    X_Const=sm.add_constant(Z[list(vIndep)],has_constant='raise')
    for i,col in enumerate(res.index.values):
        loopProgress(i,N,lbl)
        kp=~Z[col].isna()
        if X_Const.loc[kp,studyVars].std().max()>0:
            res.loc[col,'# Missing']=sum(~kp)
            mod=sm.OLS(Z.loc[kp,col],X_Const.loc[kp]).fit()
            res.loc[col,['beta, '+v for v in list(vIndep)]]=mod.params[list(vIndep)].values
            res.loc[col,['P-Value, '+v for v in list(vIndep)]]=mod.pvalues[list(vIndep)].values
            res.loc[col,['95 %tile, '+v for v in studyVars]]=mod.conf_int().loc[studyVars].apply(lambda x:'['+str(round(x[0],3))+', '+str(round(x[1],3))+']',axis=1).values
            if showIntercept:
                res.loc[col,'beta, intercept']=mod.params['const']
                res.loc[col,'P-Value, intercept']=mod.pvalues['const']
                x=mod.conf_int().loc['const']
                res.loc[col,'95 %tile, intercept']='['+str(round(x[0],3))+', '+str(round(x[1],3))+']'
        else:
            res.drop(col,inplace=True)
    
    for col in studyVars:
        ll=~res['P-Value, '+col].isna()
        res.loc[ll,'FDR, '+col]=multipletests(res.loc[ll,'P-Value, '+col],method='fdr_bh')[1]
    if showIntercept:
        res.loc[ll,'FDR, intercept']=multipletests(res.loc[ll,'P-Value, intercept'],method='fdr_bh')[1]
    return(res)

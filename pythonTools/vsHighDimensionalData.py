"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""

import numpy as np
import pandas as pd
try:
    from tqdm import tqdm
except Exception:
    def tqdm(iterable=None, total=None, desc=None):
        return iterable

# def imputeMissing(X):
#     ### EM-PCA / low-rank matrix completion
#     import numpy as np

#     w, Q = np.linalg.eigh(X.corr())
#     X0 = X.to_numpy()
#     mu = np.nanmean(X0, axis=0)
#     sd = np.nanstd(X0, axis=0, ddof=1)
    
#     # avoid divide-by-zero
#     sd[sd == 0] = 1.0
    
#     Z = (X0 - mu) / sd          # standardized data (so corr(Z)=R)
    
#     # sort eigenpairs descending
#     idx = np.argsort(w)[::-1]
#     w = w[idx]
#     Q = Q[:, idx]
    
#     k = 5                       # pick rank
#     Qk = Q[:, :k]
    
#     # scores (latent coordinates)
#     T = Z @ Qk                  # shape (n, k)
    
#     # reconstruct Z using k components
#     Zhat = T @ Qk.T             # (n, p)
    
#     # back to original scale
#     Xhat = Zhat * sd + mu
#     return(Xhat)


def imputeMissing(X, k=20, n_iter=20):
    """
    Impute missing values via iterative low-rank reconstruction.

    Parameters
    ----------
    X : pandas.DataFrame
        Data matrix with possible NaNs.
    k : int, default 20
        Rank for low-rank reconstruction.
    n_iter : int, default 20
        Number of reconstruction iterations.

    Returns
    -------
    pandas.DataFrame
        Imputed data with original index and columns.
    """
    X0 = X.to_numpy(dtype=float)
    mask = np.isnan(X0)  # True where missing

    # nan-aware standardization
    mu = np.nanmean(X0, axis=0)
    sd = np.nanstd(X0, axis=0, ddof=1)
    sd[~np.isfinite(sd)] = 1.0
    sd[sd == 0] = 1.0

    Z = (X0 - mu) / sd  # still has NaNs

    # initialize missing entries at 0 in standardized space (column mean)
    Z_imp = Z.copy()
    Z_imp[mask] = 0.0

    # build a consistent correlation matrix from the filled standardized data
    R = np.corrcoef(Z_imp, rowvar=False)
    R = 0.5 * (R + R.T)  # enforce symmetry

    # eigen-decomposition (more natural than svd for symmetric)
    w, Q = np.linalg.eigh(R)

    # sort descending
    idx = np.argsort(w)[::-1]
    w = w[idx]
    Q = Q[:, idx]

    Qk = Q[:, :k]

    # iterative low-rank reconstruction: only update missing entries
    for _ in range(n_iter):
        T = Z_imp @ Qk          # (n, k)
        Zhat = T @ Qk.T         # (n, p)
        Z_imp[mask] = Zhat[mask]

    # back to original scale
    Xhat = pd.DataFrame(Z_imp * sd + mu,index=X.index)
    Xhat.columns=list(X)
    
    return Xhat

def regressAll(vIndep,Y,studyVars=None,lbl=None,showIntercept=False):
    """
    Regress each column of Y on vIndep and summarize coefficients and p-values.

    Parameters
    ----------
    vIndep : pandas.DataFrame
        Independent variables (covariates).
    Y : pandas.DataFrame
        Dependent variables; each column is modeled separately.
    studyVars : list[str] | None
        Subset of vIndep columns to report FDR for; defaults to all columns.
    lbl : str | None
        Label for progress display.
    showIntercept : bool, default False
        Whether to include intercept estimates and FDR.

    Returns
    -------
    pandas.DataFrame
        Summary table with coefficients, p-values, confidence intervals, and FDR.
    """
    import statsmodels.api as sm
    from statsmodels.sandbox.stats.multicomp import multipletests
    
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
    for col in tqdm(res.index.values, total=N, desc=lbl or 'Regressing'):
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

def regressAllNA(vIndep,Y,studyVars=None,lbl=None,showIntercept=False):
    """
    Regress each column of Y with OLS and logistic models, aggregating p-values.

    Parameters
    ----------
    vIndep : pandas.DataFrame
        Independent variables (covariates).
    Y : pandas.DataFrame
        Dependent variables; each column is modeled separately.
    studyVars : list[str] | None
        Subset of vIndep columns to report aggregate FDR for.
    lbl : str | None
        Label for progress display.
    showIntercept : bool, default False
        Whether to include intercept estimates and FDR.

    Returns
    -------
    pandas.DataFrame
        Summary table with OLS/logistic results and aggregate p-values/FDR.
    """
    from vsPathways import computeGAGE
    import statsmodels.api as sm
    from statsmodels.sandbox.stats.multicomp import multipletests
    
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
    for prf in ['ols','logistic']:
        for v in studyVars:
            blank[prf+' beta, '+v]=np.nan
            blank[prf+' P-Value, '+v]=np.nan
            blank[prf+' 95 %tile, '+v]=''
        for v in [b for b in list(vIndep) if b not in studyVars]:
            blank[prf+' beta, '+v]=np.nan
            blank[prf+' P-Value, '+v]=np.nan
    res=pd.DataFrame(blank,index=list(Y))
    N=len(res)
    Z=vIndep.join(Y)
    X_Const=sm.add_constant(Z[list(vIndep)],has_constant='raise')
    for col in tqdm(res.index.values, total=N, desc=lbl or 'Regressing'):
        kp=~Z[col].isna()
        if X_Const.loc[kp,studyVars].std().max()>0:
            res.loc[col,'# Missing']=sum(~kp)
            mod=sm.OLS(Z.loc[kp,col],X_Const.loc[kp]).fit()
            res.loc[col,['ols beta, '+v for v in list(vIndep)]]=mod.params[list(vIndep)].values
            res.loc[col,['ols P-Value, '+v for v in list(vIndep)]]=mod.pvalues[list(vIndep)].values
            res.loc[col,['ols 95 %tile, '+v for v in studyVars]]=mod.conf_int().loc[studyVars].apply(lambda x:'['+str(round(x[0],3))+', '+str(round(x[1],3))+']',axis=1).values
            if showIntercept:
                res.loc[col,'ols beta, intercept']=mod.params['const']
                res.loc[col,'ols P-Value, intercept']=mod.pvalues['const']
                x=mod.conf_int().loc['const']
                res.loc[col,'ols 95 %tile, intercept']='['+str(round(x[0],3))+', '+str(round(x[1],3))+']'
            try:
                mod = sm.Logit(kp.astype(int), X_Const).fit()
                res.loc[col,['logistic beta, '+v for v in list(vIndep)]]=mod.params[list(vIndep)].values
                res.loc[col,['logistic P-Value, '+v for v in list(vIndep)]]=mod.pvalues[list(vIndep)].values
                res.loc[col,['logistic 95 %tile, '+v for v in studyVars]]=mod.conf_int().loc[studyVars].apply(lambda x:'['+str(round(x[0],3))+', '+str(round(x[1],3))+']',axis=1).values
                if showIntercept:
                    res.loc[col,'logistic beta, intercept']=mod.params['const']
                    res.loc[col,'logistic P-Value, intercept']=mod.pvalues['const']
                    x=mod.conf_int().loc['const']
                    res.loc[col,'logistic 95 %tile, intercept']='['+str(round(x[0],3))+', '+str(round(x[1],3))+']'
            except Exception as e:
                # optional: store the error so you can see why it failed
                res.loc[col, 'logistic error'] = str(e)
                

        else:
            res.drop(col,inplace=True)
    
    for v in studyVars:
        pvs=res[['ols P-Value, '+v,'logistic P-Value, '+v]]
        res['Aggregate P-Value, '+v]=pvs.apply(lambda x:computeGAGE(x.to_list()),axis=1)
        ll=~res['Aggregate P-Value, '+v].isna()
        res['Aggregate FDR, '+v]=multipletests(res.loc[ll,'Aggregate P-Value, '+v],method='fdr_bh')[1]
    if showIntercept:
        res.loc[ll,'FDR, intercept']=multipletests(res.loc[ll,'P-Value, intercept'],method='fdr_bh')[1]
    return(res)

def testIndexMatch(A,B,showAB=False,showBA=False):
    """
    Compare index membership between two DataFrames/Series and report differences.

    Parameters
    ----------
    A : pandas.DataFrame | pandas.Series
        First dataset with index to compare.
    B : pandas.DataFrame | pandas.Series
        Second dataset with index to compare.
    showAB : bool, default False
        If True, print elements in A not in B.
    showBA : bool, default False
        If True, print elements in B not in A.

    Returns
    -------
    None
        Prints summary to stdout.
    """
    a=pd.Series(A.index.values)
    b=pd.Series(B.index.values)
    for k,t in {'A':a,'B':b}.items():
        n=len(t)
        u=len(t.unique())
        print('Data set '+k+': '+str(n)+' rows with '+str(u)+' unique indices.')
        if u!=n:
            print(t.loc[t.duplicated(keep=False)].value_counts())
        if t.isna().sum()>0:
            print('Data set '+k+' has '+str(t.isna().sum())+' NaNs in its index')
            print('')
    
    ab=set(a.unique())-set(b.unique())
    ba=set(b.unique())-set(a.unique())
    union=set(b.unique()) & set(a.unique())
    print()
    print('Data set A has '+str(len(ab))+' unique elements that are not in B.')
    if showAB:
        print('Elements in A that are not in B:')
        print(', '.join(list(ab)))
    print()
    print('Data set B has '+str(len(ba))+' unique elements that are not in A.')
    if showBA:
        print('Elements in B that are not in A:')
        print(', '.join(list(ba)))
    print('There are '+str(len(union))+ ' unique elements in both indices.')
    
def v_buildSparse(x,y,delta=None):
    """
    Build a sparse contingency matrix from categorical indices.

    Parameters
    ----------
    x : pandas.Categorical
        Row categories.
    y : pandas.Categorical
        Column categories.
    delta : list[float] | None
        Values for each (x, y) pair; defaults to 1 for all.

    Returns
    -------
    pandas.DataFrame
        Sparse DataFrame with rows indexed by x categories and columns by y categories.
    """
    import scipy
    if delta is None:
        delta=[1]*len(x)
    i=x.codes
    j=y.codes
    rLbl=pd.Series(x.categories)
    cLbl=pd.Series(y.categories)
    M=scipy.sparse.coo_matrix((np.array(delta),(i,j)),shape=(len(rLbl),len(cLbl)))
    M=M.tocsr()
    
    return(pd.DataFrame.sparse.from_spmatrix(M,columns=cLbl,index=rLbl))

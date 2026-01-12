# Auto-generated from vsFeatures/ package modules.

# ---- source: vsFeatures/buildSparse.py ----
"""
Created on Wed Jun  7 10:29:12 2017

@author: jel2
"""

import numpy as np
import scipy

def buildSparse(df,a,b,colFilter=None,rowFilter=None):
    ct=df[[a,b]].groupby([a,b]).size()
    i=np.array(ct.index.labels[0])
    j=np.array(ct.index.labels[1])
    rLbl=list(ct.index.levels[0])
    cLbl=list(ct.index.levels[1])
    M=scipy.sparse.coo_matrix((np.array(ct),(i,j)),shape=(len(rLbl),len(cLbl)))
    M=M.tocsr()
    
    if colFilter:
        kp=np.asarray(M.sum(axis=0)).reshape(-1)>=colFilter
        M=M[:,kp]
        cLbl=list(np.asarray(cLbl)[kp])

    if rowFilter:
        kp=np.asarray(M.sum(axis=1)).reshape(-1)>=rowFilter
        M=M[kp,:]
        cLbl=list(np.asarray(rLbl)[kp])
    
    return([rLbl,cLbl,M])

# ---- source: vsFeatures/featureExtraction.py ----
"""
Created on Thu Mar  1 10:10:49 2018

@author: jel2
"""

def featureExtraction(eList,function, X=None, params=None):
    return X.join(function(eList,params),how='left').fillna(0)

# ---- source: vsFeatures/matrixWdw.py ----
"""
Created on Thu Mar  1 09:49:53 2018

@author: jel2
"""


def matrixWdw(eList,params):
    import pandas as pd
    
    if not params:
        params=dict()
    if 'wdw' not in params:
        params['wdw']=365.25
    if 'SOURCE' not in params:
        params['SOURCE']={'DIAGNOSIS'}
    if 'colName' not in params:
        params['colName']='LABEL'
    if 'rowName' not in params:
        params['rowName']='PATIENT_ID'
    if 'tmCol' not in params:
        params['tmCol']='delta'
    if 'minPt' not in params:
        params['minPt']=10
    
    ll=(eList[params['tmCol']]<0) & (eList[params['tmCol']]> -params['wdw']) & (eList.SOURCE.isin(params['SOURCE']))
    r,c,M=buildSparse(eList[ll],params['rowName'],params['colName'])
    M=pd.DataFrame(M.todense(),columns=c,index=r) ###M.todense() to convert sparse matrix to full matrix###
    M=M.loc[:,M.sum(axis=0)>=params['minPt']]
    return(M)

# ---- source: vsFeatures/sparseFactors.py ----
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

from vsVisualizations import loopProgress

def sparseFactors(Z,nFac=10,nSteps=500):
    import random
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    from scipy.stats import gamma,norm,binom,beta,multivariate_normal
    
    nSamp,nMeas=Z.shape
    
    
    phi0=1
    mu0=0
        
    u,s,v=np.linalg.svd(Z)
    nkp=100
    u=u*math.sqrt(nSamp)
    s=s/math.sqrt(nSamp)
    v=v[:nkp,:]*s[:nkp].reshape(-1,1)
    
    beta=v[:nFac,:]
    X=u[:,:nFac]
    q=np.array([.95]*nFac).reshape(-1,1)
    
    
    betaTrace=np.zeros((nFac,nMeas,nSteps))
    rTrace=np.zeros((nFac,nMeas,nSteps))
    XTrace=np.zeros((nFac,nSamp,nSteps))
    thetaTrace=np.zeros((nMeas,nSteps))
    for i in range(nSteps):
        if i%math.floor(nSteps/20)==0:
            loopProgress(math.floor(i*20/nSteps),20,'Step')
            
        ### Update theta
        epsilon=Z-X@beta
        theta=gamma.rvs(nSamp,scale=1/(epsilon**2).sum(axis=0))
        thetaTrace[:,i]=theta
        
        for j in range(nFac):
            ### Update beta and r
            cl=[v for v in range(nFac) if v!=j]
            Y=Z-X[:,cl]@beta[cl,:]
            
            nu=phi0+theta*np.square(X[:,[j]]).sum(axis=0)
            m=(mu0*phi0+theta*(Y*X[:,[j]]).sum(axis=0))/nu
        
            r=norm.pdf(0)/np.maximum(.000001,norm.pdf(0,loc=m,scale=1/np.sqrt(nu)))
            r=r/(q[j]+r)
    
            beta[j,:]=norm.rvs(loc=m,scale=1/np.sqrt(nu))
            betaTrace[j,:,i]=beta[j,:]
        
            zz=binom.rvs(1,r)
            beta[j,:]=beta[j,:]*zz
        
            rTrace[j,:,i]=r
        
            ### Update X - prior variance on X is always 1
            ll=zz!=0
            phi=(np.power(beta[j,ll],2)*theta[ll]).sum()
            m= (Y[:,ll]*(beta[j,ll]*theta[ll]).reshape(1,-1)).sum(axis=1)/phi
            X[:,j]=multivariate_normal.rvs(mean=m,cov=1/phi)
            XTrace[j,:,i]=X[:,j]
            

    return((XTrace,betaTrace,rTrace,thetaTrace))

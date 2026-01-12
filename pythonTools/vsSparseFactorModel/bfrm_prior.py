# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as st

from vsVisualizations.loopProgress import loopProgress

def pcImpute(M,nFac=2,mask=None,nSteps=1):
    if mask is None:
        mask=M.isna()
    for col in list(M):
        M.loc[M[col].isna(),col]=M[col].mean()
    for k in range(nSteps):
        U, S, VT = np.linalg.svd(M)
        M[mask]=pd.DataFrame(((U[:,:nFac]*S[:nFac])@VT[:nFac,:]),index=M.index,columns=list(M))[mask]
    return(M,mask)

def optimalBeta(Z,H,pThreshold=.05): 
    from sklearn.linear_model import Ridge
    rr=Ridge(alpha=.1)
    bList=list()
    for i in list(Z.index.values):
        t=sm.OLS(Z.loc[i],H).fit()
        bList+=[t.params*(t.pvalues<pThreshold)]
    return(pd.DataFrame(np.array(pd.concat(bList,axis=1)).T,index=Z.index,columns=list(H)))
    #     t=rr.fit(H,Z.loc[i])
    #     bList+=[t.coef_]
    # return(pd.DataFrame(np.stack(bList),index=Z.index,columns=list(H)))
    
def updateFactor(Y,bb,theta,h=None,theta0=None,q=None):
    if h is None:
        h=np.array([0]*Y.shape[1])
    else:
        h=np.array(h)
    if theta0 is None:
        theta0=np.array([.0001]*Y.shape[1])
    nu=theta.T@ np.array(bb**2)+theta0
    m=( (h*theta0)+ (np.array(Y)* (theta*np.array(bb).reshape(-1,1) )).sum(axis=0)  ) / nu
    res=st.norm.rvs(loc=m,scale=1/np.sqrt(nu))
    if q is not None:
        pp=updatePtMass(nu,m,q)
        return(res*pp,pp)
    else:
        return(res)

def updatePtMass(nu,m,q):
    r=st.norm.pdf(0)/np.maximum(.000001,st.norm.pdf(0,loc=m,scale=1/np.sqrt(nu)))
    return(r/(q+r))


def sampTheta(epsilon,pt_alpha=1,pt_beta=1):
    nSamp=epsilon.shape[1]
    tt=((epsilon**2).sum(axis=1)/2+pt_alpha)
    return(st.gamma.rvs(nSamp/2+pt_beta,scale=1/tt))


def sampZ(mu,sigma):    
    return(st.norm.rvs(loc=mu,scale=sigma))



def bfrm_prior(Z,H,z_mask=None,nFac=0,nSteps=500,burnin=0,hpStrength=50):
    import random
    import numpy as np
    import math
    import scipy.stats as st
    # from scipy.stats import gamma,norm,binom,beta,multivariate_normal
    from sklearn.decomposition import PCA
    import statsmodels.api as sm
    
    # ## add constant to H if it doesn't have one
    # if (H.std()==0).sum()==0:
    #     H=sm.add_constant(H)
    
    ## impute missing values in H
    H,h_mask=pcImpute(H,nSteps=10)
    
    ## Impute missing values in Z 
    if z_mask is not None:
        Z,z_mask=pcImpute(Z,mask=z_mask,nSteps=10)
    
    nMeas,nSamp=Z.shape
    
    beta=optimalBeta(Z,H)

    if nFac>0:
        ### Compute residuals
        R=Z-beta@H.T
        #### INITIALIZE R
        ### Instead of initializing with PCA, we will pass known variables.  These will be
        ### treated as priors rather than as fixed
        pca = PCA(n_components = nFac)
        R=pca.fit_transform(Z.T-H@beta.T)
        pcc=['pc'+str(nn).zfill(2) for nn in range(nFac)]
        H=pd.concat([H,pd.DataFrame(R,index=H.index,columns=pcc)],axis=1)
        h_mask=pd.concat([h_mask,pd.DataFrame(True,index=h_mask.index,columns=pcc)],axis=1)
        ### Re-estimate beta since we have added principal component columns to H
        beta=optimalBeta(Z,H)
        
    H0=H.copy()
    
    q=(beta!=0).mean()
    
    nVar=H.shape[1]
    betaTrace=np.zeros((nVar,nMeas,nSteps))
    pTrace=np.zeros((nVar,nMeas,nSteps))
    HTrace=np.zeros((nVar,nSamp,nSteps))
    thetaTrace=np.zeros((nMeas,nSteps))
    for i in range(nSteps+burnin):
        loopProgress(math.floor(i*20/(nSteps+burnin)),20,'Step')
        
        ### Update theta
        mn=beta@H.T
        
        theta=sampTheta(Z-mn)
        theta=np.tile(theta,(nSamp,1)).T
        
        ### sample the missing Z values
        if z_mask is not None:
            Z=pd.DataFrame(np.where(z_mask, sampZ(mn, 1/np.sqrt(theta)), Z),index=Z.index,columns=list(Z))
        
        if i>=burnin: thetaTrace[:,i-burnin]=theta[:,0]
        
        for cNum,j in enumerate(list(H)):
            ### Update beta and r
            cl=[v for v in list(H) if v!=j]
            Y=Z-beta[cl]@H[cl].T

            theta0=np.where(h_mask[j],.01,hpStrength)   # We use high precision if H was measured and low precision of it was missing
            H[j]=updateFactor(Y,beta[j],theta,h=H0[j],theta0=theta0)
            
            beta[j],p=updateFactor(Y.T,H[j],theta.T,q=q[j])
            q[j]=p.mean()
            
            if i>=burnin:
                betaTrace[cNum,:,i-burnin]=beta[j]
                HTrace[cNum,:,i-burnin]=H[j]
                pTrace[cNum,:,i-burnin]=p

    return((betaTrace, HTrace, thetaTrace, pTrace))

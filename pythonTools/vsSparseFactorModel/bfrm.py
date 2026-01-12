# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

import pandas as pd

from vsUtilityFunctions.standardize import standardize
from vsVisualizations.loopProgress import loopProgress

def bfrm(Z,nFac=None,F=None,nSteps=500,phi0=.01,mu0=0,burnin=0,updateF=True):
    import random
    import numpy as np
    import math
    import scipy.stats as st
    # from scipy.stats import gamma,norm,binom,beta,multivariate_normal
    from sklearn.decomposition import PCA
    import statsmodels.api as sm
    
    nSamp,nMeas=Z.shape
    if F is None:
        
        if type(Z)==np.ndarray:
            Z=pd.DataFrame(Z)
        Z=standardize(Z)
        if nFac is None:
            nFac=math.floor(math.sqrt(len(Z)))
    
        q0=list(np.arange(.8,.95,.15/nFac))[:nFac]
    
        #### INITIALIZE R
        ### Instead of initializing with PCA, we will pass known variables.  These will be
        ### treated as priors rather than as fixed
        pca = PCA(n_components = nFac)
        F=pca.fit_transform(Z)
    else:
        nFac=F.shape[1]
        q0=[.85]*nFac
    F=np.array(F)
    
    bList=list()
    ### THis is slow and has many redundant calculations
    for i in list(Z):
        t=sm.OLS(Z[i],F).fit()
        bList+=[t.params*(t.pvalues<.05)]
    # ### replace with this?    
    # Y = Z.values              # shape (n, m)
    # Ft = F.T
    # beta = np.linalg.solve(Ft @ F, Ft @ Y)   # shape (k, m)
    
    # # residuals
    # resid = Y - F @ beta
    # n, k = F.shape
    # sigma2 = (resid**2).sum(axis=0) / (n - k)
    
    # cov = np.linalg.inv(Ft @ F)
    # se = np.sqrt(np.diag(cov))[:, None]
    # tvals = beta / se
    # pvals = 2 * (1 - stats.t.cdf(np.abs(tvals), df=n-k))
    
    # bList = beta * (pvals < 0.05)

    beta=np.array(pd.concat(bList,axis=1))
    
    ss=F.std(axis=0)
    F=F/ss
    beta=beta*ss.reshape(-1,1)
    
    q=(beta==0).mean(axis=1)
    
    # alpha=(mu0+np.random.normal(size=(nFac,nMeas) ))/np.sqrt(phi0)
    # zz0=binom.rvs(1,q0,size=(nFac,nMeas))
    # beta0=alpha*zz0
    
    # F0=norm.rvs(size=(nSamp,nFac))
    
    # theta0=(gamma.rvs(4,size=nMeas)/2)
    
    # epsilon=np.random.normal(size=(nSamp,nMeas))*(1/np.sqrt(theta0))
    
    # Z=F0@beta0+epsilon
    
    # u,s,v=np.linalg.svd(Z)
    # nkp=100
    # u=u*math.sqrt(nSamp)
    # s=s/math.sqrt(nSamp)
    # v=v[:nkp,:]*s[:nkp].reshape(-1,1)
    
    # beta=v[:nFac,:]
    # F=u[:,:nFac]
    # q=np.array([.5]*nFac).reshape(-1,1)
    
    
    betaTrace=np.zeros((nFac,nMeas,nSteps))
    rTrace=np.zeros((nFac,nMeas,nSteps))
    FTrace=np.zeros((nFac,nSamp,nSteps))
    thetaTrace=np.zeros((nMeas,nSteps))
    for i in range(nSteps+burnin):
        if i%math.floor(nSteps/20)==0:
            loopProgress(math.floor(i*20/nSteps),20,'Step')
            
        ### Update theta
        epsilon=Z-F@beta
        theta=st.gamma.rvs(nSamp,scale=1/(epsilon**2).sum(axis=0))
        if i>=burnin: thetaTrace[:,i-burnin]=theta
        
        for j in range(nFac):
            ### Update beta and r
            cl=[v for v in range(nFac) if v!=j]
            Y=np.array(Z)-F[:,cl]@beta[cl,:]
            
            nu=phi0+theta*np.square(F[:,[j]]).sum(axis=0)
            m=(mu0*phi0+theta*(Y*F[:,[j]]).sum(axis=0))/nu
            
            r=st.norm.pdf(0)/np.maximum(.000001,st.norm.pdf(0,loc=m,scale=1/np.sqrt(nu)))
            r=r/(q0[j]+r)
    
            beta[j,:]=st.norm.rvs(loc=m,scale=1/np.sqrt(nu))
            zz=st.binom.rvs(1,r)
            beta[j,:]=beta[j,:]*zz
            
            ### Update F - prior variance on F is always 1
            ### Modify to include prior
            ll=zz!=0
            phi=(np.power(beta[j,ll],2)*theta[ll]).sum()
            m= (Y[:,ll]*(beta[j,ll]*theta[ll]).reshape(1,-1)).sum(axis=1)/phi
            if updateF:
                F[:,j]=st.multivariate_normal.rvs(mean=m,cov=1/phi)                
                ss=F[:,j].std()
                F[:,j]=F[:,j]/ss
            else:
                ss=1
            betaTrace[j,:,i]=betaTrace[j,:,i]*ss
            beta[j,:]=beta[j,:]*ss
            
            if i>=burnin:
                rTrace[j,:,i-burnin]=r
                betaTrace[j,:,i-burnin]=beta[j,:]
                FTrace[j,:,i-burnin]=F[:,j]
                
    return((betaTrace, rTrace, FTrace, thetaTrace))


# for i in range(nSteps):
#     if i%math.floor(nSteps/20)==0:
#         loopProgress(math.floor(i*20/nSteps),20,'Step')
        
#     ### Update theta
#     epsilon=Z-F@beta
#     theta=gamma.rvs(nSamp,scale=1/(epsilon**2).sum(axis=0))
#     thetaTrace[:,i]=theta
    
#     for j in range(nFac):
#         ### Update beta and r
#         cl=[v for v in range(nFac) if v!=j]
#         Y=Z-F[:,cl]@beta[cl,:]
        
#         nu=phi0+theta*np.square(F[:,[j]]).sum(axis=0)
#         m=(mu0*phi0+theta*(Y*F[:,[j]]).sum(axis=0))/nu
    
#         r=norm.pdf(0)/np.maximum(.000001,norm.pdf(0,loc=m,scale=1/np.sqrt(nu)))
#         r=r/(q[j]+r)

#         beta[j,:]=norm.rvs(loc=m,scale=1/np.sqrt(nu))
#         betaTrace[j,:,i]=beta[j,:]
    
#         zz=binom.rvs(1,r)
#         beta[j,:]=beta[j,:]*zz
    
#         rTrace[j,:,i]=r
    
#         ### Update F - prior variance on F is always 1
#         ll=zz!=0
#         phi=(np.power(beta[j,ll],2)*theta[ll]).sum()
#         m= (Y[:,ll]*(beta[j,ll]*theta[ll]).reshape(1,-1)).sum(axis=1)/phi
#         F[:,j]=multivariate_normal.rvs(mean=m,cov=1/phi)
#         FTrace[j,:,i]=F[:,j]

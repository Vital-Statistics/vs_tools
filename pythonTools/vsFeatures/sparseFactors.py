# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

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


# Auto-generated from vsSparseFactorModel/ package modules.

# ---- source: vsSparseFactorModel/bfrm.py ----
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

import pandas as pd

from vsUtilityFunctions import standardize
from vsVisualizations import loopProgress

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

# ---- source: vsSparseFactorModel/bfrm_prior.py ----
"""
Created on Sat Apr 23 23:28:10 2022

@author: joest
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
import scipy.stats as st

from vsVisualizations import loopProgress

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

# ---- source: vsSparseFactorModel/findFactorModes.py ----
#!/usr/bin/env python3
"""
Created on Sat Dec  6 10:57:14 2025

@author: rudy + ChatGPT
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def removeDups(omList, cos_tol=0.9):
    """
    Remove nearly-duplicate directions from omList based on absolute cosine similarity.
    Assumes all omegas are approximately unit vectors.
    It probably makes sense to adjust cos_tol based on the spread of the component
    """
    if len(omList) <= 1:
        return omList

    O = np.vstack(omList)              # shape (K, N)
    C = np.abs(O @ O.T)                # pairwise |cosine|
    rw, cl = np.where(C > cos_tol)

    # Keep the lower index, drop the higher one in each near-duplicate pair
    drp = {a for a, b in zip(rw, cl) if a > b}
    return [v for i, v in enumerate(omList) if i not in drp]


def initializePtmEstimate(X):
    """
    Initialize:
      - omList: unit vectors in the direction of each data point (row of X)
      - Xnorm2: squared norms of each row of X
    We no longer build full xi matrices here; they are implicit.
    """
    P, N = X.shape
    Xnorm2 = np.sum(X**2, axis=1)       # shape (P,)
    omList = [v / np.linalg.norm(v) for v in X]
    return omList, X, Xnorm2


def getWeights(omega, X, Xnorm2, beta, eps=1e-12):
    """
    Compute weights θ_i(ω, β) ∝ 1 / (ωᵀ M_i(β) ω), where
    M_i(β) = (2β + |X_i|^2) I - X_i X_iᵀ

    Using the scalar identity:
      ωᵀ M_i(β) ω = 2β + |X_i|^2 - (ωᵀ X_i)²
    """
    # proj_i = X_i ⋅ ω
    proj = X @ omega                    # shape (P,)
    quad = 2 * beta + Xnorm2 - proj**2  # ωᵀ M_i(β) ω for each i

    # Numerical safety
    quad = np.maximum(quad, eps)

    # θ_i ∝ (ωᵀ M_i ω)^(-1); alpha parameter in the gamma prior cancels in normalization
    w_raw = 1.0 / quad
    theta = w_raw / np.sum(w_raw)
    return theta


def build_M_from_theta(theta, X, Xnorm2, beta):
    """
    Build M(β, θ) = Σ_i θ_i M_i(β)
                   = 2β I + Σ_i θ_i (|X_i|² I - X_i X_iᵀ)
                   = (2β + Σ_i θ_i |X_i|²) I - Σ_i θ_i X_i X_iᵀ
    """
    P, N = X.shape

    # Σ_i θ_i |X_i|²
    s = np.dot(theta, Xnorm2)

    # Σ_i θ_i X_i X_iᵀ = (X_weighted)ᵀ (X_weighted) with sqrt(theta_i)
    Xw = np.sqrt(theta)[:, None] * X    # shape (P, N)
    S = Xw.T @ Xw                       # N x N

    M = (2 * beta + s) * np.eye(N) - S
    return M


def updatePtmEstimate(beta, omList, X, Xnorm2):
    """
    One EM-like update sweep over all current modes, at a fixed β.

    For each omega_k:
      1. Compute θ_i(omega_k).
      2. Build M_k = Σ_i θ_i M_i(β).
      3. Update omega_k to the eigenvector of M_k with smallest eigenvalue.
    """
    _, N = X.shape
    # alpha = (N - 1) / 2 + alpha0  # alpha does not affect the location of the modes

    newOmList = []

    for omega in omList:
        omega0 = omega.copy()

        # weights based on current direction
        theta = getWeights(omega, X, Xnorm2, beta)

        # weighted matrix and smallest-eigenvector update
        M = build_M_from_theta(theta, X, Xnorm2, beta)
        eigenvalues, eigenvectors = np.linalg.eigh(M)
        omega_new = eigenvectors[:, np.argmin(eigenvalues)]

        # fix sign for continuity
        omega_new *= np.sign(omega_new @ omega0)

        newOmList.append(omega_new)

    return newOmList

def findModes(X,showPlots=False,cos_tol=.999):
    omList, X_used, Xnorm2 = initializePtmEstimate(X)

    ################### variable step size
    beta = 1e-3
    beta_max = 1.0
    d_beta = 1e-4
    min_dBeta=1e-4

    omTrace=list()
    while beta < beta_max:
        # one EM sweep at current beta
        nn=len(omList)
        omList = removeDups(omList, cos_tol=cos_tol)
        if len(omList)<nn:
            omTrace+=[(beta,omList)]
        if nn==1:
            break
        prev_omList = [o.copy() for o in omList]
        omList = updatePtmEstimate(beta, omList, X, Xnorm2)

        # # measure movement
        # max_move = max([np.linalg.norm(o - p) for o, p in zip(omList, prev_omList)])

        # # measure min distance between distinct modes (for merging resolution)
        # if len(omList) > 1:
        #     O = np.vstack(omList)
        #     cos = np.abs(O @ O.T)
        #     np.fill_diagonal(cos, 0.0)
        #     min_dist = 1 - np.max(cos)  # 1 - cos(angle)
        # else:
        #     min_dist = 1.0

        # # adapt step
        # if min_dist > 0.1 and max_move < 1e-3:
        #     d_beta *= 1.5   # modes well-separated and stable → jump faster
        # elif min_dist < 0.02:
        #     d_beta *= 0.5   # approaching merge → refine β

        # d_beta = max(min(d_beta, beta_max - beta),min_dBeta)
        beta += d_beta

        if showPlots:
            plt.scatter(X[:, 0], X[:, 1], alpha=0.3)
            for i, omega in enumerate(omList):
                plt.plot([0, omega[0]], [0, omega[1]], color="gray")
                plt.text(omega[0], omega[1], str(i), fontsize=8)
            plt.title(f"beta={beta:.4f}, num modes={len(omList)}")
            plt.show()
            plt.clf()
            plt.close()
            
    return(omTrace)

def fltTbl(tbl):
    ta=tbl.copy()
    ta['mm']=cSlope(ta)
    tb=ta.loc[ta.mm<0].copy()
    while(set(ta.index.values)!=set(tb.index.values)):
        ta=tb
        ta['mm']=cSlope(ta)
        tb=ta.loc[ta.mm<0].copy()
    return(ta)

def cSlope(tbl):
    dN=[a-b for a, b in zip(tbl.N.iloc[1:],tbl.N.iloc[:-1])]
    dd=[a-b for a, b in zip(tbl.delta.iloc[1:],tbl.delta.iloc[:-1])]
    return([-np.inf]+[a/b if b>0 else np.inf for a,b in zip(dN,dd)])

def filterModes(omTrace,showPlots=False):
    beta=[a for a,b in omTrace]
    tbl=pd.DataFrame({'delta':[a-b for a,b in zip(beta[1:],beta[:-1])]})
    tbl['N']=[len(b) for a,b in omTrace[:-1]]
    tbl=tbl.loc[(tbl.N<tbl.N.max()/4) & (tbl.N>1)]

    tbl=tbl.sort_values('N',ascending=False).reset_index(drop=True)
    tbl['mm']=cSlope(tbl)

    tt=fltTbl(tbl)
    if showPlots:
        plt.scatter(tbl.delta,tbl.N)
        plt.scatter(tt.delta,tt.N)
    kpn=set(tt.N.unique())
    return([(beta,V) for beta,V in omTrace if len(V) in kpn])

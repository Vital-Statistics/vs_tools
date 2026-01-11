#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 09:51:29 2025

@author: rudy
"""

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
    import numpy as np
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

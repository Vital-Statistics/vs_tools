#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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



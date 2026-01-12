# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 10:29:12 2017

@author: jel2
"""

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

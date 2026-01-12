# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 09:25:16 2018

@author: jel2
"""

import scipy.stats as stat

def multiTest(q,M,f=stat.fisher_exact,lbl='feature',cLbl=None,verbose=False):
    if not cLbl:
        cLbl=['Feature '+str(x) for x in list(range(0,M.shape[1]))]
    if f not in {stat.fisher_exact,stat.chi2_contingency}:
        print('Error: Only scipy.stats functions fisher_exact and chi2_contingency are supported')
        return
    res=pd.DataFrame({lbl:cLbl, 'p-value':np.ones(len(cLbl)), 'odds ratio':np.zeros(len(cLbl)),'N':np.zeros(len(cLbl)),'PPV':np.zeros(len(cLbl)),'NPV':np.zeros(len(cLbl)),'Sensitivity':np.zeros(len(cLbl)),'Specificity':np.zeros(len(cLbl))})
    
    #### this test could be done in the loop to avoid copying code
    if(scipy.sparse.issparse(M)):
        for i in range(0,M.shape[1]):
            t=np.asarray(pd.crosstab(q,np.asarray(M[:,i].todense()).reshape(-1)>0))
            res.loc[i,'p-value']=f(t)[1]
            res.loc[i,'odds ratio']=t[0,0]*t[1,1]/(t[1,0]*t[0,1])
            res.loc[i,'N']=sum(t[:,1])
            res.loc[i,'PPV']=t[1,1]/sum(t[:,1])
            res.loc[i,'NPV']=t[0,0]/sum(t[:,0])
            res.loc[i,'Sensitivity']=t[1,1]/sum(t[1,:])
            res.loc[i,'Specificity']=t[0,0]/sum(t[0,:])
#            if i%100==0 & verbose:
#                print(str(i)+' of '+str(len(cLbl))+'\n')
    else:
        for i in range(0,M.shape[1]):
            t=np.asarray(pd.crosstab(q,np.asarray(M[:,i]).reshape(-1)>0))
            res.loc[i,'p-value']=f(t)[1]
            res.loc[i,'odds ratio']=t[0,0]*t[1,1]/(t[1,0]*t[0,1])
            res.loc[i,'N']=sum(t[:,1])
            res.loc[i,'PPV']=t[1,1]/sum(t[:,1])
            res.loc[i,'NPV']=t[0,0]/sum(t[:,0])
            res.loc[i,'Sensitivity']=t[1,1]/sum(t[1,:])
            res.loc[i,'Specificity']=t[0,0]/sum(t[0,:])
#            if i%100==0 & verbose:
#                print(str(i)+' of '+str(len(cLbl))+'\n')
    
    return(res)
    
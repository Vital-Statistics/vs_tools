# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 13:29:53 2018

@author: jel2
"""

def xval(x,y,g,mdl,params):
    M=pd.DataFrame(g).join(y,how='outer').join(x,how='outer')
#    Y=np.array(M[list(y)])
    #### if Y is two columns with one binary column then treat it as survival data
#    if(len(Y.shape)>1):
#        surv=True
#        ix=np.where(np.apply_along_axis(lambda z: len(set(z)),0,Y)==2)[0].item(0)
#        tx=({0,1}-{ix}).pop()
#        ix=Y[:,ix]
#        tm=Y[:,tx]
#    else:
#        surv=False
#    G=np.array(M[list(g)]).transpose()[0]
#    X=np.array(M[list(x)])
    M['xVal']=0
    for i in set(M.grp):
        ll= M.grp==i
        mod=lsReg(M.loc[~ll,:].drop(['grp','xVal'],axis=1),params)
        M.loc[ll,'xVal']=mod.predict_expectation(M.loc[ll,:].drop(['grp','xVal'],axis=1)).iloc[:,0]
        
    return(M)
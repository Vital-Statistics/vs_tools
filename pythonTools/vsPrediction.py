# Auto-generated from vsPrediction/ package modules.

# ---- source: vsPrediction/lsReg.py ----
"""
Created on Thu Mar  1 14:14:51 2018

@author: jel2
"""

def lsReg(M,params):
    import lifelines as lf
    
    cph = lf.CoxPHFitter()
    return(cph.fit(M, duration_col=params['tmCol'], event_col=params['ixCol']))
        

# ---- source: vsPrediction/predict_proba.py ----
#!/usr/bin/env python3
"""
Created on Sat Sep  9 12:46:51 2023

@author: rudy
"""

import numpy as np
import pandas as pd

def predict_proba(mod,Z):
    ll=(Z.isna().sum(axis=1)==0)
    yhat=pd.Series(np.nan,index=Z.index)
    yhat[ll]=mod.predict_proba(Z.loc[ll])[:,1]
    return(yhat)

# ---- source: vsPrediction/selectAndTrain.py ----
#!/usr/bin/env python3
"""
Created on Sat Sep  9 09:07:54 2023

@author: rudy
"""

import numpy as np
import pandas as pd

from vsVisualizations import loopProgress

def selectAndTrain(X,y,lbl='Setting regression penalty',returnYhat=False,verbose=False,leaveOneOut=False):
    from sklearn.metrics import roc_auc_score
    from scipy import stats as stat
    from scipy.special import logit
    from statsmodels.sandbox.stats.multicomp import multipletests
    from sklearn.linear_model import LogisticRegression

    ### Set Regression Penalty with cross validation
    nNan=pd.concat([X,y],axis=1).isna().sum(axis=1)
    nDrop=(nNan>0).sum()
    if nDrop>0:
        print('Dropping '+str(nDrop)+' observations due to missingness.')
        X=X.loc[nNan==0]
        y=y.loc[nNan==0]
    xvp=pd.DataFrame({'tt':np.arange(.025,2.0,.025),'cstat':.5,'P-Value':1.0,'FDR':1.0})
    for i,tt in enumerate(xvp.tt):
        loopProgress(i,len(xvp),lbl)
        mod=LogisticRegression(max_iter=1000,penalty='l1',solver='liblinear',C=tt)
        yhat=xval_train(mod,X,y,leaveOneOut=leaveOneOut)
        
        xvp.loc[xvp.tt==tt,'P-Value']=stat.ranksums(yhat.loc[y==0],yhat.loc[y==1],alternative='less')[1]  ## One sided test since the model is designed to assign higher values when outcome is 1
        
        xvp.loc[xvp.tt==tt,'cstat']=roc_auc_score(y,yhat)
        if xvp.loc[i,'cstat']==xvp.cstat.max():
            bestYHat=yhat
    xvp.sort_values('cstat',inplace=True,ascending=False)
    
    ### Retrain the model using the identified regression penalty
    xvp['FDR']=multipletests(xvp['P-Value'],method='fdr_bh')[1]
    tt=xvp.iloc[0].tt
    
    mod=LogisticRegression(max_iter=1000,penalty='l1',solver='liblinear',C=tt)
    mod.fit(X,y)
    
    ## identify the variables that were used
    aUsd=pd.DataFrame({'Analyte':list(X),'beta':mod.coef_.squeeze()})
    kpZ=list(aUsd.loc[aUsd.beta!=0,'Analyte'])
    
    if verbose:
        print(xvp.head(10))
    
    if returnYhat:
        res=(mod,kpZ,bestYHat)
    else:
        res=(mod,kpZ)

    ########## Accuracy on validation data
    return(res)

# We used logistic regression with a lasso penalty to simultaneously perform variable selection
# and model fitting.  We tested a range of penalties, in each case testing accuracy on hold outs
# using leave-one-out cross-validation.  All performance statistics are on the highest performing
# model, but are computed on the held out samples only.  Variables identified as included in the
# model are those that had non-zero regression coefficients when using the optimal penalty on the
# full data set.

# ---- source: vsPrediction/xval.py ----
"""
Created on Thu Mar  1 13:29:53 2018

@author: jel2
"""

import pandas as pd


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

# ---- source: vsPrediction/xval_train.py ----
"""
Created on Tue May 26 13:52:01 2020

@author: jolucas
"""

import pandas as pd

def xval_train(mod,X,Y,grp=None,folds=20,leaveOneOut=False,verbose=False):
    import random
    import numpy as np
    
    # X=np.array(X)
    # Y=np.array(Y)
    if type(X)==np.ndarray:
        X=pd.DataFrame(X,index=range(len(X)),columns=['X'+v for v in range(X.shape[1])])
        Y=pd.DataFrame(Y,index=range(len(Y)),columns=['Y'+v for v in range(Y.shape[1])])
    
    S=pd.concat([Y.rename('y'),X],axis=1)
    S['xval']=0.0
    cols=list(X)
    if grp is None:
        if not leaveOneOut:
            S['grp']=pd.Series(random.choices(list(range(folds)),k=len(S)),index=S.index)
        else:
            S['grp']=range(len(S))
    else:
        S['grp']=grp
    for i,g in enumerate(S.grp.unique()):
        if verbose:
            if i==0:
                print('Step '+str(i)+' of '+str(S.grp.unique().size),end='')
            else:
                print('\rStep '+str(i)+' of '+str(S.grp.unique().size),end='')
        ll=S.grp!=g
        mod.fit(S.loc[ll,cols],S.loc[ll,'y'])
        S.loc[~ll,'xval']=mod.predict_proba(S.loc[~ll,cols])[:,1]
    if verbose: print('\n')
    return(S['xval'])

def xval_train_xgb(mod,X,y,folds=20,verbose=False):
    import random
    
    S=pd.concat([y.rename('y'),X],axis=1)
    S['xval']=0
    cols=list(X)
    grp=pd.Series(random.choices(list(range(folds)),k=len(y)),index=y.index)
    yy=pd.DataFrame({'y':list(y),'grp':grp})
    for i in range(folds):
        if verbose:
            print(str(i)+' of '+str(folds))
        ll=yy.grp!=i
        A=S.loc[ll,cols]
        mu=A.mean()
        sig=A.std()
        # mod.fit(S.loc[ll,cols],S.loc[ll,'y'])
        mod.fit((A-mu)/sig,S.loc[ll,'y'])
        # betaTrace.append(mod.coef_[0])
        S.loc[~ll,'xval']=mod.predict((S.loc[~ll,cols]-mu)/sig)
    return(S['xval'])

def xval_train_penalty(mod,X,y):
    import random
    
    S=pd.merge(y.rename('y'),X,left_index=True,right_index=True)
    S['xval']=0
    cols=list(X)
    for i in S.index.values:
        ll=S.index.values!=i
        A=S.loc[ll,cols]
        mu=A.mean()
        sig=A.std()
        # mod.fit(S.loc[ll,cols],S.loc[ll,'y'])
        mod.fit((A-mu)/sig,S.loc[ll,'y'])
        # S.loc[~ll,'xval']=mod.predict_proba((S.loc[~ll,cols]-mu)/sig)[:,1]
        S.loc[~ll,'xval']=mod.predict((S.loc[~ll,cols]-mu)/sig)[0]
    return(S['xval'])

def xval_train_deep(mod,X,y,ep=150,bs=50):
    import random
    
    S=pd.merge(y.rename('y'),X,left_index=True,right_index=True)
    S['xval']=0
    cols=list(X)
    ww = mod.get_weights()
    for i in S.index.values:
        ll=S.index.values!=i
        mod.set_weights(ww)
        mod.fit(S.loc[ll,cols],S.loc[ll,'y'],epochs=ep,batch_size=bs,verbose=0)
        S.loc[~ll,'xval']=mod.predict_proba(S.loc[~ll,cols])
    return(S['xval'])

# ---- source: vsPrediction/xval_train_Other.py ----
"""
Created on Tue May 26 13:52:01 2020

@author: jolucas
"""

import pandas as pd

def xval_trainV2(mod,X,y):
    S=pd.merge(y.rename('y'),X,left_index=True,right_index=True)
    S['xval']=0
    cols=list(X)
    betaTrace=[]
    for i in S.index.values:
        ll=S.index.values!=i
        A=S.loc[ll,cols]
        mu=A.mean()
        sig=A.std()
        # mod.fit(S.loc[ll,cols],S.loc[ll,'y'])
        mod.fit((A-mu)/sig,S.loc[ll,'y'])
        betaTrace.append(mod.coef_[0])
        S.loc[~ll,'xval']=mod.predict_proba((S.loc[~ll,cols]-mu)/sig)[:,1]
    return(S['xval'],pd.DataFrame(dict(zip(cols,zip(*betaTrace)))))

def xval_train_penalty(mod,X,y):
    S=pd.merge(y.rename('y'),X,left_index=True,right_index=True)
    S['xval']=0
    cols=list(X)
    for i in S.index.values:
        ll=S.index.values!=i
        A=S.loc[ll,cols]
        mu=A.mean()
        sig=A.std()
        # mod.fit(S.loc[ll,cols],S.loc[ll,'y'])
        mod.fit((A-mu)/sig,S.loc[ll,'y'])
        # S.loc[~ll,'xval']=mod.predict_proba((S.loc[~ll,cols]-mu)/sig)[:,1]
        S.loc[~ll,'xval']=mod.predict((S.loc[~ll,cols]-mu)/sig)[0]
    return(S['xval'])

def xval_train_deep(mod,X,y,ep=150,bs=50):
    S=pd.merge(y.rename('y'),X,left_index=True,right_index=True)
    S['xval']=0
    cols=list(X)
    ww = mod.get_weights()
    for i in S.index.values:
        ll=S.index.values!=i
        mod.set_weights(ww)
        mod.fit(S.loc[ll,cols],S.loc[ll,'y'],epochs=ep,batch_size=bs,verbose=0)
        S.loc[~ll,'xval']=mod.predict_proba(S.loc[~ll,cols])
    return(S['xval'])

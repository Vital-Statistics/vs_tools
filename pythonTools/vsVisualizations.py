# Auto-generated from vsVisualizations/ package modules.

# ---- source: vsVisualizations/boxSwarm.py ----
#!/usr/bin/env python3
"""
Created on Sun Jan  4 11:51:49 2026

@author: rudy
"""

import pandas as pd


def boxSwarm(group,value):
    import seaborn as sns
    df=pd.DataFrame({'value':value,'group':group})

    ax=sns.boxplot(data=df, x="group", y="value",hue="group",showfliers=False, width=0.6,palette=vs_palette())
    sns.swarmplot(data=df, x="group", y="value",hue="group", size=4, alpha=0.7,palette=vs_palette())
    # ax.legend_.remove()

# ---- source: vsVisualizations/histogram.py ----
"""
Created on Fri Mar  9 08:07:47 2018

@author: jel2
"""

def stratifiedHistogram(t,x,group=None,bins=10):
    
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    y=t[x]
    if(group is None):
        grp=pd.Series('Population',index=t.index)
    else:
        grp=t[group]
    
    b=np.arange(min(y),max(y),(max(y)-min(y))/bins)
    clr=[(1,0,0,.3),(0,1,0,.3),(0,0,1,.3)]
    for i,g in enumerate(set(grp)):
        plt.hist(y[grp==g],bins=b,label=g,normed=True,color=clr[i])
        
    plt.legend()
    plt.xlabel(x)
    plt.ylabel('Probability')
#    plt.title('Age distribution')

def temporalHistogram(t):
    import matplotlib.pyplot as plt
    import pandas as pd
    #### t=t.astype("datetime64")
    ### t should be a series

    m=(t-pd.to_timedelta(t.dt.dayofweek, unit='d')).value_counts()
    plt.scatter(m.index.values,m)
    plt.xlabel('Time')
    plt.ylabel('Count')
    

# ---- source: vsVisualizations/loopProgress.py ----
"""
Created on Mon Dec 20 12:16:57 2021

@author: JoeLucas
"""

def loopProgress(i,N,lbl='',nSteps=20):
    import math
    i=i+1
    if i==1:
        print(lbl+': |'+'-'*i+' '*(nSteps-i)+'|',end='')
    i=math.floor(i*nSteps/N)
    N=nSteps
    if i<N:
        print('\r'+lbl+': |'+'-'*i+' '*(N-i)+'|',end='')
    else:
        print('\r'+lbl+': |'+'-'*i+' '*(N-i)+'|')
    
    

# ---- source: vsVisualizations/multiCorrPlot.py ----
"""
Created on Fri Jun  3 10:08:42 2022

@author: joest
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compSig(dat):
    dat=dat.loc[dat.isna().sum(axis=1)==0]
    dat=np.array(dat)
    mn=dat.mean(axis=0)
    dat=dat-mn
    res=np.array(dat.transpose())@np.array(dat)/len(dat)
    return(res)

def mvnOval(a,b,color=None,scale=2.65,alpha=.2):
    from scipy.linalg import sqrtm
    t=pd.concat([a,b],axis=1)
    mn=np.array(t).mean(axis=0)
    sig=scale*compSig(t)
    x=np.array([[math.cos(theta),math.sin(theta)] for theta in np.arange(0,2*math.pi,2*math.pi/40)])
    res=x@sqrtm(sig+1e-6*np.identity(2))
    if color is not None:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],color=color,edgecolor=None,alpha=alpha)
    else:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],edgecolor=None,alpha=alpha)
    return res

def multiCorrPlot(u,x,y,grp,xlbl='',ylbl=''):
    u=u.loc[u[[x,y]].isna().sum(axis=1)==0]
    nclr=len(u[grp].unique())-1
    for i,gg in enumerate(u[grp].unique()):
        plt.scatter(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],label=gg,color=vs_clr(nclr-i))
        mvnOval(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],color=vs_clr(nclr-i,alpha=.2))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend()
    plt.tight_layout()


# ---- source: vsVisualizations/pivotalEncounter.py ----
"""
Created on Mon Mar 12 10:04:12 2018

@author: jel2
"""

import datetime
import pandas as pd

def pivotalEncounter(evt,M,colName='LABEL',wdw=365.25/2,buffer=1,flt=None,countPatients=True):
    ############## Examples of diagnoses that occur more often after the initial diagnosis than before
    import math
    import scipy.stats as stat
    
    M=M.join(evt,'PATIENT_ID')

    q=datetime.timedelta(days=buffer)
    w=datetime.timedelta(days=wdw)
    afterLbl='Subsequent '+str(math.floor(wdw))+' days'
    beforeLbl='Prior '+str(math.floor(wdw))+' days'
    durLbl='During encounter'
    
    kp=(M['ENCOUNTER_ID']!=M['eID']) & ((M['EVENT_TIME']<M['eStartTime']-q) & (M['EVENT_TIME']>M['eStartTime']-w))
    u=M.loc[kp,['PATIENT_ID',colName]]
    if countPatients:
        u=u.drop_duplicates()
    v=pd.DataFrame(u[colName].value_counts().rename(beforeLbl))
    
    kp=(M['ENCOUNTER_ID']==M['eID']) | ((M['EVENT_TIME']>=M['eStartTime']-q) & (M['EVENT_TIME']<=M['eStartTime']+q))
    u=M.loc[kp,['PATIENT_ID',colName]]
    if countPatients:
        u=u.drop_duplicates()
    v=pd.merge(v,pd.DataFrame(u[colName].value_counts().rename(durLbl)),how='outer',left_index=True,right_index=True)
    
    kp=(M['ENCOUNTER_ID']!=M['eID']) & ((M['EVENT_TIME']>M['eEndTime']+q) & (M['EVENT_TIME']<M['eEndTime']+w))
    u=M.loc[kp,['PATIENT_ID',colName]]
    if countPatients:
        u=u.drop_duplicates()
    v=pd.merge(v,pd.DataFrame(u[colName].value_counts().rename(afterLbl)),how='outer',left_index=True,right_index=True)
    
    v=v.fillna(0)    

    if flt is None:
        flt=min(.005*len(set(M['PATIENT_ID'])),10)

    v=v.loc[v.sum(axis=1)>flt,:]
    n=sum(v[beforeLbl])
    m=sum(v[afterLbl])
    
    v['OR']=1
    v['p value']=1
    
    for i in v.index.values:
        v.loc[i,'OR'],v.loc[i,'p value']=stat.fisher_exact([[v.loc[i,afterLbl],v.loc[i,beforeLbl]],[m-v.loc[i,afterLbl],n-v.loc[i,beforeLbl]]])
    
#    priorWeight=6
#    priorMean=.5
#    v['ratio']=(v[afterLbl]+priorWeight*priorMean)/(v[beforeLbl]+priorWeight*(1-priorMean))
#    v['p value']=v.apply(lambda x: stat.beta.cdf(.5,x[afterLbl]+priorWeight*priorMean,x[beforeLbl]+priorWeight*(1-priorMean)),axis=1)
#    v['p value']=v['p value'].apply(lambda x: min(x,1-x)*2)
    v=v.sort_values('p value')
    
    return(v)
    

# ---- source: vsVisualizations/stratifiedHistogram.py ----
"""
Created on Fri Mar  9 08:07:47 2018

@author: jel2
"""

def stratifiedHistogram(t,x,group=None,bins=10):
    
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    y=t[x]
    if(group is None):
        grp=pd.Series('Population',index=t.index)
    else:
        grp=t[group]
    
    b=np.arange(min(y),max(y),(max(y)-min(y))/bins)
    clr=[(1,0,0,.3),(0,1,0,.3),(0,0,1,.3)]
    for i,g in enumerate(set(grp)):
        plt.hist(y[grp==g],bins=b,label=g,normed=True,color=clr[i])
        
    plt.legend()
    plt.xlabel(x)
    plt.ylabel('Probability')
#    plt.title('Age distribution')



# ---- source: vsVisualizations/stratifiedSurvival.py ----
"""
Created on Fri Mar  9 08:07:47 2018

@author: jel2
"""

def stratifiedSurvival(t,eventTime,eventIndicator=None,followupTime=None,group=None):
    
    import matplotlib.pyplot as plt
    import lifelines as lf
    from lifelines.plotting import add_at_risk_counts
    import pandas as pd
    import copy

    tm=t[eventTime].copy()
    
    if(group is None):
        grp=pd.Series('Population',index=t.index)
    else:
        grp=t[group]
    
    if(eventIndicator is None):
        ev=~t[eventTime].isnull()
        tm[tm.isnull()]=t.loc[tm.isnull(),followupTime]
        
    ######### Kaplan Meier curves stratified by sex
    kl=list()
    kmf = lf.KaplanMeierFitter()
    fig,ax=plt.subplots()
    for g in set(grp):
        kmf.fit(tm[grp==g],ev[grp==g],label=g)
        kmf.plot(ax=ax)
        kl.append(copy.deepcopy(kmf))
        
    add_at_risk_counts(*kl, ax=ax)
        
    plt.legend(loc='lower left')
    plt.ylim([0,1])
    plt.xlabel('Time (years)')
    plt.ylabel('Survival')
    plt.title('Kaplan-Meier survival curve')
#    add_at_risk_counts(kmf1,kmf2, ax=ax)


#ax.spines['bottom'].set_position(('axes', -0.15 * 6.0 / fig.get_figheight()))

# ---- source: vsVisualizations/twoCorrPlot.py ----
"""
Created on Fri Jun  3 10:08:42 2022

@author: joest
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def compSig(dat):
    dat=dat.loc[dat.isna().sum(axis=1)==0]
    dat=np.array(dat)
    mn=dat.mean(axis=0)
    dat=dat-mn
    res=np.array(dat.transpose())@np.array(dat)/len(dat)
    return(res)

def mvnOval(a,b,color=None,scale=2.65,alpha=.2):
    from scipy.linalg import sqrtm
    t=pd.concat([a,b],axis=1)
    mn=np.array(t).mean(axis=0)
    sig=scale*compSig(t)
    x=np.array([[math.cos(theta),math.sin(theta)] for theta in np.arange(0,2*math.pi,2*math.pi/40)])
    res=x@sqrtm(sig+1e-6*np.identity(2))
    if color is not None:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],color=color,edgecolor=None,alpha=alpha)
    else:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],edgecolor=None,alpha=alpha)
    return res

def twoCorrPlot(x0,y0,x1,y1,lbl0='',lbl1='',xlbl='',ylbl=''):
    plt.scatter(x0,y0,color=vs_clr('dark gray'),label=lbl0)
    mvnOval(x0,y0,color=vs_clr('dark gray',alpha=.2))
    plt.scatter(x1,y1,color=vs_clr(),label=lbl1)
    mvnOval(x1,y1,color=vs_clr(alpha=.2))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend()
    plt.tight_layout()

def multiCorrPlot(u,x,y,grp,xlbl='',ylbl=''):
    u=u.loc[u[[x,y]].isna().sum(axis=1)==0]
    nclr=len(u[grp].unique())-1
    for i,gg in enumerate(u[grp].unique()):
        plt.scatter(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],label=gg,color=vs_clr(nclr-i))
        mvnOval(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],color=vs_clr(nclr-i,alpha=.2))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend()
    plt.tight_layout()


# ---- source: vsVisualizations/vs_clr.py ----
"""
Created on Mon Feb  1 10:47:06 2021

@author: JoeLucas
"""

def vs_clr(cc='blue',alpha=1):
    r={'blue':(15/256,98/256,160/256)
       ,'dark blue':(48/256,48/256,95/256)
       ,'grey':(182/256,183/256,182/256)
       ,'gray':(182/256,183/256,182/256)
       ,'dark grey':(84/256,85/256,84/256)
       ,'dark gray':(84/256,85/256,84/256)
       }
    if type(cc)==int:
        cc=['blue','dark gray','gray'][cc]
    res=r[cc]
    if alpha<1:
        res=tuple(list(res)+[alpha])
    return(res)

def vs_clr_V1(cc='blue',alpha=1):
    r=(15/256,98/256,160/256)
    if(cc=='dark blue'):
        r=(48/256,48/256,95/256)
    elif((cc=='grey') | (cc=='gray')):
        r=(182/256,183/256,182/256)
    elif((cc=='dark grey') | (cc=='dark gray')):
        r=(84/256,85/256,84/256)
    if alpha<1:
        r=tuple(list(r)+[alpha])
    return(r)

# ---- source: vsVisualizations/vs_palette.py ----
#!/usr/bin/env python3
"""
Created on Thu Feb 15 15:52:55 2024

@author: rudy
"""

def vs_palette():
    pp=['#0e619f',
     '#b22222',
     '#006400',
     '#b5b6b5',
     '#2f2f5e',
     '#e57373',
     '#98fb98',
     '#535453']
    return(pp)

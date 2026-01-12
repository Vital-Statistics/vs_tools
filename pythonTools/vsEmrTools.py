# Auto-generated from vsEmrTools/ package modules.

# ---- source: vsEmrTools/emrTools.py ----
"""
Created on Wed May  3 15:18:36 2017

@author: jel2
"""

import datetime
from collections import Counter
import numpy as np
import pandas as pd
import scipy
import operator
from numpy import ones

def allValsList(t,pth,res=None):
    if(res==None):
        res=list()
    if (type(t)==list) | (type(t)==set):
        for s in t:
            allValsList(s,pth,res)
    elif len(pth)==0:
        res.append(t)
    elif pth[0] in t:
        allValsList(t[pth[0]],pth[1:],res)
    return(res)

def allValsSet(t,pth,res=None):
    if(res==None):
        res=set()
    if (type(t)==list) | (type(t)==set):
        for s in t:
            allValsSet(s,pth,res)
    elif len(pth)==0:
        res.update([t])
    elif pth[0] in t:
        allValsSet(t[pth[0]],pth[1:],res)
    return(res)

###### failing.  not sure how to fix
#def allValsGen(t,pth):
#    if type(t)==list:
#        for s in t:
#            yield allValsGen(s,pth)
#    elif len(pth)==0:
#        yield t
#    elif pth[0] in t:
#        yield allValsGen(t[pth[0]],pth[1:])

def detailsToDF(pts):
    for t in pts:
        t['EVENTS']=pd.DataFrame.from_records(t['EVENTS'],index='EVENT_TIME')

#def dateFilter(pts,wdw=datetime.timedelta(-180),st=None):
#    if type(wdw)==int:
#        wdw=datetime.timedelta(days=wdw)
#    if type(st)==int:
#        st=datetime.timedelta(days=st)
#    if st is None:
#        st=[x['STUDY_START'] for x in pts]
#    mxd=datetime.timedelta(0)
#    i=0
#    while True:
#        delta=pts[i]['EVENTS'].index-st[i]
#        yield pts[i]['EVENTS'][(min(wdw,mxd)<=delta) & (delta<max(wdw,mxd))]
#        i+=1

def computeSparse(pts,pth,src=None,dmn=datetime.datetime(1,1,1),dmx=datetime.datetime(9999,1,1)):
    if type(src)=='str':
        src=set([src])
    ########## sparse matrices computed only from the events dictionary
    if pth[0]=='EVENTS':
        pth.pop(0)
    #### compute the dictionary of elements contained in path across all patients
    #### this creates a label for every type of source.  We really want to restrict
    #### this to the sources in src.  How to do?
    lbl=allValsSet(pts,['EVENTS']+pth)
    lbl=list(lbl)
    lbl=dict((value, idx) for idx,value in enumerate(lbl))
    
    drc=list()
    for i,s in enumerate(pts):
        if src is not None:
            try:
                t=[x for x in s['EVENTS'] if x['SOURCE'] in src and x['EVENT_TIME'] is not None and x['EVENT_TIME']>=dmn[i] and x['EVENT_TIME']<dmx[i]]
            except KeyError:
                continue
        else:
            t=s['EVENTS']
        U=allValsList(t,pth)
        if len(U)>0:
            ix=list()
            for u in U:
                ix.append(lbl[u])
            a,b=zip(*Counter(ix).items())
            drc+=zip(ones(i),a,b)
        if i%100==0:
            print(i)
    a,b,c=zip(*drc)
    M=scipy.sparse.coo_matrix( (np.array(c),(np.array(a),np.array(b))),shape=(len(pts),len(lbl)) )
    M=M.tocsr()
    sm=np.array(M.sum(0))[0]
    M=M[:,sm>0]
    lbl=list(lbl.items())
    lbl.sort(key=operator.itemgetter(1))
    lbl=[s[0] for s,t in zip(lbl,sm) if t>0]
    
    return lbl, M


# ---- source: vsEmrTools/encounterSummary.py ----
"""
Created on Mon Mar 12 11:19:12 2018

@author: jel2
"""

def encounterSummary(M,grp=dict()):
    ### add functionality: column identifying presence of group variables
    sTime=M.loc[M.groupby('ENCOUNTER_ID')['EVENT_TIME'].idxmin(),['PATIENT_ID','EVENT_TIME','ENCOUNTER_ID']].rename(columns={'EVENT_TIME':'eStartTime'})
    sTime=sTime.join(M.groupby('ENCOUNTER_ID')['EVENT_TIME'].max().rename('eEndTime'),on='ENCOUNTER_ID')
    
    return(sTime)

# ---- source: vsEmrTools/eventExpand.py ----

import pandas as pd

def eventExpand(eList,lbl,collapseToEncounter=False):
#    set(eList.SOURCE)
    cn=list(eList)
    cn.remove('DETAILS')
    a=eList.loc[eList.SOURCE==lbl,cn]
    dtl=eList.loc[eList.SOURCE==lbl,'DETAILS'].apply(lambda x: {} if x is None else x)
    a=a.join(pd.DataFrame.from_records(list(dtl),index=a.index))
    
    ##### collapse all event times to the date of first observation associated with encounter
    if collapseToEncounter:
        a=a.join(a.groupby(by='ENCOUNTER_ID')['EVENT_TIME'].agg(min).rename('enc_time'),on='ENCOUNTER_ID')
        a.loc[a.enc_time.isnull(),'EVENT_TIME']= a.loc[a.enc_time.isnull(),'EVENT_TIME']
    
    return(a)

# ---- source: vsEmrTools/firstEvent.py ----
"""
Created on Mon Mar 12 10:44:56 2018

@author: jel2
"""

def firstEvent(M,col,pat):
    
    import re
    
    sTime=M.loc[M[M[col].fillna('').str.contains(re.compile(pat,re.IGNORECASE))].groupby('PATIENT_ID')['EVENT_TIME'].idxmin(),['PATIENT_ID','ENCOUNTER_ID']].set_index('PATIENT_ID')
    return(sTime)

# ---- source: vsEmrTools/getStartTable.py ----
"""
Created on Thu May 17 13:44:25 2018

@author: jel2
"""

import datetime
import math
import pandas as pd



def getStartTable(fev,pts,eList):

    def computeFollowup(x):
        if (x['survival'] is not None) & (not math.isnan(x['survival'])):
            res=x['survival']
        elif(x['fu_collectionDate']<x['fu_lastRecord']):
            res=x['fu_collectionDate']
        else:
            res=.75*x['fu_collectionDate']+.25*x['fu_lastRecord']
#            res=x['fu_lastRecord']
        return(res)

    #    sTime=firstEvent(edd,'DISPOSITION','LWBS')
    sTime=encounterSummary(pd.merge(fev,eList,on='ENCOUNTER_ID',how='left'))
    sTime.set_index('PATIENT_ID',inplace=True)
    
    ########## Date (time) of the last time each patient was observed by the health system
    sTime=sTime.join(eList.groupby('PATIENT_ID')['EVENT_TIME'].max().rename('maxTime'))
    sTime=sTime.join(eList.groupby('PATIENT_ID')['EVENT_TIME'].min().rename('minTime'))
    sTime=sTime.join(pts[['DEATH_DATE','SEX','DOB','RACE']])
    sTime.RACE=sTime.RACE.apply(cleanRaceVariable)
    sTime['survival']=((sTime.DEATH_DATE-sTime.eStartTime)/datetime.timedelta(days=1))/365.25
    sTime['fu_lastRecord']=((sTime.maxTime-sTime.eStartTime)/datetime.timedelta(days=1))/365.25
    sTime['fu_collectionDate']=((mxt-sTime.eStartTime)/datetime.timedelta(days=1))/365.25
    
    sTime.loc[sTime.DOB>sTime.minTime,'minTime']=sTime.DOB[sTime.DOB>sTime.minTime]
    sTime['priorObsTime']=((sTime.eStartTime-sTime.minTime)/datetime.timedelta(days=1))/365.25
    
    sTime['DIED']= ~sTime.survival.isnull()
    sTime['followup']=sTime.apply(computeFollowup,axis=1)
    
    sTime['age']=(sTime.eStartTime-sTime.DOB)/datetime.timedelta(days=365.25)
    
    #### throw out patients who died before event start time
    sTime=sTime[sTime.followup>=0]
    
    return(sTime)

# ---- source: vsEmrTools/getSubtable.py ----
"""
Created on Wed Jun 14 10:05:36 2017

@author: jel2
"""

import pandas as pd

#### this version resets the index of X in the returned table
def getSubtable(X,dictColumn='DETAILS'):
    X=X.reset_index()
    t=list(X)
    t.remove(dictColumn)
    return(X[t].join(pd.DataFrame.from_records(X.DETAILS.tolist(),index=X.index)))
    

# ---- source: vsEmrTools/obsWindow.py ----


def obsWindow(eList):
    u=eList.groupby(by='PATIENT_ID')['EVENT_TIME'].aggregate([min,max])
    u.rename(columns={'min':'First Observation','max':'Last Observation'},inplace=True)
    return(u)
    

# ---- source: vsEmrTools/revenueByWeek.py ----
"""
Created on Thu May 17 13:44:25 2018

@author: jel2
"""

import math


def revenueByWeek(rev,sTime,low=-1,hi=1,w=None):
    R=timeFilter(rev,sTime,low,hi)
    R['week']=R.delta.apply(lambda x: math.floor(x*52))
    if not w:
        sTime['w']=1
        w='w'
    R=R.join(sTime[w],on='PATIENT_ID')
    t=R.groupby('week')[['AMOUNT',w]].apply(lambda x: sum(x['AMOUNT']*x[w])).rename('AMOUNT').reset_index()
    t.AMOUNT=-t.AMOUNT
    
    t.set_index('week',inplace=True)
    t['N']=0
    for _,a,b,q in sTime[['priorObsTime','followup','w']].itertuples():
        a=-min(52,math.floor(a*52))
        b=min(52,math.floor(b*52))
        t.loc[a:b,'N']+=q
    
#    t['N']=len(set(R.PATIENT_ID))
#    drp=[math.floor(i*52)+1 for i in sTime.followup[sTime.followup<hi]]
#    for i in drp:
#        t.loc[i:,'N']=t.loc[i:,'N'].apply(lambda x:x-1)
    t.reset_index(inplace=True)
    t.AMOUNT=t.AMOUNT/t.N
    t.drop('N',axis=1,inplace=True)
    
    return(t)

# ---- source: vsEmrTools/timeFilter.py ----
"""
Created on Tue Apr  3 08:36:56 2018

@author: jel2
"""

import datetime

def timeFilter(D,tm,low,hi,onVar='PATIENT_ID',tmVar='eStartTime',relTimeVar='delta',dTimeVar='EVENT_TIME'):
    D=D.join(tm[[tmVar]],on=onVar,how='inner')
    D[relTimeVar]=(D[dTimeVar]-D[tmVar])/datetime.timedelta(days=365.25)
    D.drop(tmVar,inplace=True,axis=1)
    D=D[(D[relTimeVar]>=low) & (D[relTimeVar]<hi)]
    return(D)

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:04:12 2018

@author: jel2
"""

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
    
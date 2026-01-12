# -*- coding: utf-8 -*-
"""
Created on Thu May 17 13:44:25 2018

@author: jel2
"""

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

# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 07:11:01 2016

@author: lucas
"""

def obsPerWeek(A,stLbl='firstVisitDate',edLbl='lastVisitDate',t0Lbl='eStartTime'):
    A=A[~A[t0Lbl].isnull()]
    a=(A[stLbl]-A[t0Lbl]).apply(lambda x: round(float(x.days)/7)).rename('start')
    b=(A[edLbl]-A[t0Lbl]).apply(lambda x: round(float(x.days)/7)).rename('end')
    d=pd.concat([a,b],axis=1)
    
    #### number of patients under observation in each week
    cts=pd.DataFrame({'week':range(min(d.start),max(d.end))})
    cts['N observed']=0
    
    for a,b,c in d.itertuples():
        cts.loc[(cts.week>=b) & (cts.week<=c),'N observed']+=1
        
    return(cts.set_index('week'))

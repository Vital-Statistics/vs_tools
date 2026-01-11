
import pandas as pd
from scipy.special.cython_special import betaln
#from scipy.stats import fisher_exact
from math import log

def llBBN(d,a,mu=.05,enObs=1):
    return(betaln(d+mu*enObs,a+(1-mu)*enObs)-betaln(mu*enObs,(1-mu)*enObs))

def beforeAfterComparison(tm,grp,cst=None,tThresh=30,wdw=1,maxPlot=10,patientID=None):
    
    tm.rename('delta',inplace=True)
    grp.rename('grp',inplace=True)
    if cst is None:
        cst=pd.Series(1,index=tm.index).rename('cst')
    else:
        cst.rename('cst',inplace=True)
    #grpr='Department'
    #grpr='LABEL'
#    a=pd.concat([tm.apply(lambda x: round(float(x.days)/7)),grp,cst,patientID],axis=1)
    a=pd.concat([tm.apply(lambda x: round(x/7)),grp,cst,patientID],axis=1)
    if(patientID is None):
        a['PATIENT_ID']=range(0,len(a))
    ####### wait times from study start
    
    smBefore=a[(a.delta>-tThresh) & (a.delta<-wdw)].groupby(by=['PATIENT_ID','grp'])['cst'].agg(max).reset_index()
    smBefore=smBefore.groupby(by=['grp'])['cst'].agg(sum).reset_index()
    smBefore['time window']='Before'
    
    smAfter=a[(a.delta<tThresh) & (a.delta>wdw)].groupby(by=['PATIENT_ID','grp'])['cst'].agg(max).reset_index()
    smAfter=smAfter.groupby(by=['grp'])['cst'].agg(sum).reset_index()
    smAfter['time window']='After'
    
    smTrans=a[(a.delta<=wdw) & (a.delta>=-wdw)].groupby(by=['PATIENT_ID','grp'])['cst'].agg(max).reset_index()
    smTrans=smTrans.groupby(by=['grp'])['cst'].agg(sum).reset_index()
    smTrans['time window']='Transition'
    
    q=pd.concat([smBefore,smAfter,smTrans])
    q.loc[q.cst<0,'cst']=0
    tbl=q.set_index(['grp','time window']).unstack('time window').fillna(0)['cst']
    tbl['OR']=((tbl['After']+1/len(tbl))/(1+sum(tbl['After'])))/((tbl['Before']+1/len(tbl))/(1+sum(tbl['Before'])))
    tbl['OR']=tbl['OR'].apply(lambda x: log(x,10))
    tbl['N']=tbl['Before']+tbl['After']
    
    sb=sum(tbl['Before'])
    sa=sum(tbl['After'])
    splt=tbl.Before.apply(lambda x: llBBN(x,sb-x))+tbl.After.apply(lambda x: llBBN(x,sa-x))
    jn=(tbl.Before+tbl.After).apply(lambda x: llBBN(x,sb+sa-x))
    tbl['BF']=splt-jn
    tbl=tbl.sort_values('BF',ascending=False)
#    tbl['Fisher exact']=tbl.apply(lambda x: fisher_exact([[x['Before'],x['After']],[sb,sa]])[1],axis=1)

    return(tbl)

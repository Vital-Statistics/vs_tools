# Auto-generated from vsOutcomeTools/ package modules.

# ---- source: vsOutcomeTools/addPlot_event.py ----
"""
Created on Thu Jul  6 16:17:05 2017

@author: jel2
"""
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

def cBands_poisGamma(u,alpha=.5,beta=.5,ivl=.95):
    a=u['Events']+alpha
    b=u['N']+beta
    delta=(1-ivl)/2
    u['low']=stats.gamma.ppf(delta,a,scale=1/b)
    u['hi']=stats.gamma.ppf(1-delta,a,scale=1/b)
    return(u)

def addPlot_event(fig,w,n,e,clr='k',lbl=None):
    u=pd.DataFrame({'week':w, 'N':n, 'Events':e})
    u=u.sort_values('week')
    u['N']=u['N'].cumsum()
    u['Events']=u['Events'].cumsum()
    u=u.apply(cBands_poisGamma,axis=1)
    plt.plot(u['week'],u['Events']/u['N'],color=clr,label=lbl)
    x=list(u.week)+list(u.week)[::-1]
    y=list(u.low)+list(u.hi)[::-1]
    fig.gca().fill(x,y,facecolor=clr,edgecolor=None,alpha=.3)

# ---- source: vsOutcomeTools/beforeAfterComparison.py ----

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

# ---- source: vsOutcomeTools/departmentUtilization.py ----
"""
Created on Mon Aug  1 07:11:01 2016

@author: lucas
"""

import re
import pandas as pd

def departmentUtilization(eList,pts,dSet):
    
    t=eList.loc[eList.SOURCE.isin({'ENCOUNTER','INPATIENT'}) & eList.DEPARTMENT_ID.isin(dSet),['PATIENT_ID','EVENT_TIME','ENCOUNTER_ID']]
    
    #### only one event per hospital encounter.  Set time to first time for hospital encounter
    t=t.groupby(by=['ENCOUNTER_ID','PATIENT_ID']).agg(min).reset_index().set_index('PATIENT_ID')
    t=pd.merge(t,pts[['STUDY_START']],left_index=True,right_index=True)
    
    #### compute event time as offsets from study start time
    t['week']=(t.EVENT_TIME-t.STUDY_START).apply(lambda x: round(float(x.days)/7))
    
    return(pd.DataFrame(t.week.groupby(t.week).size().rename('events')))

# ---- source: vsOutcomeTools/dictToExcel.py ----
#!/usr/bin/env python3
"""
Created on Mon Nov 27 09:02:53 2023

@author: rudy
"""

import pandas as pd

def dictToExcel(d,outFile):
    with pd.ExcelWriter(outFile, engine='xlsxwriter') as writer:
        for sheet_name, df in d.items():
            df.to_excel(writer, sheet_name=sheet_name)
        # writer.save()

# ---- source: vsOutcomeTools/multiTest.py ----
"""
Created on Mon Mar 19 09:25:16 2018

@author: jel2
"""

import numpy as np
import pandas as pd
import scipy
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
    

# ---- source: vsOutcomeTools/obsPerWeek.py ----
"""
Created on Mon Aug  1 07:11:01 2016

@author: lucas
"""

import pandas as pd

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

# ---- source: vsOutcomeTools/paidBarChart.py ----

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def paidBarChart(tbl,maxPlot=10):
    
    if(len(tbl)>maxPlot):
        kp=tbl.iloc[:maxPlot,:]
    kp=set(tbl.index)
    q=pd.DataFrame(tbl[['After','Transition','Before']].stack().rename('cst')).reset_index()
    q=q[q['grp'].apply(lambda x: x in kp)]
    
    fig, ax = plt.subplots(figsize=(8, round(maxPlot/2)+1), dpi=80)
    fig=sns.barplot(x='cst',y='grp',hue='time window',data=q,ax=ax)
    return(fig)

# ---- source: vsOutcomeTools/tableOne.py ----
"""
Created on Mon Feb 26 12:10:57 2018

@author: gsl8
"""

###############################################################################
def columns(df,d,i,out):

    import numpy as np

    if df[i].dtype is np.dtype('O'):

        c = {x: {d: {} for d in sorted(d.keys())} for x in sorted(set(df.loc[~df[i].isnull(),i]))}

        for c_key in sorted(c.keys()):
            a=[]
            for d_key in sorted(d.keys()):
                c[c_key][d_key]['num'] = sum(d[d_key]['df'][i]==c_key)
                c[c_key][d_key]['perc'] = '{:.1f}'.format(round(c[c_key][d_key]['num']/d[d_key]['num'],3)*100)
                a.append(str(c[c_key][d_key]['num'])+' ('+str(c[c_key][d_key]['perc'])+')')
            out[str(c_key.capitalize())+' - n.(%)'] = a

    elif df[i].dtype is np.dtype('bool'):

        c = {True: {d: {} for d in sorted(d.keys())}}
        a=[]
        for d_key in sorted(d.keys()):
            c[True][d_key]['num'] = sum(d[d_key]['df'][i]==True)
            c[True][d_key]['perc'] = '{:.1f}'.format(round(c[True][d_key]['num']/d[d_key]['num'],3)*100)
            a.append(str(c[True][d_key]['num'])+' ('+str(c[True][d_key]['perc'])+')')
        out[i.capitalize()+' - n.(%)'] = a

    elif df[i].dtype is np.dtype('float'):

        c = {i:{d:{} for d in sorted(d.keys())}}
        a=[]
        for d_key in sorted(d.keys()):
            c[i][d_key]['avg'] = round(np.mean(d[d_key]['df'][i]),1)
            c[i][d_key]['sd'] = round(np.std(d[d_key]['df'][i]),1)
            a.append(str(c[i][d_key]['avg'])+' ('+str(c[i][d_key]['sd'])+')')
        out[str(i.capitalize())+' - avg.(sd.)'] = a

###############################################################################
def demTable(df, splitColumn=None, output=None, filePath=None):

    '''
    Creates 'Table 1' or demographics table.

    df = Pandas dataframe
    splitColumn = Column used to group variables.
    output = If none, output is dataframe.  Other formats: latex, html, and csv.
    '''

    import numpy as np
    import pandas as pd
    import datetime

    if not splitColumn:
        df['ukjdsfg']=''
        splitColumn='ukjdsfg'

    pd.options.display.max_colwidth=200
    out = pd.DataFrame()

    d = {x: {'df':pd.DataFrame()} for x in sorted(set(df[splitColumn]))}

    a = []
    for key in sorted(d.keys()):
        d[key]['df'] = df.loc[df[splitColumn]==key,:] # populate each dataframe as subset of df
        d[key]['num'] = d[key]['df'].shape[0]           # N of each group
        a.append('N = '+str(d[key]['num']))

    out['Chacteristics'] = a # display N of each group


    for i in sorted(list(df.drop(splitColumn,axis=1))):
        columns(df,d,i,out=out)

    out = out.transpose()

    out_col = [x for x in sorted(d.keys())]

    for i in range(0,len(d.keys())):
        out.rename(columns={i:str(out_col[i].capitalize())},inplace=True)

    if output=='csv':
        out.to_csv(filePath)
    elif output=='latex':
        return print(out.to_latex(escape=False))
    elif output=='html':
        return print(out.to_html(escape=False))
    elif output==None:
        return out
###############################################################################
#### Example:
#test = eList[eList['SOURCE']=='MEDICATION']
#test = test.merge(pts, left_on='PATIENT_ID', right_index=True)
#test['DRUG'] = 'NEITHER'
#test.loc[test['LABEL'].str.contains(re.compile('metformin',re.IGNORECASE),re.IGNORECASE),'DRUG'] = 'METFORMIN'
#test.loc[test['LABEL'].str.contains(re.compile('sulf',re.IGNORECASE),re.IGNORECASE),'DRUG'] = 'SULFONYLUREA'
#Counter(test['DRUG'])
#test['DEATH'] = ~test['DEATH_DATE'].isna()
#Counter(test['DEATH'])
## This does not actually calculate age, just using this as an example of numerical data.
#test['AGE'] = (test['EVENT_TIME'] - test['DOB'])/datetime.timedelta(days=365.25)
#test = test[['DRUG','SEX','AGE','DEATH']]
#
#demTable(df=test,splitColumn='DRUG',output=None)
#demTable(df=test,splitColumn='DRUG',output='latex')
#demTable(df=test,splitColumn='DRUG',output='html')
#demTable(df=test,splitColumn='DRUG',output='csv',filePath='/Users/samlusk/Desktop/test.csv')

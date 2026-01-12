
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

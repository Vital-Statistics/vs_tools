# -*- coding: utf-8 -*-
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

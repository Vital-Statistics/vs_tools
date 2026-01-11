# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:44:56 2018

@author: jel2
"""

def firstEvent(M,col,pat):
    
    import re
    
    sTime=M.loc[M[M[col].fillna('').str.contains(re.compile(pat,re.IGNORECASE))].groupby('PATIENT_ID')['EVENT_TIME'].idxmin(),['PATIENT_ID','ENCOUNTER_ID']].set_index('PATIENT_ID')
    return(sTime)
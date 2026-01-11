# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 11:19:12 2018

@author: jel2
"""

def encounterSummary(M,grp=dict()):
    ### add functionality: column identifying presence of group variables
    sTime=M.loc[M.groupby('ENCOUNTER_ID')['EVENT_TIME'].idxmin(),['PATIENT_ID','EVENT_TIME','ENCOUNTER_ID']].rename(columns={'EVENT_TIME':'eStartTime'})
    sTime=sTime.join(M.groupby('ENCOUNTER_ID')['EVENT_TIME'].max().rename('eEndTime'),on='ENCOUNTER_ID')
    
    return(sTime)

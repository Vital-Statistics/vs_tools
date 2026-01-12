# -*- coding: utf-8 -*-
"""
Created on Thu May 17 13:44:25 2018

@author: jel2
"""

import datetime
import math
import pandas as pd

from vsEmrTools.encounterSummary import encounterSummary


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

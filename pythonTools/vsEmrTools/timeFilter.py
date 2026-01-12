# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 08:36:56 2018

@author: jel2
"""

import datetime

def timeFilter(D,tm,low,hi,onVar='PATIENT_ID',tmVar='eStartTime',relTimeVar='delta',dTimeVar='EVENT_TIME'):
    D=D.join(tm[[tmVar]],on=onVar,how='inner')
    D[relTimeVar]=(D[dTimeVar]-D[tmVar])/datetime.timedelta(days=365.25)
    D.drop(tmVar,inplace=True,axis=1)
    D=D[(D[relTimeVar]>=low) & (D[relTimeVar]<hi)]
    return(D)
    

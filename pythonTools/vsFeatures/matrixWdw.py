# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 09:49:53 2018

@author: jel2
"""

from vsFeatures.buildSparse import buildSparse

def matrixWdw(eList,params):
    import pandas as pd
    
    if not params:
        params=dict()
    if 'wdw' not in params:
        params['wdw']=365.25
    if 'SOURCE' not in params:
        params['SOURCE']={'DIAGNOSIS'}
    if 'colName' not in params:
        params['colName']='LABEL'
    if 'rowName' not in params:
        params['rowName']='PATIENT_ID'
    if 'tmCol' not in params:
        params['tmCol']='delta'
    if 'minPt' not in params:
        params['minPt']=10
    
    ll=(eList[params['tmCol']]<0) & (eList[params['tmCol']]> -params['wdw']) & (eList.SOURCE.isin(params['SOURCE']))
    r,c,M=buildSparse(eList[ll],params['rowName'],params['colName'])
    M=pd.DataFrame(M.todense(),columns=c,index=r) ###M.todense() to convert sparse matrix to full matrix###
    M=M.loc[:,M.sum(axis=0)>=params['minPt']]
    return(M)

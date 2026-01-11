#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 15:10:29 2023

@author: rudy
"""

def tOne_cleanBinary(tbl):
    import pandas as pd
    import math
    if 'Variable' not in list(tbl):
        tbl=tbl.reset_index()
    if 'Variable' not in list(tbl):
        print('Unable to find "Variable" column. Perhaps run tOne_toDF to convert default tableone output?')
        return
    if 'Value' not in list(tbl):
        print('Unable to find "Value" column.')
        return
    
    addRows=list()
    for rw in tbl.Variable.unique():
        t=tbl.loc[tbl.Variable==rw]
        if (len(t)==2) and (set(t.Value) in [{'1.0','0.0'},{'True','False'},{'0','1'},{'Y','N'}]):
            nMissing=[v for v in t.Missing if v!='']
            if 'P-Value' in list(t):
                pv=[v for v in t['P-Value'] if v!='']
                if len(pv)!=1:
                    print('P-Values are strange for variable: '+rw)
                    return
                else:
                    pv=pv[0]
            if len(nMissing)!=1:
                print('Missing counts are strange for variable: '+rw)
                return
            else:
                nMissing=nMissing[0]
            newRow=t.loc[t.Value.isin(['1.0','True','1','Y'])].copy()
            newRow.Missing=nMissing
            if 'P-Value' in list(t):
                newRow['P-Value']=pv
            newRow.Value=''
            addRows+=[newRow]
            tbl.drop(t.index.values,inplace=True)
    tbl=pd.concat([tbl]+addRows)
    tbl=tbl.set_index(['Variable','Value'])
    return(tbl)

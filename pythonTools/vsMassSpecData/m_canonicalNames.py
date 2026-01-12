#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 20:56:27 2025

@author: rudy
"""

import os
import re
import pandas as pd

from vsMassSpecData.apiLookups import apiLookups

def t_tg(s):
    res=re.sub(r'^TG\.(\d+)\.(\d+)_([\d\.]+)\.(\d+)', r'TG \1:\2_\3:\4', s)
    return(res)

def t_lpc(s):
    res=re.sub(r'^LPC[\s\.](\d+)[\s\.:](\d+)',r'LPC \1:\2',s)
    return(res)

def t_PC(s):
    res=re.sub(r'^PC[\s\.](\d+)[\s\.:](\d+)',r'PC \1:\2',s)
    return(res)

def t_PCO(s):
    res=re.sub(r'^PC[\s\.]O[\s\.](\d+)[\s\.:](\d+)',r'PC O-\1:\2',s)
    return(res)

def t_Hex2(s):
    res=re.sub(r'^Hex2Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex2Cer \1:\2/\3:\4',s)
    return(res)

def t_Hex3(s):
    res=re.sub(r'^Hex3Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex3Cer \1:\2/\3:\4',s)
    return(res)

def t_Hex(s):
    res=re.sub(r'^Hex[\s\.]Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex-Cer \1:\2/\3:\4',s)
    return(res)

def t_DG(s):
    res=re.sub(r'^DG[\s\.](\d+)[\s\.:](\d+)[\s\.:_](\d+)[\s\.:](\d+)',r'DG \1:\2_\3:\4',s)
    return(res)

def t_Cer(s):
    res=re.sub(r'^Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Cer \1:\2/\3:\4',s)
    return(res)

def t_CE(s):
    res=re.sub(r'^CE[\s\.](\d+)[\s\.:](\d+)',r'CE \1:\2',s)
    return(res)

def t_simple(s):
    res=re.sub(r'^(SM|FA)[\s\.](\d+)[\s\.:](\d+)',r'\1 \2:\3',s)
    return(res)

def t_PLA(s):
    res=re.sub(r'^PLA2[\s\.]Activity[\s\.:]+(\d+)[\s\.:]*',r'PLA2 Activity (\1)',s)
    return(res)

def m_canonicalNames(res,colName='oName',metaAppend=[]):
    import numpy as np
    syn=pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/synonyms.parquet')
    syn=syn.rename(columns={'Analyte':'oldName','New Name':'Analyte'})
    if isinstance(res, pd.Series):
        res=pd.DataFrame(res.rename('oName'))
    else:
        res=res.rename(columns={colName:'oName'})
    nov=pd.Series(index=list(syn))
    for col in nov.index.values:
        nov.loc[col]=res.oName.isin(syn[col]).sum()/len(res)
    if nov.max()>.8:
        res=res.join(syn.set_index(nov.sort_values(ascending=False).index.values[0]),on='oName')
    else:
        res['Analyte']=res.oName.apply(lambda x:t_tg(x))
        res.Analyte=np.where(res.Analyte.str.startswith('Ratio'),res.Analyte.str.replace('.',' '),res.Analyte)
        res.Analyte=np.where(res.Analyte.str.startswith('Sum'),res.Analyte.str.replace('.',' '),res.Analyte)
        res.Analyte=res.Analyte.apply(lambda x:t_PC(x))
        res.Analyte=res.Analyte.apply(lambda x:t_PCO(x))
        res.Analyte=res.Analyte.apply(lambda x:t_Cer(x))
        res.Analyte=res.Analyte.apply(lambda x:t_Hex2(x))
        res.Analyte=res.Analyte.apply(lambda x:t_Hex3(x))
        res.Analyte=res.Analyte.apply(lambda x:t_Hex(x))
        res.Analyte=res.Analyte.apply(lambda x:t_CE(x))
        res.Analyte=res.Analyte.apply(lambda x:t_simple(x))    
        res.Analyte=res.Analyte.apply(lambda x:t_lpc(x))
        res.Analyte=res.Analyte.apply(lambda x:t_PLA(x))
        res.Analyte=res.Analyte.apply(lambda x:t_DG(x))
        res.Analyte=res.Analyte.str.replace('.',' ')
        res.Analyte=res.Analyte.str.strip()
        res=res.join(syn.set_index('Analyte'),on='Analyte')
    
    metaData=list()
    if 'ros' in metaAppend: 
        metaData+=[pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/rosetta.parquet')]
    if 'api' in metaAppend:
        print('Searching pubchem and metabolomics workbench')
        metaData+=[apiLookups(s).drop(columns='Search Date')]
    if 'fnins-16-804216' in metaAppend:
        metaData+=[pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/papers/fnins-16-804216.parquet')]
    if 'hmdb' in metaAppend:
        metaData+=[pd.read_parquet(os.environ['DATA_PATH']+'HMDB/hmdbDirect.parquet')]
    if len(metaData)>0:
        ros=pd.concat(metaData,ignore_index=True).fillna('')
        rosList=list()
        for col in [v for v in list(ros) if v !='Analyte']:
            rosList+=[ros.groupby('Analyte')[col].apply(lambda x:list(set([v for v in x if v!=''])))]
        ros=pd.concat(rosList,axis=1)
        res=res.join(ros[[v for v in list(ros) if v not in {'Class'}]],on='Analyte').fillna('')
    
    return(res)

# #########
# print(res.Analyte.isin(syn.Analyte).sum())
# print(res.Analyte.nunique())

# for col in list(syn):
#     print(col+': '+str(syn[col].isin(lbl.index.values).sum()))

# st='Cer'
# for tbl in [res,syn]:
#     print('\n')
#     print(tbl.loc[tbl.Analyte.str.startswith(st)].head())

# res.loc[~res.Analyte.isin(list(syn.Analyte)) & res.Analyte.str.startswith(st)]

# # res.loc[~res.Analyte.isin(lbl.Analyte),'Analyte'].str[:3].value_counts().head()
# res.loc[~res.Analyte.isin(syn['Analyte']),'Analyte'].str[:3].value_counts().head()

# res.loc[~res.Analyte.isin(syn['Analyte']),'Analyte']

# st='PC '
# res.loc[~res.Analyte.isin(syn.Analyte) & res.Analyte.str.startswith(st)]

# lbl.loc[lbl['Short Name'].str.startswith('PC '),'Short Name'].str[:6].value_counts()

# for a in list(ros):
#     for b in list(syn):
#         print('ros['+a+'], syn['+b+']: '+str(len(set(ros[a]) & set(syn[b]))))
    

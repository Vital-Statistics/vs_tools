#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 08:49:06 2025

@author: rudy
"""

def loadSynonyms(meta):
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    
    syn=pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/synonyms.parquet')
    
    return(meta.join(lbl,on='Compound'))

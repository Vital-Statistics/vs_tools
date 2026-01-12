#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 15:48:11 2023

@author: rudy
"""

import os
import pandas as pd

def addHMDB(meta):
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    return(meta.join(lbl,on='Compound'))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 18:25:56 2025

@author: rudy
"""

def listKeggPathways(org='hsa'):
    import requests
    import pandas as pd
    BASE = "http://rest.kegg.jp"
    rr = requests.get(f"{BASE}/list/pathway/{org}", timeout=20)
    rr.raise_for_status()
    res = {}
    for ln in rr.text.strip().splitlines():
        pid, name = ln.split("\t", 1)
        res[pid] = name
    res=pd.DataFrame(res,index=[0]).T.reset_index()
    res.columns=['Pathway ID','Pathway Name']
    res['Pathway ID']='path:'+res['Pathway ID']
    return(res.set_index('Pathway ID'))

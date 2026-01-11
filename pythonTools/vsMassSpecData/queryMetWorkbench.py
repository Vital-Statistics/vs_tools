#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 15:15:46 2025

@author: rudy
"""


# List of metabolite names to standardize
def queryMetWorkbench(mwl):
    import requests
    res=list()
    for i,name in enumerate(mwl):
        loopProgress(i,len(mwl),'matching')
        # Call the RefMet API's match endpoint for each name
        url = f"https://www.metabolomicsworkbench.org/rest/refmet/match/{name}"
        response = requests.get(url)
        res += [response.json()]  # Parse the JSON response
    
    res=pd.DataFrame(res,index=mwl)
    
    nmList=res.loc[res.refmet_name!='-','refmet_name'].unique()
    dList=list()
    for n,std_name in enumerate(nmList):
        rw=v[1]
        loopProgress(n,len(res),'querying details')
        # print(std_name)
        dList += [requests.get(f"https://www.metabolomicsworkbench.org/rest/refmet/name/{std_name}/all").json()]
    dets=pd.DataFrame(dList)
    
    dups=set(res) & set(dets)
    res=res.join(dets.set_index('name').drop(columns=list(dups)),on='refmet_name')
    return(res)

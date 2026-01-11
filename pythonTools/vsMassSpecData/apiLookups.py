#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 15:35:53 2025

@author: rudy
"""

import requests
import time


def lookupKeggPathways(keggID,reload=False):
    pwFilePath=os.environ['DATA_PATH']+'KEGG/keggPathways.parquet'
    pw=pd.read_parquet(pwFilePath)
    if keggID in set(pw.kegg) and not reload:
        return(pw.set_index('kegg').at[keggID,'pathways'])
    else:
        url = f"https://rest.kegg.jp/link/pathway/{keggID}"
        r = requests.get(url)
        if not r.text.strip():
            return []
        time.sleep(.5)
        res=[line.split("\t")[1] for line in r.text.strip().split("\n")]
        pw=pd.concat([pw.loc[pw.kegg!=keggID],pd.DataFrame({'kegg':keggID,'pathways':[res]},index=[0])],ignore_index=True)
        pw.to_parquet(pwFilePath)
        return(res)

def uniprot_to_kegg(uniprot_id):
    ak=pd.read_parquet(os.environ['DATA_PATH']+'KEGG/allKegg.parquet')
    if uniprot_id in list(ak.uniprot):
        return(ak.set_index('uniprot').loc[[uniprot_id],'kegg'].tolist())
    else:
        url = "https://rest.kegg.jp/conv/hsa/uniprot:" + uniprot_id
        r = requests.get(url)
        r.raise_for_status()
    
        kegg = [v.strip().split("\t")[-1] for v in r.text.split('\n') if v !='']
        ak=pd.concat([ak,pd.DataFrame({'kegg':kegg,'uniprot':uniprot_id},index=range(len(kegg)))],ignore_index=True)
        ak.to_parquet(os.environ['DATA_PATH']+'KEGG/allKegg.parquet')
        time.sleep(.75)
        return(kegg)

def lookup_hmdb(hmdb_id):
    url = f"https://www.metabolomicsworkbench.org/rest/compound/hmdb_id/{hmdb_id}/"
    r = requests.get(url)
    return r.json()

def map_names_via_metaboanalyst(names):
    url = "https://rest.xialab.ca/api/mapcompounds"
    payload = {
        "queryList": ";".join(names),
        "inputType": "name"
    }
    return requests.post(url, json=payload).json()

import pubchempy as pcp

def lookup_pubchem_by_name(name, max_hits=3):
    comps = pcp.get_compounds(name, 'name')
    results = []
    for c in comps[:max_hits]:
        results.append({
            "query_name": name,
            "pubchem_cid": c.cid,
            "iupac_name": c.iupac_name,
            "formula": c.molecular_formula,
            "inchikey": c.inchikey
        })
    return results

def mw_lookup_pubchem(cid):
    url = f"https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/{cid}/all/json"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()

def apiLookups(s,rerunSearch=180):
    from datetime import date,timedelta
    import os, time
    import pandas as pd
    
    psPath=os.environ['DATA_PATH']+'metabolomics/apiLookups.parquet'
    priorSearch=pd.read_parquet(psPath)
    priorSearch["Search Date"] = pd.to_datetime(priorSearch["Search Date"], errors="coerce")

    if rerunSearch:
        cutoff = pd.Timestamp.today().normalize() - pd.Timedelta(days=rerunSearch)
        tbl=priorSearch.loc[priorSearch.query_name.isin(list(s)) & (priorSearch['Search Date']>=cutoff)].copy()
        print(tbl)
        s=s.loc[~s.isin(tbl.query_name)]
        pcList=list()
        mwList=list()
        for v in list(s):
            print(v)
            out=lookup_pubchem_by_name(v)
            if len(out)>0:
                pcList+=out
                for v in out:
                    dets=mw_lookup_pubchem(v['pubchem_cid'])
                    if len(dets)>0:
                        mwList+=dets
            time.sleep(2)
        if len(pcList)>0:
            pc=pd.DataFrame(pcList)
            pc.pubchem_cid=pc.pubchem_cid.astype(str)
            if len(mwList)>0:
                mw=pd.DataFrame(mwList)
                mw.pubchem_cid=mw.pubchem_cid.astype(str)
                pc=pd.merge(pc,mw[[v for v in list(mw) if v not in ['formula','inchi_key']]],on='pubchem_cid',how='left')
            pc['Search Date']=pd.Timestamp.today().normalize()
            pc=pc.fillna('')
            
            tbl=pd.concat([tbl,pc],ignore_index=True)
            priorSearch=pd.concat([priorSearch,pc],ignore_index=True)
            priorSearch=priorSearch.fillna('')
            
            priorSearch.to_parquet(psPath)
        tbl=tbl.fillna('')
    else:
        tbl=priorSearch.loc[priorSearch.query_name.isin(list(s))]
    tbl=tbl.rename(columns={'query_name':'Analyte'})
    return(tbl)

# for res in lookup_pubchem_by_name("alanine"):
#     print(res)

# for res in lookup_pubchem_by_name("choline"):
#     print(res)



# # Example usage
# print(map_names_via_metaboanalyst(["alanine","choline"]))
# print(lookup_hmdb("HMDB0000161"))

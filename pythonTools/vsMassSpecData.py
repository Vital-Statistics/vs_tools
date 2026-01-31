"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""

import os
import re
import time
from datetime import date, timedelta

import numpy as np
import pandas as pd
import pubchempy as pcp
import requests

from vsVisualizations import loopProgress

def addHMDB(meta):
    """
    Add HMDB IDs to a metadata table using a Rosetta Stone mapping.

    Parameters
    ----------
    meta : pandas.DataFrame
        Metadata table with a 'Compound' column.

    Returns
    -------
    pandas.DataFrame
        Metadata joined with HMDB ID lists.
    """
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    return(meta.join(lbl,on='Compound'))

def lookupKeggPathways(keggID,reload=False):
    """
    Look up KEGG pathways for a KEGG compound ID with local caching.

    Parameters
    ----------
    keggID : str
        KEGG ID (e.g., 'C00031' or 'hsa:1234').
    reload : bool, default False
        If True, force a fresh API lookup and update the cache.

    Returns
    -------
    list[str]
        Pathway IDs linked to the KEGG entry.
    """
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
    """
    Convert a UniProt ID to KEGG IDs with local caching.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession.

    Returns
    -------
    list[str]
        KEGG IDs for the UniProt accession.
    """
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
    """
    Look up HMDB details from Metabolomics Workbench.

    Parameters
    ----------
    hmdb_id : str
        HMDB identifier.

    Returns
    -------
    dict
        JSON response from the API.
    """
    url = f"https://www.metabolomicsworkbench.org/rest/compound/hmdb_id/{hmdb_id}/"
    r = requests.get(url)
    return r.json()

def map_names_via_metaboanalyst(names):
    """
    Map compound names via MetaboAnalyst REST API.

    Parameters
    ----------
    names : list[str]
        Compound names to map.

    Returns
    -------
    dict
        JSON response from the API.
    """
    url = "https://rest.xialab.ca/api/mapcompounds"
    payload = {
        "queryList": ";".join(names),
        "inputType": "name"
    }
    return requests.post(url, json=payload).json()


def lookup_pubchem_by_name(name, max_hits=3):
    """
    Query PubChem by compound name and return basic annotations.

    Parameters
    ----------
    name : str
        Compound name.
    max_hits : int, default 3
        Maximum number of hits to return.

    Returns
    -------
    list[dict]
        Annotated PubChem hits (CID, IUPAC name, formula, InChIKey).
    """
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
    """
    Look up Metabolomics Workbench details by PubChem CID.

    Parameters
    ----------
    cid : str | int
        PubChem CID.

    Returns
    -------
    dict
        JSON response from the API.
    """
    url = f"https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/{cid}/all/json"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()

def apiLookups(s,rerunSearch=180):
    """
    Perform cached PubChem and Metabolomics Workbench lookups for analytes.

    Parameters
    ----------
    s : pandas.Series
        Series of analyte names.
    rerunSearch : int, default 180
        Re-run searches older than this many days.

    Returns
    -------
    pandas.DataFrame
        Lookup results for the analytes.
    """
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

def loadMetabolomics(filePath,lastSampCol='Injection Number',sid='Customer Sample Identification',maxMissing=1
                     ,header=1,lodCol=None,nHeader=None,log2Transform=True,normalize=False,platform=None,stdFilter=False
                     ,imputeMissing=True):
    """
    Load metabolomics data from a vendor export and perform basic cleaning.

    Parameters
    ----------
    filePath : str
        Path to the input Excel file.
    lastSampCol : str, default 'Injection Number'
        Column name for the last sample metadata column.
    sid : str, default 'Customer Sample Identification'
        Sample ID column name.
    maxMissing : int, default 1
        Maximum missing values per analyte before dropping.
    header : int, default 1
        Excel header row.
    lodCol : str | None
        Optional LOD column name.
    nHeader : int | None
        Number of header rows to skip.
    log2Transform : bool | None
        Whether to log2-transform data (must be explicitly provided).
    normalize : bool, default False
        Whether to normalize analyte levels.
    platform : str | None
        Optional prefix for analyte names.
    stdFilter : float | False
        Optional standard deviation filter threshold.
    imputeMissing : bool, default True
        Whether to impute missing values.

    Returns
    -------
    tuple
        (ss, aa, t) where ss is sample metadata, aa is analyte annotation,
        and t is the analyte matrix.
    """
    
    # read in data
    t=pd.read_excel(filePath,header=header,dtype=str,na_filter=False)
    if lastSampCol.upper() in [v.upper() for v in list(t)]:
        nn=[i for i,v in enumerate(list(t)) if v==lastSampCol][0]+1
        print('Found "'+lastSampCol+'" in column '+str(nn-1)+'.  First analyte found is '+list(t)[nn])
        ss=t.iloc[:,:nn]
        mm=ss.loc[ss.iloc[:,0]!=''].index.values[0]
        aa=t.iloc[:mm,(nn-1):].set_index(lastSampCol)
        if nHeader is None:
            print('First non-null value in first column, "'+list(t)[0]+'", is in row '+str(mm))
        else:
            print('Skipping '+str(nHeader)+' rows of header.')
            mm=nHeader
        
        ss=ss.iloc[mm:,:]
        t=t.iloc[mm:,:]

        drpRow=(ss[sid]=='') | ss[sid].str.contains('pool|NIST SRM|Golden West|SPQC',flags=re.IGNORECASE)
        t=t.loc[~drpRow]
        ss=ss.loc[~drpRow]
        t.index=t[sid]
        ss.set_index(sid,inplace=True)
        t=t.iloc[:,nn:]
        drpCol=[v for v in list(t) if re.search('Status',v,re.IGNORECASE)]
        #### Consider modifying expression levels based on < LLOQ or > ULOQ
        kpCol=[v for v in list(t) if v not in drpCol]
        if len(drpCol)!=len(kpCol):
            print('Mismatch in the number of status and expression columns.')
        t=t[kpCol]
        aa=aa[kpCol].transpose()
        aa.columns=[str(v) for v in aa]
        if sum([a!=b for a,b in zip(aa.index.values,list(t))])>0:
            print('Mismatch in analyte metadata and expression labels')
            
        if lodCol is None:
            lodCol=[v for v in list(aa) if re.search('LOD',v)]
            if len(lodCol)==0:
                lodCol=[v for v in list(aa) if re.search('Lowest Calibration Standard',v,re.IGNORECASE)]
            print('Found '+str(len(lodCol))+' column(s) with LOD label.  Using "'+lodCol[0]+'".')
            lodCol=lodCol[0]
        aa[lodCol]=pd.to_numeric(aa[lodCol],errors='coerce')
        aa.loc[(aa[lodCol]<=0) | aa[lodCol].isna(),lodCol]=aa.loc[aa[lodCol]>0,lodCol].min()
        
        if (maxMissing>0) & (maxMissing<=1): maxMissing*=len(t)
        print('Dropping columns with >='+str(round(100*maxMissing/len(t),2))+'% missing data:')
        nDrop=0
        for cc in list(t):
            t[cc]=pd.to_numeric(t[cc],errors='coerce').astype(float)
            if imputeMissing:
                t.loc[t[cc]<aa.loc[cc,lodCol],cc]=aa.loc[cc,lodCol]/2
            if t[cc].isna().sum()>=maxMissing:
                t.drop(columns=cc,inplace=True)
                aa.drop(cc,inplace=True)
                print('\t'+cc)
                nDrop+=1
            else:
                if imputeMissing: t.loc[t[cc].isna(),cc]=float(aa.loc[cc,lodCol])/2
                
        print('\t'+str(t.shape[1])+' out of '+str(t.shape[1]+nDrop)+' metabolites kept.')
        if 'Compound Class' in list(aa) and 'Class' not in list(aa):
            aa.rename(columns={'Compound Class':'Class'},inplace=True)
            
        if normalize:
            print('Normalizing to total metabolites')
            t=pd.DataFrame(np.array(t)/np.array(t.sum(axis=1)).reshape(-1,1),columns=t.columns,index=t.index)
            
        if log2Transform:
            t = np.log2(t.where(t > 0))
        
        if platform:
            aa.reset_index(inplace=True)
            aa['Compound']=aa['index']
            aa['index']=platform+aa['index']
            aa.set_index('index',inplace=True)
            t.columns=[platform+v for v in list(t)]
            
        if stdFilter:
            kp=t.std()>stdFilter
            kp=kp.loc[kp].index.values
            print('Dropping '+str(t.shape[1]-len(kp))+' analytes due to standard deviation filter.')
            t=t[kp]
    else:
        print('Format unknown.  Returning unmodified.')
        ss=None
        aa=None
    return((ss,aa,t))
        
def loadSynonyms(meta):
    """
    Load synonym mappings and join with metadata.

    Parameters
    ----------
    meta : pandas.DataFrame
        Metadata table with a 'Compound' column.

    Returns
    -------
    pandas.DataFrame
        Metadata joined with synonym mappings.
    """
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    
    syn=pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/synonyms.parquet')
    
    return(meta.join(lbl,on='Compound'))

def t_tg(s):
    """Normalize triglyceride naming format."""
    res=re.sub(r'^TG\.(\d+)\.(\d+)_([\d\.]+)\.(\d+)', r'TG \1:\2_\3:\4', s)
    return(res)

def t_lpc(s):
    """Normalize LPC naming format."""
    res=re.sub(r'^LPC[\s\.](\d+)[\s\.:](\d+)',r'LPC \1:\2',s)
    return(res)

def t_PC(s):
    """Normalize PC naming format."""
    res=re.sub(r'^PC[\s\.](\d+)[\s\.:](\d+)',r'PC \1:\2',s)
    return(res)

def t_PCO(s):
    """Normalize PC-O naming format."""
    res=re.sub(r'^PC[\s\.]O[\s\.](\d+)[\s\.:](\d+)',r'PC O-\1:\2',s)
    return(res)

def t_Hex2(s):
    """Normalize Hex2Cer naming format."""
    res=re.sub(r'^Hex2Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex2Cer \1:\2/\3:\4',s)
    return(res)

def t_Hex3(s):
    """Normalize Hex3Cer naming format."""
    res=re.sub(r'^Hex3Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex3Cer \1:\2/\3:\4',s)
    return(res)

def t_Hex(s):
    """Normalize HexCer naming format."""
    res=re.sub(r'^Hex[\s\.]Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Hex-Cer \1:\2/\3:\4',s)
    return(res)

def t_DG(s):
    """Normalize DG naming format."""
    res=re.sub(r'^DG[\s\.](\d+)[\s\.:](\d+)[\s\.:_](\d+)[\s\.:](\d+)',r'DG \1:\2_\3:\4',s)
    return(res)

def t_Cer(s):
    """Normalize Cer naming format."""
    res=re.sub(r'^Cer[\s\.](d\d+)[\s\.:](\d+)[\s\.:](\d+)[\s\.:](\d+)',r'Cer \1:\2/\3:\4',s)
    return(res)

def t_CE(s):
    """Normalize CE naming format."""
    res=re.sub(r'^CE[\s\.](\d+)[\s\.:](\d+)',r'CE \1:\2',s)
    return(res)

def t_simple(s):
    """Normalize simple lipid naming formats."""
    res=re.sub(r'^(SM|FA)[\s\.](\d+)[\s\.:](\d+)',r'\1 \2:\3',s)
    return(res)

def t_PLA(s):
    """Normalize PLA2 naming format."""
    res=re.sub(r'^PLA2[\s\.]Activity[\s\.:]+(\d+)[\s\.:]*',r'PLA2 Activity (\1)',s)
    return(res)

def m_canonicalNames(res,colName='oName',metaAppend=[]):
    """
    Standardize analyte names and optionally append metadata.

    Parameters
    ----------
    res : pandas.DataFrame | pandas.Series
        Input analyte table or series.
    colName : str, default 'oName'
        Column name containing original analyte names.
    metaAppend : list[str]
        Optional metadata sources to append (e.g., 'ros', 'api', 'hmdb').

    Returns
    -------
    pandas.DataFrame
        Table with standardized analyte names and optional metadata.
    """
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

def queryMetWorkbench(mwl):
    """
    Query Metabolomics Workbench RefMet matching and details for names.

    Parameters
    ----------
    mwl : list[str]
        List of metabolite names.

    Returns
    -------
    pandas.DataFrame
        RefMet match results joined with details.
    """
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

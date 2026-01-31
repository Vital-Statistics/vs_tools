"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
# Auto-generated from vsMassSpecData/ package modules.

# ---- source: vsMassSpecData/addHMDB.py ----
#!/usr/bin/env python3


import os
import pandas as pd

def addHMDB(meta):
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    return(meta.join(lbl,on='Compound'))

# ---- source: vsMassSpecData/apiLookups.py ----
#!/usr/bin/env python3


import time
import os


def lookupKeggPathways(keggID,reload=False):
    import requests
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
    import requests
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
    import requests
    url = f"https://www.metabolomicsworkbench.org/rest/compound/hmdb_id/{hmdb_id}/"
    r = requests.get(url)
    return r.json()

def map_names_via_metaboanalyst(names):
    import requests
    url = "https://rest.xialab.ca/api/mapcompounds"
    payload = {
        "queryList": ";".join(names),
        "inputType": "name"
    }
    return requests.post(url, json=payload).json()


def lookup_pubchem_by_name(name, max_hits=3):
    import pubchempy as pcp
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
    import requests
    url = f"https://www.metabolomicsworkbench.org/rest/compound/pubchem_cid/{cid}/all/json"
    r = requests.get(url)
    r.raise_for_status()
    return r.json()

def apiLookups(s,rerunSearch=180):
    from datetime import date,timedelta
    import os, time
    
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

# ---- source: vsMassSpecData/loadMetabolomics.py ----


import re

def loadMetabolomics(filePath,lastSampCol='Injection Number',sid='Customer Sample Identification',maxMissing=1
                     ,header=1,lodCol=None,nHeader=None,log2Transform=None,normalize=False,platform=None,stdFilter=False
                     ,imputeMissing=True):
    import numpy as np
    import math
    
    if log2Transform is None:
        print('Code as been modified to require explicit log2Transform.')
        return()
    
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
            # t=t.apply(np.vectorize(math.log2))
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
        
    # bileAcids = pd.read_excel(filePath,
    #                           sheet_name='Status Columns Removed',
    #                           header=[1])
        
    # # extract limits of detection
    # LOD = pd.to_numeric(bileAcids.iloc[1,:],errors='coerce')
    
    # # drop extra rows
    # # bileAcids = bileAcids.drop([0,1,2,3])
    # drpRow=bileAcids['Customer Sample Identification'].str.contains('study pool',re.IGNORECASE).fillna(True)
    # bileAcids=bileAcids.loc[~drpRow]
    
    # # drop extra columns
    # drops = ['Plate Bar Code','Sample Bar Code','Sample Type','Sample Identification',
    #          'Species','Material','Well Position','Sample Volume','Run Number','Injection Number']
    # bileAcids=bileAcids.drop(drops,axis=1).rename(columns={'Customer Sample Identification':'Sample ID'})
    
    # # drop columns with greater than 40% <LOD
    # keeps=['Sample ID']
    # for k in list(bileAcids)[1:]:
    #     if sum((bileAcids[k]=='<LOD') | (bileAcids[k]==0) | (bileAcids[k].isnull()))/len(bileAcids[k])<0.40:
    #         keeps.append(k)
    
    # bileAcids = bileAcids[keeps]
    
    # # replace <LOD and NA with limits of detection
    # for k in list(bileAcids)[1:]:
    #     bileAcids.loc[(bileAcids[k]=='<LOD') | (bileAcids[k].isnull()) | (bileAcids[k]==0),k] = LOD[k]
    
    # # convert to float
    # bileAcids[list(bileAcids)[1:]] = bileAcids[list(bileAcids)[1:]].astype('float')
    # return(bileAcids)


# ---- source: vsMassSpecData/loadSynonyms.py ----
#!/usr/bin/env python3


import os

def loadSynonyms(meta):
    lbl=pd.concat([pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='LC Part',header=1),
                   pd.read_excel(os.environ['DATA_PATH']+'metabolomics/q500 metaboliteNames/Rosetta Stone Quant 500_BioIDs_20190121.xlsx',sheet_name='FIA Part',header=1)])
    lbl=lbl.loc[~lbl.HMDB.isna() & ~lbl['Short Name'].isna()][['Short Name','HMDB']].groupby('Short Name')['HMDB'].apply(lambda x:list(x))
    
    syn=pd.read_parquet(os.environ['DATA_PATH']+'metabolomics/synonyms.parquet')
    
    return(meta.join(lbl,on='Compound'))

# ---- source: vsMassSpecData/m_canonicalNames.py ----
#!/usr/bin/env python3


import os
import re


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
    

# ---- source: vsMassSpecData/queryMetWorkbench.py ----
#!/usr/bin/env python3





# List of metabolite names to standardize
def queryMetWorkbench(mwl):
    from vsVisualizations import loopProgress
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

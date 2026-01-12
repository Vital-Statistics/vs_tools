# Auto-generated from vsPathways/ package modules.

# ---- source: vsPathways/apiUniprot.py ----
#!/usr/bin/env python3
"""
Created on Sun Aug 31 07:31:24 2025

@author: rudy
"""

import pandas as pd

def parseUniProt(v):
    op={'Assession':v['primaryAccession']}
    op['UniProt ID']=v.get('uniProtkbId',{})
    op['Gene']=v.get("genes", [{}])[0].get("geneName", {}).get("value")

    protDesc=v.get("proteinDescription",{})
    if 'recommendedName' in protDesc:
        op['Protein name']=protDesc["recommendedName"]["fullName"]["value"]
    else:
        op['Protein name']=protDesc["submissionNames"][0]["fullName"]["value"]
        
    op['Organism']=v["organism"]["scientificName"]
    xrefs = v.get("uniProtKBCrossReferences", [])
    kid=list()
    for x in xrefs:
        if x["database"] == "KEGG":
            kid+=[x["id"]]
    op['KEGG ID']=frozenset(kid)
    return(op)

def apiUniprot(analyte,organism_id='9606'):
    import requests
    
    # Example: query UniProt for "IL10" (human)
    
    query = f'"{analyte}" AND organism_id:{organism_id}'  # 9606 = Homo sapiens
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    params = {
        "query": query,
        "fields": "id,accession,gene_names,protein_name,organism_name,go_p,xref_kegg",  # select fields
        "format": "json"
    }
    
    response = requests.get(url, params=params)
    upOut = response.json()["results"]
    
    tbl=pd.DataFrame([parseUniProt(v) for v in upOut],index=range(len(upOut)))

    return((tbl,upOut))

# ---- source: vsPathways/computeGAGE.py ----
"""
Created on Mon Dec 13 09:11:43 2021

@author: joest
"""

import math
import pandas as pd

def computeGAGE(pv):
    ### pv should be a list of p-values for analytes in the pathway
    import scipy.stats as stats
    pv=[max(1e-30,v) for v in pv if not pd.isna(v)]
    gam=stats.gamma(len(pv))
    return(max(1e-30,1-gam.cdf(sum([-math.log(v) for v in pv]))))


# ############## Testing for accuracy
# nSamp=10000
# sampSize=10
# gam=stats.gamma(sampSize)
# res=pd.Series(0,index=list(range(nSamp)))
# for i in res.index.values:
#     rr=np.random.uniform(size=sampSize)
#     # rr=-np.log(np.random.uniform(size=sampSize))
#     # rr=sum(rr)
#     # rr=np.mean(gam.rvs(100))
#     # res.loc[i]= 1-gam.cdf(rr)
#     res.loc[i]=computeGAGE(rr)
# res.hist(bins=np.arange(0,1.01,.01))


# np.mean(gam.rvs(10000))

# ---- source: vsPathways/keggPathways.py ----
#!/usr/bin/env python3
"""
Created on Sun Aug 31 14:51:35 2025

@author: rudy
"""

import pandas as pd



def keggPathways(kegg_id: str) -> pd.DataFrame:
    import requests
    import pandas as pd
    BASE = "http://rest.kegg.jp"
    # 1) get pathway IDs for this entry
    r = requests.get(f"{BASE}/link/pathway/{kegg_id}", timeout=20)
    r.raise_for_status()
    lines = [ln for ln in r.text.strip().splitlines() if ln.strip()]
    path_ids = sorted({ ln.split("\t")[1] for ln in lines })  # e.g., 'path:hsa04060'

    if not path_ids:
        return pd.DataFrame(columns=["pathway_id", "pathway_name"])

    # # 2) build a name map (organism-specific if gene; generic for compounds)
    # org = None
    # if ":" in kegg_id:
    #     prefix = kegg_id.split(":")[0]
    #     if prefix not in {"cpd", "drug", "glycan"}:
    #         org = prefix  # e.g., 'hsa', 'mmu'

    # name_map = {}

    # # organism-specific names (e.g., path:hsa04060)
    # if org:
    #     rr = requests.get(f"{BASE}/list/pathway/{org}", timeout=20)
    #     rr.raise_for_status()
    #     for ln in rr.text.strip().splitlines():
    #         pid, name = ln.split("\t", 1)
    #         name_map[pid] = name

    # # generic map names (e.g., path:map00010)
    # else:
    #     rr = requests.get(f"{BASE}/list/pathway", timeout=20)
    #     rr.raise_for_status()
    #     for ln in rr.text.strip().splitlines():
    #         pid, name = ln.split("\t", 1)
    #         name_map[pid] = name

    rows = [{"pathway_id": pid, "pathway_name": name} for pid in path_ids]
    return pd.DataFrame(rows)

# ---- source: vsPathways/listKeggPathways.py ----
#!/usr/bin/env python3
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

# ---- source: vsPathways/pwSMPDB.py ----
"""
Created on Mon Dec 13 10:05:56 2021

@author: joest
"""

import glob
import os
import numpy as np
import pandas as pd

from vsVisualizations import loopProgress

def pwSMPDB(M,mLbl='HMDB', pthLbl='HMDB ID',sigCol='p-value',st='t-statistic',minMeasured=3,dataType='metabolomics'
            ,pathwayType=None,dropPathwayList=None,computation='GAGE'):
    import math
    from collections import defaultdict
    import gseapy as gp
    from statsmodels.stats.multitest import multipletests
    
    M=M.sort_values(sigCol)
    rnk = pd.Series(data=np.linspace(len(M), 1, len(M)),index=M[mLbl])
    pth=os.environ['DATA_PATH']+dataType+'/pathways/smpdb/'
    V=pd.read_csv(pth+'smpdb_pathways.csv').set_index('SMPDB ID')
    if pathwayType is not None:
        V=V.loc[V.Subject.isin(pathwayType)]
    if dropPathwayList is not None:
        print('Removing pathways: '+', '.join(dropPathwayList))
        V=V.loc[~V.Subject.isin(dropPathwayList)]
    V['pathway ID']=V.index.values
    V['Measured Components']=''
    V['p-value']=np.nan
    V['Decrease']=0
    V['No change']=0
    V['% up']=0.0
    V['Increase']=0
    V['N Measured']=0
    V['Pathway Size']=0
    V['Average Change']=0.0
    V['Median Change']=0.0
    V['Increase + P<.05']=0
    V['Decrease + P<.05']=0
    V['Lead Genes']=''
    drpList=list()
    
    fList=pd.DataFrame({'filePath':glob(pth+'pathwayLists/*.csv')})
    fList['SMPDB ID']=fList.filePath.str.split(os.sep).str[-1].str.split('_').str[0]
    fList=fList.set_index('SMPDB ID')
    V=V.join(fList)
    V=V.loc[~V.filePath.isna()]

    pthKeys=defaultdict(lambda:None)
    sSize=math.ceil(len(V)/30)
    for i,(pathway,fl) in enumerate(V.filePath.items()):
        if i%sSize==0:
            loopProgress(int(i/sSize),30,'Pathway Analysis')
        t=pd.read_csv(fl)
        V.loc[pathway,'Pathway Size']=len(t)
        V.loc[pathway,'N Measured']=t[pthLbl].isin(M[mLbl]).sum()
        if V.loc[pathway,'N Measured']>=minMeasured:
            q=M.loc[M[mLbl].isin(t[pthLbl])]
            origPw=pthKeys['|'.join(q.index.values)]
            if not origPw:
                V.loc[pathway,'Measured Components']='|'.join(q.index.values)
                pthKeys['|'.join(q.index.values)]=pathway
                if computation=='GAGE':
                    V.loc[pathway,'p-value']=computeGAGE(q[sigCol])
                elif computation=='GSEA':
                    pathway_set = {pathway: set(t[pthLbl].unique())}
                    gsea = gp.prerank(rnk=rnk,gene_sets=pathway_set,min_size=minMeasured,permutation_num=5000,seed=42,outdir=None,verbose=False).res2d
                    V.loc[pathway,'p-value']=gsea.at[0,'NOM p-val']
                    V.loc[pathway,'Lead Genes']=gsea.at[0,'Lead_genes']
                V.loc[pathway,'Increase']=((q[st]>0)).sum()
                V.loc[pathway,'No change']=((q[st]==0)).sum()
                V.loc[pathway,'Decrease']=((q[st]<0)).sum()
                V.loc[pathway,'% up']=((q[st]>0)).mean()
                V.loc[pathway,'Increase + P<.05']=((q[st]>0) & (q[sigCol]<.05)).sum()
                V.loc[pathway,'Decrease + P<.05']=((q[st]<0) & (q[sigCol]<.05)).sum()
                V.loc[pathway,'Average Change']=q[st].mean()
                V.loc[pathway,'Median Change']=q[st].median()
            else:
                V.loc[origPw,'Name']=V.loc[origPw,'Name']+'|'+V.loc[pathway,'Name']
                V.loc[origPw,'Subject']=V.loc[origPw,'Subject']+'|'+V.loc[pathway,'Subject']
                V.loc[origPw,'Description']=V.loc[origPw,'Description']+'|'+V.loc[pathway,'Description']
                V.loc[origPw,'pathway ID']=V.loc[origPw,'pathway ID']+'|'+V.loc[pathway,'pathway ID']
                drpList+=[pathway]
        else:
            drpList+=[pathway]
    
    V.drop(drpList,inplace=True)
    # V['N Measured']=V.Decrease+V.Increase+V['No change']
    # V['N Measured']=V.Decrease+V.Increase+V['No change']
    V.sort_values('p-value',inplace=True)
    _, qvals, _, _ = multipletests(V['p-value'], method="fdr_bh")
    V['FDR']=qvals
    
    return(V)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 31 07:31:24 2025

@author: rudy
"""

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


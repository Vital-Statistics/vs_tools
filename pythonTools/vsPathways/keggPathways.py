#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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

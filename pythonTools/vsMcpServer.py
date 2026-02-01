"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
from __future__ import annotations

import os
from typing import Any, Iterable

import pandas as pd
from mcp.server.fastmcp import FastMCP

import vsHighDimensionalData as hdd
import vsMassSpecData as ms
import vsPathways as pathways
import vsTableFormatting as tablefmt

def _env_int(name: str, default: int) -> int:
    try:
        return int(os.getenv(name, str(default)))
    except ValueError:
        return default


_host = os.getenv("MCP_HOST", "127.0.0.1")
_port = _env_int("MCP_PORT", 8000)
_log_level = os.getenv("MCP_LOG_LEVEL", "INFO").upper()
_mount_path = os.getenv("MCP_MOUNT_PATH", "/")
_sse_path = os.getenv("MCP_SSE_PATH", os.getenv("MCP_PATH", "/sse"))
_message_path = os.getenv("MCP_MESSAGE_PATH", "/messages/")

mcp = FastMCP(
    "vs-tools",
    host=_host,
    port=_port,
    log_level=_log_level,
    mount_path=_mount_path,
    sse_path=_sse_path,
    message_path=_message_path,
)


def _df_from_any(data: Any) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data
    if isinstance(data, list):
        return pd.DataFrame.from_records(data)
    if isinstance(data, dict):
        if "records" in data:
            return pd.DataFrame.from_records(data["records"])
        if "data" in data and "columns" in data:
            df = pd.DataFrame(data["data"], columns=data["columns"])
            if "index" in data:
                df.index = data["index"]
            return df
    raise ValueError("Unsupported table format; use list[dict] or dict with records/split keys.")


def _df_to_records(df: pd.DataFrame, reset_index: bool = True) -> list[dict]:
    if reset_index:
        df = df.reset_index()
    return df.to_dict(orient="records")


@mcp.tool()
def tOne_cleanBinary(tbl: Any) -> list[dict]:
    """
    Clean binary variables in a TableOne-style summary table.

    Parameters
    ----------
    tbl : list[dict] | dict
        Table rows as list-of-dicts or pandas split/records dict.

    Returns
    -------
    list[dict]
        Cleaned table as list-of-dicts (index reset).
    """
    df = _df_from_any(tbl)
    out = tablefmt.tOne_cleanBinary(df)
    if out is None:
        return []
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def tOne_toDF(tableone: Any) -> list[dict]:
    """
    Convert a TableOne-style object to a DataFrame with standardized index.

    Parameters
    ----------
    tableone : dict
        TableOne .tableone DataFrame in split/records format.

    Returns
    -------
    list[dict]
        Converted table as list-of-dicts (index reset).
    """
    df = _df_from_any(tableone)
    # Mimic TableOne object structure by wrapping a lightweight container.
    class _TableOne:
        def __init__(self, t):
            self.tableone = t
    out = tablefmt.tOne_toDF(_TableOne(df))
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def imputeMissing(X: Any, k: int = 20, n_iter: int = 20) -> list[dict]:
    """
    Impute missing values via iterative low-rank reconstruction.
    """
    df = _df_from_any(X)
    out = hdd.imputeMissing(df, k=k, n_iter=n_iter)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def regressAll(vIndep: Any, Y: Any, studyVars: Iterable[str] | None = None,
               lbl: str | None = None, showIntercept: bool = False) -> list[dict]:
    """
    Regress each column of Y on vIndep and summarize coefficients and p-values.
    """
    df_x = _df_from_any(vIndep)
    df_y = _df_from_any(Y)
    out = hdd.regressAll(df_x, df_y, studyVars=studyVars, lbl=lbl, showIntercept=showIntercept)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def regressAllNA(vIndep: Any, Y: Any, studyVars: Iterable[str] | None = None,
                 lbl: str | None = None, showIntercept: bool = False) -> list[dict]:
    """
    Regress each column of Y with OLS and logistic models, aggregating p-values.
    """
    df_x = _df_from_any(vIndep)
    df_y = _df_from_any(Y)
    out = hdd.regressAllNA(df_x, df_y, studyVars=studyVars, lbl=lbl, showIntercept=showIntercept)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def testIndexMatch(A: Any, B: Any, showAB: bool = False, showBA: bool = False) -> str:
    """
    Compare index membership between two DataFrames/Series and report differences.
    """
    df_a = _df_from_any(A)
    df_b = _df_from_any(B)
    # Capture prints as a single string.
    import io
    import sys
    buf = io.StringIO()
    old = sys.stdout
    try:
        sys.stdout = buf
        hdd.testIndexMatch(df_a, df_b, showAB=showAB, showBA=showBA)
    finally:
        sys.stdout = old
    return buf.getvalue()


@mcp.tool()
def v_buildSparse(x: list[str], y: list[str], delta: list[float] | None = None) -> list[dict]:
    """
    Build a sparse contingency matrix from categorical indices.
    """
    xc = pd.Categorical(x)
    yc = pd.Categorical(y)
    out = hdd.v_buildSparse(xc, yc, delta=delta)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def parseUniProt(v: dict) -> dict:
    """
    Parse a UniProt JSON record into a flat dict of key fields.
    """
    return pathways.parseUniProt(v)


@mcp.tool()
def apiUniprot(analyte: str, organism_id: str = "9606") -> dict:
    """
    Query UniProt for a gene/protein symbol and return parsed results.
    """
    tbl, raw = pathways.apiUniprot(analyte, organism_id=organism_id)
    return {"table": _df_to_records(tbl, reset_index=True), "raw_results": raw}


@mcp.tool()
def computeGAGE(pv: list[float]) -> float:
    """
    Compute a GAGE-style pathway p-value from analyte p-values.
    """
    return pathways.computeGAGE(pv)


@mcp.tool()
def keggPathways(kegg_id: str) -> list[dict]:
    """
    Fetch KEGG pathways linked to a KEGG entry.
    """
    out = pathways.keggPathways(kegg_id)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def listKeggPathways(org: str = "hsa") -> list[dict]:
    """
    List KEGG pathways for an organism.
    """
    out = pathways.listKeggPathways(org=org)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def pwSMPDB(M: Any, mLbl: str = "HMDB", pthLbl: str = "HMDB ID", sigCol: str = "p-value",
            st: str = "t-statistic", minMeasured: int = 3, dataType: str = "metabolomics",
            pathwayType: list[str] | None = None, dropPathwayList: list[str] | None = None,
            computation: str = "GAGE") -> list[dict]:
    """
    Run SMPDB pathway analysis on a ranked metabolomics/proteomics table.
    """
    df = _df_from_any(M)
    out = pathways.pwSMPDB(
        df,
        mLbl=mLbl,
        pthLbl=pthLbl,
        sigCol=sigCol,
        st=st,
        minMeasured=minMeasured,
        dataType=dataType,
        pathwayType=pathwayType,
        dropPathwayList=dropPathwayList,
        computation=computation,
    )
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def addHMDB(meta: Any) -> list[dict]:
    """
    Add HMDB IDs to a metadata table using a Rosetta Stone mapping.
    """
    df = _df_from_any(meta)
    out = ms.addHMDB(df)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def lookupKeggPathways(keggID: str, reload: bool = False) -> list[str]:
    """
    Look up KEGG pathways for a KEGG compound ID with local caching.
    """
    return ms.lookupKeggPathways(keggID, reload=reload)


@mcp.tool()
def uniprot_to_kegg(uniprot_id: str) -> list[str]:
    """
    Convert a UniProt ID to KEGG IDs with local caching.
    """
    return ms.uniprot_to_kegg(uniprot_id)


@mcp.tool()
def lookup_hmdb(hmdb_id: str) -> dict:
    """
    Look up HMDB details from Metabolomics Workbench.
    """
    return ms.lookup_hmdb(hmdb_id)


@mcp.tool()
def map_names_via_metaboanalyst(names: list[str]) -> dict:
    """
    Map compound names via MetaboAnalyst REST API.
    """
    return ms.map_names_via_metaboanalyst(names)


@mcp.tool()
def lookup_pubchem_by_name(name: str, max_hits: int = 3) -> list[dict]:
    """
    Query PubChem by compound name and return basic annotations.
    """
    return ms.lookup_pubchem_by_name(name, max_hits=max_hits)


@mcp.tool()
def mw_lookup_pubchem(cid: str) -> dict:
    """
    Look up Metabolomics Workbench details by PubChem CID.
    """
    return ms.mw_lookup_pubchem(cid)


@mcp.tool()
def apiLookups(names: list[str], rerunSearch: int = 180) -> list[dict]:
    """
    Perform cached PubChem and Metabolomics Workbench lookups for analytes.
    """
    s = pd.Series(names)
    out = ms.apiLookups(s, rerunSearch=rerunSearch)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def loadMetabolomics(filePath: str, lastSampCol: str = "Injection Number",
                     sid: str = "Customer Sample Identification", maxMissing: int = 1,
                     header: int = 1, lodCol: str | None = None, nHeader: int | None = None,
                     log2Transform: bool | None = None, normalize: bool = False,
                     platform: str | None = None, stdFilter: float | bool = False,
                     imputeMissing: bool = True) -> dict:
    """
    Load metabolomics data from a vendor export and perform basic cleaning.
    """
    ss, aa, t = ms.loadMetabolomics(
        filePath,
        lastSampCol=lastSampCol,
        sid=sid,
        maxMissing=maxMissing,
        header=header,
        lodCol=lodCol,
        nHeader=nHeader,
        log2Transform=log2Transform,
        normalize=normalize,
        platform=platform,
        stdFilter=stdFilter,
        imputeMissing=imputeMissing,
    )
    return {
        "samples": _df_to_records(ss, reset_index=True) if ss is not None else [],
        "analytes": _df_to_records(aa, reset_index=True) if aa is not None else [],
        "table": _df_to_records(t, reset_index=True) if t is not None else [],
    }


@mcp.tool()
def loadSynonyms(meta: Any) -> list[dict]:
    """
    Load synonym mappings and join with metadata.
    """
    df = _df_from_any(meta)
    out = ms.loadSynonyms(df)
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def t_tg(s: str) -> str:
    """Normalize triglyceride naming format."""
    return ms.t_tg(s)


@mcp.tool()
def t_lpc(s: str) -> str:
    """Normalize LPC naming format."""
    return ms.t_lpc(s)


@mcp.tool()
def t_PC(s: str) -> str:
    """Normalize PC naming format."""
    return ms.t_PC(s)


@mcp.tool()
def t_PCO(s: str) -> str:
    """Normalize PC-O naming format."""
    return ms.t_PCO(s)


@mcp.tool()
def t_Hex2(s: str) -> str:
    """Normalize Hex2Cer naming format."""
    return ms.t_Hex2(s)


@mcp.tool()
def t_Hex3(s: str) -> str:
    """Normalize Hex3Cer naming format."""
    return ms.t_Hex3(s)


@mcp.tool()
def t_Hex(s: str) -> str:
    """Normalize HexCer naming format."""
    return ms.t_Hex(s)


@mcp.tool()
def t_DG(s: str) -> str:
    """Normalize DG naming format."""
    return ms.t_DG(s)


@mcp.tool()
def t_Cer(s: str) -> str:
    """Normalize Cer naming format."""
    return ms.t_Cer(s)


@mcp.tool()
def t_CE(s: str) -> str:
    """Normalize CE naming format."""
    return ms.t_CE(s)


@mcp.tool()
def t_simple(s: str) -> str:
    """Normalize simple lipid naming formats."""
    return ms.t_simple(s)


@mcp.tool()
def t_PLA(s: str) -> str:
    """Normalize PLA2 naming format."""
    return ms.t_PLA(s)


@mcp.tool()
def m_canonicalNames(res: Any, colName: str = "oName", metaAppend: list[str] | None = None) -> list[dict]:
    """
    Standardize analyte names and optionally append metadata.
    """
    df = _df_from_any(res)
    out = ms.m_canonicalNames(df, colName=colName, metaAppend=metaAppend or [])
    return _df_to_records(out, reset_index=True)


@mcp.tool()
def queryMetWorkbench(mwl: list[str]) -> list[dict]:
    """
    Query Metabolomics Workbench RefMet matching and details for names.
    """
    out = ms.queryMetWorkbench(mwl)
    return _df_to_records(out, reset_index=True)


def main() -> None:
    transport = os.getenv("MCP_TRANSPORT", "sse")
    mount_path = os.getenv("MCP_MOUNT_PATH", None)
    if transport == "stdio":
        mcp.run()
        return
    mcp.run(transport=transport, mount_path=mount_path)


if __name__ == "__main__":
    main()

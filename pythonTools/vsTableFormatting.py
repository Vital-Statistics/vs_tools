"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
def tOne_cleanBinary(tbl):
    """
    Clean binary variables in a TableOne-style summary table.

    This converts two-row binary variables (e.g., 0/1, True/False, Y/N) into a
    single row while preserving the missing count and (if present) the P-Value.
    The function expects a DataFrame with 'Variable' and 'Value' columns; if
    needed it will attempt a reset_index to surface those columns.

    Parameters
    ----------
    tbl : pandas.DataFrame
        TableOne-style summary table (or output of `tOne_toDF`) with
        columns including 'Variable', 'Value', and 'Missing'. May also include
        'P-Value'.  This object is what is returned by the tOne_toDF function.

    Returns
    -------
    pandas.DataFrame
        Cleaned table indexed by ['Variable', 'Value'] with binary variables
        collapsed to a single row.

    Raises
    ------
    None
        Errors are reported via print statements and the function returns None.
    """
    import pandas as pd
    import math
    if 'Variable' not in list(tbl):
        tbl=tbl.reset_index()
    if 'Variable' not in list(tbl):
        print('Unable to find "Variable" column. Perhaps run tOne_toDF to convert default tableone output?')
        return
    if 'Value' not in list(tbl):
        print('Unable to find "Value" column.')
        return
    
    addRows=list()
    for rw in tbl.Variable.unique():
        t=tbl.loc[tbl.Variable==rw]
        if (len(t)==2) and (set(t.Value) in [{'1.0','0.0'},{'True','False'},{'0','1'},{'Y','N'}]):
            nMissing=[v for v in t.Missing if v!='']
            if 'P-Value' in list(t):
                pv=[v for v in t['P-Value'] if v!='']
                if len(pv)!=1:
                    print('P-Values are strange for variable: '+rw)
                    return
                else:
                    pv=pv[0]
            if len(nMissing)!=1:
                print('Missing counts are strange for variable: '+rw)
                return
            else:
                nMissing=nMissing[0]
            newRow=t.loc[t.Value.isin(['1.0','True','1','Y'])].copy()
            newRow.Missing=nMissing
            if 'P-Value' in list(t):
                newRow['P-Value']=pv
            newRow.Value=''
            addRows+=[newRow]
            tbl.drop(t.index.values,inplace=True)
    tbl=pd.concat([tbl]+addRows)
    tbl=tbl.set_index(['Variable','Value'])
    return(tbl)


def tOne_toDF(tbl):
    """
    Convert a TableOne object to a DataFrame with a standardized index.

    Parameters
    ----------
    tbl : tableone.TableOne
        TableOne object whose `.tableone` attribute is a MultiIndex DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by ['Variable', 'Value'] with flattened columns.
    """
    tbl=tbl.tableone
    tbl.columns=[v[1] for v in list(tbl)]
    tbl.reset_index(inplace=True)
    tbl.rename(columns={'level_0':'Variable','level_1':'Value'},inplace=True)
    return(tbl.set_index(['Variable','Value']))

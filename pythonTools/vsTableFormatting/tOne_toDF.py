#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 15:32:38 2023

@author: rudy
"""

def tOne_toDF(tbl):
    tbl=tbl.tableone
    tbl.columns=[v[1] for v in list(tbl)]
    tbl.reset_index(inplace=True)
    tbl.rename(columns={'level_0':'Variable','level_1':'Value'},inplace=True)
    return(tbl.set_index(['Variable','Value']))
    
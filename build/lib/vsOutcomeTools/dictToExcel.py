#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 09:02:53 2023

@author: rudy
"""

def dictToExcel(d,outFile):
    with pd.ExcelWriter(outFile, engine='xlsxwriter') as writer:
        for sheet_name, df in d.items():
            df.to_excel(writer, sheet_name=sheet_name)
        # writer.save()

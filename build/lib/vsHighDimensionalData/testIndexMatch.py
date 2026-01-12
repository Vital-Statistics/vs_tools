#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 05:56:24 2023

@author: rudy
"""

def testIndexMatch(A,B,showAB=False,showBA=False):
    a=pd.Series(A.index.values)
    b=pd.Series(B.index.values)
    for k,t in {'A':a,'B':b}.items():
        n=len(t)
        u=len(t.unique())
        print('Data set '+k+': '+str(n)+' rows with '+str(u)+' unique indices.')
        if u!=n:
            print(t.loc[t.duplicated(keep=False)].value_counts())
        if t.isna().sum()>0:
            print('Data set '+k+' has '+str(t.isna().sum())+' NaNs in its index')
            print('')
    
    ab=set(a.unique())-set(b.unique())
    ba=set(b.unique())-set(a.unique())
    union=set(b.unique()) & set(a.unique())
    print()
    print('Data set A has '+str(len(ab))+' unique elements that are not in B.')
    if showAB:
        print('Elements in A that are not in B:')
        print(', '.join(list(ab)))
    print()
    print('Data set B has '+str(len(ba))+' unique elements that are not in A.')
    if showBA:
        print('Elements in B that are not in A:')
        print(', '.join(list(ba)))
    print('There are '+str(len(union))+ ' unique elements in both indices.')
    
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 15:18:36 2017

@author: jel2
"""

import datetime
from collections import Counter
import numpy as np
import pandas as pd
import scipy
import operator

def allValsList(t,pth,res=None):
    if(res==None):
        res=list()
    if (type(t)==list) | (type(t)==set):
        for s in t:
            allValsList(s,pth,res)
    elif len(pth)==0:
        res.append(t)
    elif pth[0] in t:
        allValsList(t[pth[0]],pth[1:],res)
    return(res)

def allValsSet(t,pth,res=None):
    if(res==None):
        res=set()
    if (type(t)==list) | (type(t)==set):
        for s in t:
            allValsSet(s,pth,res)
    elif len(pth)==0:
        res.update([t])
    elif pth[0] in t:
        allValsSet(t[pth[0]],pth[1:],res)
    return(res)

###### failing.  not sure how to fix
#def allValsGen(t,pth):
#    if type(t)==list:
#        for s in t:
#            yield allValsGen(s,pth)
#    elif len(pth)==0:
#        yield t
#    elif pth[0] in t:
#        yield allValsGen(t[pth[0]],pth[1:])

def detailsToDF(pts):
    for t in pts:
        t['EVENTS']=pd.DataFrame.from_records(t['EVENTS'],index='EVENT_TIME')

#def dateFilter(pts,wdw=datetime.timedelta(-180),st=None):
#    if type(wdw)==int:
#        wdw=datetime.timedelta(days=wdw)
#    if type(st)==int:
#        st=datetime.timedelta(days=st)
#    if st is None:
#        st=[x['STUDY_START'] for x in pts]
#    mxd=datetime.timedelta(0)
#    i=0
#    while True:
#        delta=pts[i]['EVENTS'].index-st[i]
#        yield pts[i]['EVENTS'][(min(wdw,mxd)<=delta) & (delta<max(wdw,mxd))]
#        i+=1

def computeSparse(pts,pth,src=None,dmn=datetime.datetime(1,1,1),dmx=datetime.datetime(9999,1,1)):
    if type(src)=='str':
        src=set([src])
    ########## sparse matrices computed only from the events dictionary
    if pth[0]=='EVENTS':
        pth.pop(0)
    #### compute the dictionary of elements contained in path across all patients
    #### this creates a label for every type of source.  We really want to restrict
    #### this to the sources in src.  How to do?
    lbl=allValsSet(pts,['EVENTS']+pth)
    lbl=list(lbl)
    lbl=dict((value, idx) for idx,value in enumerate(lbl))
    
    drc=list()
    for i,s in enumerate(pts):
        if src is not None:
            try:
                t=[x for x in s['EVENTS'] if x['SOURCE'] in src and x['EVENT_TIME'] is not None and x['EVENT_TIME']>=dmn[i] and x['EVENT_TIME']<dmx[i]]
            except KeyError:
                continue
        else:
            t=s['EVENTS']
        U=allValsList(t,pth)
        if len(U)>0:
            ix=list()
            for u in U:
                ix.append(lbl[u])
            a,b=zip(*Counter(ix).items())
            drc+=zip(ones(i),a,b)
        if i%100==0:
            print(i)
    a,b,c=zip(*drc)
    M=scipy.sparse.coo_matrix( (np.array(c),(np.array(a),np.array(b))),shape=(len(pts),len(lbl)) )
    M=M.tocsr()
    sm=np.array(M.sum(0))[0]
    M=M[:,sm>0]
    lbl=list(lbl.items())
    lbl.sort(key=operator.itemgetter(1))
    lbl=[s[0] for s,t in zip(lbl,sm) if t>0]
    
    return lbl, M


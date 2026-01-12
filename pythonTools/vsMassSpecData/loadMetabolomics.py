# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 14:48:04 2018

@author: joseph lucas
"""

import re

def loadMetabolomics(filePath,lastSampCol='Injection Number',sid='Customer Sample Identification',maxMissing=1
                     ,header=1,lodCol=None,nHeader=None,log2Transform=None,normalize=False,platform=None,stdFilter=False
                     ,imputeMissing=True):
    import pandas as pd
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


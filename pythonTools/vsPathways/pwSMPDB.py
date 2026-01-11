# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 10:05:56 2021

@author: joest
"""

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

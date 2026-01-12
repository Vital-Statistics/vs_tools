# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 12:10:57 2018

@author: gsl8
"""

###############################################################################
def columns(df,d,i,out):

    import numpy as np

    if df[i].dtype is np.dtype('O'):

        c = {x: {d: {} for d in sorted(d.keys())} for x in sorted(set(df.loc[~df[i].isnull(),i]))}

        for c_key in sorted(c.keys()):
            a=[]
            for d_key in sorted(d.keys()):
                c[c_key][d_key]['num'] = sum(d[d_key]['df'][i]==c_key)
                c[c_key][d_key]['perc'] = '{:.1f}'.format(round(c[c_key][d_key]['num']/d[d_key]['num'],3)*100)
                a.append(str(c[c_key][d_key]['num'])+' ('+str(c[c_key][d_key]['perc'])+')')
            out[str(c_key.capitalize())+' - n.(%)'] = a

    elif df[i].dtype is np.dtype('bool'):

        c = {True: {d: {} for d in sorted(d.keys())}}
        a=[]
        for d_key in sorted(d.keys()):
            c[True][d_key]['num'] = sum(d[d_key]['df'][i]==True)
            c[True][d_key]['perc'] = '{:.1f}'.format(round(c[True][d_key]['num']/d[d_key]['num'],3)*100)
            a.append(str(c[True][d_key]['num'])+' ('+str(c[True][d_key]['perc'])+')')
        out[i.capitalize()+' - n.(%)'] = a

    elif df[i].dtype is np.dtype('float'):

        c = {i:{d:{} for d in sorted(d.keys())}}
        a=[]
        for d_key in sorted(d.keys()):
            c[i][d_key]['avg'] = round(np.mean(d[d_key]['df'][i]),1)
            c[i][d_key]['sd'] = round(np.std(d[d_key]['df'][i]),1)
            a.append(str(c[i][d_key]['avg'])+' ('+str(c[i][d_key]['sd'])+')')
        out[str(i.capitalize())+' - avg.(sd.)'] = a

###############################################################################
def demTable(df, splitColumn=None, output=None, filePath=None):

    '''
    Creates 'Table 1' or demographics table.

    df = Pandas dataframe
    splitColumn = Column used to group variables.
    output = If none, output is dataframe.  Other formats: latex, html, and csv.
    '''

    import numpy as np
    import pandas as pd
    import datetime

    if not splitColumn:
        df['ukjdsfg']=''
        splitColumn='ukjdsfg'

    pd.options.display.max_colwidth=200
    out = pd.DataFrame()

    d = {x: {'df':pd.DataFrame()} for x in sorted(set(df[splitColumn]))}

    a = []
    for key in sorted(d.keys()):
        d[key]['df'] = df.loc[df[splitColumn]==key,:] # populate each dataframe as subset of df
        d[key]['num'] = d[key]['df'].shape[0]           # N of each group
        a.append('N = '+str(d[key]['num']))

    out['Chacteristics'] = a # display N of each group


    for i in sorted(list(df.drop(splitColumn,axis=1))):
        columns(df,d,i,out=out)

    out = out.transpose()

    out_col = [x for x in sorted(d.keys())]

    for i in range(0,len(d.keys())):
        out.rename(columns={i:str(out_col[i].capitalize())},inplace=True)

    if output=='csv':
        out.to_csv(filePath)
    elif output=='latex':
        return print(out.to_latex(escape=False))
    elif output=='html':
        return print(out.to_html(escape=False))
    elif output==None:
        return out
###############################################################################
#### Example:
#test = eList[eList['SOURCE']=='MEDICATION']
#test = test.merge(pts, left_on='PATIENT_ID', right_index=True)
#test['DRUG'] = 'NEITHER'
#test.loc[test['LABEL'].str.contains(re.compile('metformin',re.IGNORECASE),re.IGNORECASE),'DRUG'] = 'METFORMIN'
#test.loc[test['LABEL'].str.contains(re.compile('sulf',re.IGNORECASE),re.IGNORECASE),'DRUG'] = 'SULFONYLUREA'
#Counter(test['DRUG'])
#test['DEATH'] = ~test['DEATH_DATE'].isna()
#Counter(test['DEATH'])
## This does not actually calculate age, just using this as an example of numerical data.
#test['AGE'] = (test['EVENT_TIME'] - test['DOB'])/datetime.timedelta(days=365.25)
#test = test[['DRUG','SEX','AGE','DEATH']]
#
#demTable(df=test,splitColumn='DRUG',output=None)
#demTable(df=test,splitColumn='DRUG',output='latex')
#demTable(df=test,splitColumn='DRUG',output='html')
#demTable(df=test,splitColumn='DRUG',output='csv',filePath='/Users/samlusk/Desktop/test.csv')

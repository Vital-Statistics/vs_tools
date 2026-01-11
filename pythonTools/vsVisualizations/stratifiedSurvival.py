# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 08:07:47 2018

@author: jel2
"""

def stratifiedSurvival(t,eventTime,eventIndicator=None,followupTime=None,group=None):
    
    import matplotlib.pyplot as plt
    import lifelines as lf
    from lifelines.plotting import add_at_risk_counts
    import pandas as pd
    import copy

    tm=t[eventTime].copy()
    
    if(group is None):
        grp=pd.Series('Population',index=t.index)
    else:
        grp=t[group]
    
    if(eventIndicator is None):
        ev=~t[eventTime].isnull()
        tm[tm.isnull()]=t.loc[tm.isnull(),followupTime]
        
    ######### Kaplan Meier curves stratified by sex
    kl=list()
    kmf = lf.KaplanMeierFitter()
    fig,ax=plt.subplots()
    for g in set(grp):
        kmf.fit(tm[grp==g],ev[grp==g],label=g)
        kmf.plot(ax=ax)
        kl.append(copy.deepcopy(kmf))
        
    add_at_risk_counts(*kl, ax=ax)
        
    plt.legend(loc='lower left')
    plt.ylim([0,1])
    plt.xlabel('Time (years)')
    plt.ylabel('Survival')
    plt.title('Kaplan-Meier survival curve')
#    add_at_risk_counts(kmf1,kmf2, ax=ax)


#ax.spines['bottom'].set_position(('axes', -0.15 * 6.0 / fig.get_figheight()))

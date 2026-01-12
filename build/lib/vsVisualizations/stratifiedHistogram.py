# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 08:07:47 2018

@author: jel2
"""

def stratifiedHistogram(t,x,group=None,bins=10):
    
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    y=t[x]
    if(group is None):
        grp=pd.Series('Population',index=t.index)
    else:
        grp=t[group]
    
    b=np.arange(min(y),max(y),(max(y)-min(y))/bins)
    clr=[(1,0,0,.3),(0,1,0,.3),(0,0,1,.3)]
    for i,g in enumerate(set(grp)):
        plt.hist(y[grp==g],bins=b,label=g,normed=True,color=clr[i])
        
    plt.legend()
    plt.xlabel(x)
    plt.ylabel('Probability')
#    plt.title('Age distribution')



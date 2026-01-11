# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:14:51 2018

@author: jel2
"""

def lsReg(M,params):
    import lifelines as lf
    
    cph = lf.CoxPHFitter()
    return(cph.fit(M, duration_col=params['tmCol'], event_col=params['ixCol']))
        
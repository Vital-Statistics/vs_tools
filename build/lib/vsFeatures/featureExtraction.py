# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:10:49 2018

@author: jel2
"""

def featureExtraction(eList,function, X=None, params=None):
    return X.join(function(eList,params),how='left').fillna(0)
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:18:44 2022

@author: joest
"""

def rhoToZ(rho):
    import math
    return(.5*math.log((1+rho)/(1-rho)))

# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:28:12 2022

@author: joest
"""

############ Reference:
###########  Hinkle DE, Wiersma W, Jurs SG (1988) Applied statistics for the behavioral sciences. 2nd ed. Boston: Houghton Mifflin Company
############ Fisher r to z plus t-test

from vsUtilityFunctions.rhoToZ import rhoToZ

def twoCorrTest(rho0,n0,rho1,n1):
    import math
    from scipy.stats import norm
    z=(rhoToZ(rho0)-rhoToZ(rho1))/math.sqrt(1/(n0-3)+1/(n1-3))
    p=(1-norm.cdf(abs(z)))*2
    # p=norm.logcdf(abs(z))
    return((z,p))



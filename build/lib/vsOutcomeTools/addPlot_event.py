# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 16:17:05 2017

@author: jel2
"""
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

def cBands_poisGamma(u,alpha=.5,beta=.5,ivl=.95):
    a=u['Events']+alpha
    b=u['N']+beta
    delta=(1-ivl)/2
    u['low']=stats.gamma.ppf(delta,a,scale=1/b)
    u['hi']=stats.gamma.ppf(1-delta,a,scale=1/b)
    return(u)

def addPlot_event(fig,w,n,e,clr='k',lbl=None):
    u=pd.DataFrame({'week':w, 'N':n, 'Events':e})
    u=u.sort_values('week')
    u['N']=u['N'].cumsum()
    u['Events']=u['Events'].cumsum()
    u=u.apply(cBands_poisGamma,axis=1)
    plt.plot(u['week'],u['Events']/u['N'],color=clr,label=lbl)
    x=list(u.week)+list(u.week)[::-1]
    y=list(u.low)+list(u.hi)[::-1]
    fig.gca().fill(x,y,facecolor=clr,edgecolor=None,alpha=.3)

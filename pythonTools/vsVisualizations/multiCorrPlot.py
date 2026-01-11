# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:08:42 2022

@author: joest
"""

def compSig(dat):
    dat=dat.loc[dat.isna().sum(axis=1)==0]
    dat=np.array(dat)
    mn=dat.mean(axis=0)
    dat=dat-mn
    res=np.array(dat.transpose())@np.array(dat)/len(dat)
    return(res)

def mvnOval(a,b,color=None,scale=2.65,alpha=.2):
    from scipy.linalg import sqrtm
    t=pd.concat([a,b],axis=1)
    mn=np.array(t).mean(axis=0)
    sig=scale*compSig(t)
    x=np.array([[math.cos(theta),math.sin(theta)] for theta in np.arange(0,2*math.pi,2*math.pi/40)])
    res=x@sqrtm(sig+1e-6*np.identity(2))
    if color is not None:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],color=color,edgecolor=None,alpha=alpha)
    else:
        plt.fill(res[:,0]+mn[0],res[:,1]+mn[1],edgecolor=None,alpha=alpha)
    return res

def multiCorrPlot(u,x,y,grp,xlbl='',ylbl=''):
    u=u.loc[u[[x,y]].isna().sum(axis=1)==0]
    nclr=len(u[grp].unique())-1
    for i,gg in enumerate(u[grp].unique()):
        plt.scatter(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],label=gg,color=vs_clr(nclr-i))
        mvnOval(u.loc[u[grp]==gg,x],u.loc[u[grp]==gg,y],color=vs_clr(nclr-i,alpha=.2))
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.legend()
    plt.tight_layout()


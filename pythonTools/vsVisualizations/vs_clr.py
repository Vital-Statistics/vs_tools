# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:47:06 2021

@author: JoeLucas
"""

def vs_clr(cc='blue',alpha=1):
    r={'blue':(15/256,98/256,160/256)
       ,'dark blue':(48/256,48/256,95/256)
       ,'grey':(182/256,183/256,182/256)
       ,'gray':(182/256,183/256,182/256)
       ,'dark grey':(84/256,85/256,84/256)
       ,'dark gray':(84/256,85/256,84/256)
       }
    if type(cc)==int:
        cc=['blue','dark gray','gray'][cc]
    res=r[cc]
    if alpha<1:
        res=tuple(list(res)+[alpha])
    return(res)

def vs_clr_V1(cc='blue',alpha=1):
    r=(15/256,98/256,160/256)
    if(cc=='dark blue'):
        r=(48/256,48/256,95/256)
    elif((cc=='grey') | (cc=='gray')):
        r=(182/256,183/256,182/256)
    elif((cc=='dark grey') | (cc=='dark gray')):
        r=(84/256,85/256,84/256)
    if alpha<1:
        r=tuple(list(r)+[alpha])
    return(r)

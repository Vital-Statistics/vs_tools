# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 12:16:57 2021

@author: JoeLucas
"""

def loopProgress(i,N,lbl='',nSteps=20):
    import math
    i=i+1
    if i==1:
        print(lbl+': |'+'-'*i+' '*(nSteps-i)+'|',end='')
    i=math.floor(i*nSteps/N)
    N=nSteps
    if i<N:
        print('\r'+lbl+': |'+'-'*i+' '*(N-i)+'|',end='')
    else:
        print('\r'+lbl+': |'+'-'*i+' '*(N-i)+'|')
    
    
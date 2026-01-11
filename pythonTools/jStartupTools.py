# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 08:10:22 2016

@author: lucas
"""

import os
from glob import glob
import re
import pandas as pd
pd.set_option('display.max_columns', 500)

def vs_palette():
    return(['#0e619f', '#b22222', '#006400', '#b5b6b5', '#2f2f5e','#e57373', '#98fb98', '#535453'])

def cloneAnalysis(pth='.',root='START_PATH'):
    import shutil
    import datetime
    os.chdir(pth)
    nm=os.getcwd().split(os.sep)[-1]
    d,nm=nm.split('-')
    newD=datetime.datetime.today().strftime('%Y%m%d')
    if os.path.isdir('../'+newD+'-'+nm):
        newD=datetime.datetime.today().strftime('%Y%m%d%H%M')
    np='../'+newD+'-'+nm+'/'
    os.makedirs(np)
    ## Identify files and subdirectories that aren't in results
    oList = [d for d in glob('*', recursive=True) if 'results' not in d.split(os.path.sep)]
    for ob in oList:
        if os.path.isdir(ob): shutil.copytree(ob,np+ob)
        else: shutil.copy2(ob,np+ob)
    # fList=glob('*.py')
    # for fl in fList:
    #     shutil.copy2(fl, np+fl)
    #     if os.path.isdir('localTools'):
    #         shutil.copytree('localTools',np+'localTools')
    os.chdir(np)
    
    np='/'.join(os.getcwd().split(os.sep)).replace(os.environ[root],"os.environ["+root+"]+'")+"'"
    # nl='os.chdir('+np+')'
    for fl in glob('**/*.py',recursive=True):
        file = open(fl, "r")
        newText=''
        for line in file:
            if re.match(r"^\s*PROJECT_PATH\s*=", line):  ### could modify to exchange any collection of parameters by passing dictionary
                line=line.split('=')[0]+'='+np+'\n'
            newText+=line
        # for line in file:
        #     if line[:9]=='os.chdir(':
        #         newText+=nl+'\n'
        #     else:
        #         newText+=line
        file.close()
        # opening the file in write mode
        fout = open(fl, "w")
        fout.write(newText)
        fout.close()

def newAnalysis(pth,resFolder='analysis',fileName='analysis',outputFolder='results'):
    import datetime
    
    fld=setFolder(pth+'/'+datetime.date.today().strftime("%Y%m%d")+'-'+resFolder)
    os.chdir(fld)
    fld=os.getcwd()
    fld='/'.join(fld.split(os.sep))
    # if fileName is None:
    #     fileName=fld.replace(os.environ['START_PATH'],'')
    #     if fileName==fld:
    #         print('Not in start_path.  Setting filename to "analysis".')
    #         fileName='analysis'
    #     else:
    #         fileName=fileName.split('/')
        
    a="os.environ['START_PATH']+"+"'"+fld.replace(os.environ['START_PATH'],'')+"'"
    fh=open(fileName+'.py','w',encoding='utf8')
    fileText="""
# created by newAnalysis script on '+datetime.date.today().strftime("%Y/%m/%d")
import math

os.chdir(<a>)

reGit()
reTool()
reLoc()

fld=setFolder('<outputFolder>')


"""
    fh.write(fileText.replace('<outputFolder>',outputFolder).replace('<a>',a))
    fh.close()
    if not os.path.exists('localTools'):
        os.makedirs('localTools')
    return '%edit '+fileName +'.py'

def setFolder(folder):
    from pathlib import Path
    Path(folder).mkdir(parents=True, exist_ok=True)
    if folder[-1]!='/]':
        folder+='/'
    return(folder)

def newLocalFunction(nm,*args):
    if not os.path.exists('localTools'):
        os.makedirs('localTools')
    if not os.path.isfile('localTools'+nm+'.py'):
        fh=open('localTools/'+nm+'.py','w',encoding='utf8')
    else:
        fh=open('localTools/'+nm+'.py','a',encoding='utf8')
    if(len(args)==0):
        fh.write('\ndef '+nm+'():\n')
    else:
        for fct in args:
            fh.write('\ndef '+fct+'()\n')
    fh.close()
    return '%edit localTools/'+nm+'.py'

def reLoad(fld):
    from glob import glob
    for fl in glob(fld+'/**/*.py',recursive=True):
        print(fl)
        exec(open(fl).read(),globals())
    return True

def reLoc():
    from glob import glob
    for fl in glob('localTools/*.py'):
        print(fl)
        exec(open(fl).read(),globals())
    return True

def reTool(folder=None):
    from glob import glob
    if folder is None:
        pth=os.environ['UTIL_PATH']+'/**/*.py'
    else:
        pth=os.environ['UTIL_PATH']+'/'+folder+'/*.py'
    for fl in glob(pth,recursive=True):
        print(fl)
        exec(open(fl).read(),globals())
    return True

def reGit(folder=None):
    from glob import glob
    if folder is None:
        pth=os.environ['GIT_PATH']+'/**/*.py'
    else:
        pth=os.environ['GIT_PATH']+'/'+folder+'/*.py'
    for fl in [v for v in glob(pth,recursive=True) if not re.search('Example',v)]:
        print(fl)
        exec(open(fl).read(),globals())
    return True

def savePKL(obj,pth):
    import pickle
    with open(pth, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def readPKL(pth):
    import pickle
    with open(pth, 'rb') as handle:
        b = pickle.load(handle)
    return(b)

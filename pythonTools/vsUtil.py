# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 08:10:22 2016

@author: Vital Statistics, LLC.
"""

import os
from glob import glob
import re
import pandas as pd
pd.set_option('display.max_columns', 500)

def vs_palette():
    """
    Return the default Vital Statistics color palette.

    Returns
    -------
    list[str]
        Hex color codes.
    """
    return(['#0e619f', '#b22222', '#006400', '#b5b6b5', '#2f2f5e','#e57373', '#98fb98', '#535453'])

def vs_clr(cc='blue',alpha=1):
    """
    Map a named color (or index) to an RGB/RGBA tuple.

    Parameters
    ----------
    cc : str | int, default 'blue'
        Color name or index (0..2 maps to blue/dark gray/gray).
    alpha : float, default 1
        Alpha channel for RGBA output.

    Returns
    -------
    tuple
        RGB or RGBA tuple with values in [0, 1].
    """
    palette = vs_palette()
    name_to_hex = {
        'blue': palette[0],
        'red': palette[1],
        'dark red': palette[1],
        'green': palette[2],
        'dark green': palette[2],
        'gray': palette[3],
        'grey': palette[3],
        'dark blue': palette[4],
        'light red': palette[5],
        'light green': palette[6],
        'dark gray': palette[7],
        'dark grey': palette[7],
    }
    if isinstance(cc, int):
        hex_color = palette[cc]
    else:
        key = str(cc).strip().lower()
        hex_color = name_to_hex[key]
    res = tuple(int(hex_color[i:i+2], 16)/255 for i in (1, 3, 5))
    if alpha < 1:
        res = tuple(list(res) + [alpha])
    return res

def cloneAnalysis(pth='.',root='START_PATH'):
    """
    Clone the current analysis folder into a date-stamped sibling directory.

    Parameters
    ----------
    pth : str, default '.'
        Path to the analysis folder to clone.
    root : str, default 'START_PATH'
        Environment variable name for the root path used to rewrite PROJECT_PATH.

    Returns
    -------
    None
    """
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
    """
    Create a new analysis starter script and local tools folder. If the path doesn't exist, create it.

    Parameters
    ----------
    pth : str
        Base directory in which to create the analysis folder.
    resFolder : str, default 'analysis'
        Base name for the dated analysis folder.
    fileName : str, default 'analysis'
        Name of the starter script (without extension).
    outputFolder : str, default 'results'
        Output subfolder name.
    """
    import datetime
    
    fld=setFolder(pth+'/'+datetime.date.today().strftime("%Y%m%d")+'-'+resFolder)
    os.chdir(fld)
    fld=os.getcwd()
    fld='/'.join(fld.split(os.sep))
    a="os.environ['START_PATH']+"+"'"+fld.replace(os.environ['START_PATH'],'')+"'"
    fh=open(fileName+'.py','w',encoding='utf8')
    fileText="""
# created by newAnalysis script on '+datetime.date.today().strftime("%Y/%m/%d")
# @author: Vital Statistics, LLC
# Copyright (c) 2026 Vital Statistics, LLC

from vsUtil import *

os.chdir(<a>)
reLoc()

fld=setFolder('<outputFolder>')

"""
    fh.write(fileText.replace('<outputFolder>',outputFolder).replace('<a>',a))
    fh.close()
    if not os.path.exists('localTools'):
        os.makedirs('localTools')

def setFolder(folder):
    """
    Ensure a folder exists and return a normalized path.

    Parameters
    ----------
    folder : str
        Folder path to create.

    Returns
    -------
    str
        Folder path with a trailing slash.
    """
    from pathlib import Path
    Path(folder).mkdir(parents=True, exist_ok=True)
    if folder[-1]!='/]':
        folder+='/'
    return(folder)

def reLoad(fld):
    """
    Reload all Python files matching a glob into the current globals.

    Parameters
    ----------
    fld : str
        Glob pattern for Python files to load.

    Returns
    -------
    bool
        True on completion.
    """
    from glob import glob
    for fl in glob(fld+'/**/*.py',recursive=True):
        print(fl)
        exec(open(fl).read(),globals())
    return True

def reLoc():
    """
    Reload all Python files from the localTools folder.

    Returns
    -------
    bool
        True on completion.
    """
    return reLoad('localTools/*.py')

def savePKL(obj,pth):
    """
    Save an object to disk using pickle.

    Parameters
    ----------
    obj : any
        Object to serialize.
    pth : str
        File path for output.
    """
    import pickle
    with open(pth, 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def readPKL(pth):
    """
    Load a pickled object from disk.

    Parameters
    ----------
    pth : str
        File path to the pickle file.

    Returns
    -------
    any
        Deserialized object.
    """
    import pickle
    with open(pth, 'rb') as handle:
        b = pickle.load(handle)
    return(b)

"""
@author: Vital Statistics, LLC
Copyright (c) 2026 Vital Statistics, LLC
"""
# Auto-generated from vsJavaScript/ package modules.

# ---- source: vsJavaScript/data2js.py ----


def data2js(data,dataName,outFile):
	### Load a csv into 'data'.  dataName is the name of the javascript package, and outFile is the location to write the javascript
    import pandas as pd

    data=data.applymap(str)
    header="""(function (global, factory) {
		typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
		typeof define === 'function' && define.amd ? define(['exports'], factory) :
		(factory((global.<dataName> = global.<dataName> || {})));
	}(this, (function (exports) { 'use strict';
	var data=[""".replace('<dataName>',dataName)

    footer="""   ]
	exports.data=data;
	Object.defineProperty(exports, '__esModule', { value: true });
	})));"""

    cols=list(data)    

    f = open(outFile, "w")
    f.write(header)

    ot='{"'+cols[0]+'":"'+data[cols[0]].str.strip().replace('"',"'")+'"'
    for lbl in cols[1:]:
        ot+=',"'+lbl+'":"'+data[lbl].str.replace('"',"'")+'"'
    ot+='}'
    ot=ot.str.cat(sep=',\n')
    f.write(ot)
    f.write(footer)
    f.close()

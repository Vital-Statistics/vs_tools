#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  4 11:51:49 2026

@author: rudy
"""

import pandas as pd

from vsVisualizations.vs_palette import vs_palette

def boxSwarm(group,value):
    import seaborn as sns
    df=pd.DataFrame({'value':value,'group':group})

    ax=sns.boxplot(data=df, x="group", y="value",hue="group",showfliers=False, width=0.6,palette=vs_palette())
    sns.swarmplot(data=df, x="group", y="value",hue="group", size=4, alpha=0.7,palette=vs_palette())
    # ax.legend_.remove()

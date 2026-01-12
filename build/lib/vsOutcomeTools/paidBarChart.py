
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def paidBarChart(tbl,maxPlot=10):
    
    if(len(tbl)>maxPlot):
        kp=tbl.iloc[:maxPlot,:]
    kp=set(tbl.index)
    q=pd.DataFrame(tbl[['After','Transition','Before']].stack().rename('cst')).reset_index()
    q=q[q['grp'].apply(lambda x: x in kp)]
    
    fig, ax = plt.subplots(figsize=(8, round(maxPlot/2)+1), dpi=80)
    fig=sns.barplot(x='cst',y='grp',hue='time window',data=q,ax=ax)
    return(fig)

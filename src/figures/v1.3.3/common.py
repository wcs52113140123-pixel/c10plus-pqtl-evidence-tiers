import sys, os, numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(__file__))
import c10style as S
from c10style import C, MSL, TIER, GRP, Fig, forest, wilson, fmt_p, FS, FS_T
import matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from matplotlib.lines import Line2D
S.setup()
OUT=os.path.join(os.path.dirname(__file__),'out'); os.makedirs(OUT,exist_ok=True)
def box(ax,x,y,w,h,title,body,ec):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=1.2',lw=0.8,ec=ec,fc='white'))
    ax.text(x+1.5,y+h-1.4,title,ha='left',va='top',fontsize=FS,fontweight='bold',color=ec)
    if body: ax.text(x+1.5,y+h-5.0,body,ha='left',va='top',fontsize=FS,color=C['text'],linespacing=1.25)
def woolf(a,b,c,d):
    o=(a*d)/(b*c); se=np.sqrt(1/a+1/b+1/c+1/d); return o,np.exp(np.log(o)-1.96*se),np.exp(np.log(o)+1.96*se)
GROUPS=['Neither','PAV only','Anchor only','Anchor + PAV']
def grp_bars(ax,df,ylab='Replicated (%)',ymax=None,labels=True,xt=('Neither','PAV\nonly','Anchor\nonly','Anchor\n+ PAV')):
    x=np.arange(4)
    for i,g in enumerate(GROUPS):
        r=df[df.group==g].iloc[0]; k,n=r['k'],r['n']; p=k/n; lo,hi=wilson(k,n)
        ax.bar(i,100*p,0.68,color=GRP[g],ec='none',zorder=2)
        ax.errorbar(i,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.5,capthick=0.7,zorder=3)
        if labels: ax.text(i,100*hi+1.2,f'{100*p:.1f}',ha='center',va='bottom',fontsize=FS)
    ax.set_xticks(x); ax.set_xticklabels(xt); ax.tick_params(axis='x',length=0); ax.set_ylabel(ylab)
    if ymax: ax.set_ylim(0,ymax); ax.set_yticks([t for t in range(0,int(ymax)+1,10) if t<=ymax-3])
    ax.set_xlim(-0.6,3.6)
def title_fit(F):
    """returns titles whose right edge exceeds their panel region by >0.5 mm"""
    r=F.fig.canvas.get_renderer(); bad=[]
    for p in F.panels:
        if p['T'] is None: continue
        bb=p['T'].get_window_extent(r); right=bb.x1/F.fig.dpi*25.4; x,y,w,h=p['region']
        if right>x+w+0.5: bad.append((p['letter'],round(right-(x+w),1)))
    return bad

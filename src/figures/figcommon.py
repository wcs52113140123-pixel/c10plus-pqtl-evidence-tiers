# Common palette, helpers and QA for C10+ manuscript figures (Genome Medicine v0)
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
# ---- threaded palette ----
C_NONE   = '#AFC3DD'   # no anchor, no PAV
C_ANCH   = '#2F5D9E'   # anchored, no PAV  (= Tier 1)
C_PAVNA  = '#E07B39'   # PAV, no anchor    (= Tier 4, caution class)
C_ANCHPAV= '#8E4A1C'   # anchored + PAV
C_PAV0   = '#5C7FB0'   # PAV-absent (binary contrasts)
C_PAV1   = C_PAVNA     # PAV-present
TIER_C   = {1:C_ANCH, 2:'#8FB3DE', 3:'#BDBDBD', 4:C_PAVNA}
TIER_L   = {1:'T1 anchored, no PAV',2:'T2 anchored + PAV, or isoform-only',3:'T3 no anchor, no PAV',4:'T4 no anchor + PAV'}
L_AFF    = '#6A5A8C'   # affinity replication layers (SomaScan)
L_MS     = '#2A8C82'   # mass spectrometry layers
L_XA     = '#9A7B2F'   # cross-ancestry (same platform)
L_DISC   = '#444444'   # discovery / neutral
GREY     = '#8A8A8A'
PH_FACE  = '#F2F2F2'   # placeholder face
def wilson(k,n,z=1.96):
    k=np.asarray(k,float); n=np.asarray(n,float)
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return c-h, c+h
def fisher_ci(r,n):
    z=np.arctanh(r); se=1/np.sqrt(n-3); return np.tanh(z-1.96*se), np.tanh(z+1.96*se)
def woolf(a,b,c,d):
    orr=(a*d)/(b*c); se=np.sqrt(1/a+1/b+1/c+1/d); return orr, np.exp(np.log(orr)-1.96*se), np.exp(np.log(orr)+1.96*se)
def fmt_p(p):
    if p is None or (isinstance(p,float) and np.isnan(p)): return ''
    if p<1e-3:
        e=int(np.floor(np.log10(p))); m=p/10**e
        return f'{m:.1f}×10$^{{{e}}}$'
    return f'{p:.2g}' if p<0.1 else f'{p:.2f}'
def bbox_check(fig):
    r=fig.canvas.get_renderer()
    texts=[(t,t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    ticks=set()
    for ax in fig.axes:
        ticks|=set(ax.get_xticklabels(which='both')+ax.get_yticklabels(which='both'))
    ov=[]
    for i,(a,ba) in enumerate(texts):
        for b,bb in texts[i+1:]:
            if ba.overlaps(bb): ov.append((a.get_text()[:30],b.get_text()[:30]))
    fb=fig.bbox
    out=[t.get_text()[:30] for t,b in texts if b.x0<fb.x0-1 or b.x1>fb.x1+1 or b.y0<fb.y0-1 or b.y1>fb.y1+1]
    return ov,out
def save(fig,name,outdir='work/figs'):
    fig.savefig(f'{outdir}/{name}.png',dpi=300)
    fig.savefig(f'{outdir}/{name}.pdf')
    ov,out=bbox_check(fig)
    print(name,'overlaps:',ov[:12],'outside:',out[:8])

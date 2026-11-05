"""C10+ v1.3.3 figure style: one palette, one font, mm layout with aligned panel letters/titles, verification."""
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, json
from matplotlib.patches import FancyBboxPatch, Rectangle, Patch
from matplotlib.lines import Line2D
MM=1/25.4
# ---------- fixed semantic palette (used identically in every figure) ----------
C=dict(
 anchor='#2F5D9E',      # eQTL anchor (and Tier 1)
 anchor_lt='#9BB5DB',   # isoform-only / anchored+PAV tier (Tier 2)
 pav='#E07B39',         # PAV (and Tier 4)
 anchorpav='#7A3B12',   # anchor + PAV
 neither='#C9C9C9',     # neither feature (and Tier 3)
 pavneg='#8C8C8C',      # PAV-negative comparator in binary contrasts
 absent='#D9D9D9',      # feature absent (MS feature bars)
 conc='#6A5A8C',        # platform concordance / affinity layers (SomaScan)
 ms='#2A8C82',          # mass spectrometry family (generic)
 xa='#9A7B2F',          # non-European ancestries
 mr='#4D4D4D',          # published MR
 ctrl='#7F7F7F',        # wrong-gene negative control
 full='#1A1A1A', b1='#8C8C8C', b0='#C9C9C9',   # grading model vs baselines
 text='#1A1A1A', sub='#555555', ref='#9A9A9A')
MSL={'MSPA':'#0F4F48','Seer':'#2A8C82','Generation Scotland':'#5FB3A9','HOLBAEK':'#9ED3CB'}  # MS layer shades (dark->light)
TIER={1:C['anchor'],2:C['anchor_lt'],3:C['neither'],4:C['pav']}
GRP={'Neither':C['neither'],'PAV only':C['pav'],'Anchor only':C['anchor'],'Anchor + PAV':C['anchorpav']}
FS_T=8; FS=7; FS_L=9   # titles; all other text; panel letters
def setup():
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','Helvetica','DejaVu Sans'],
      'font.size':FS,'axes.titlesize':FS_T,'axes.labelsize':FS,'xtick.labelsize':FS,'ytick.labelsize':FS,'legend.fontsize':FS,
      'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,'ytick.major.size':2.5,
      'xtick.direction':'out','ytick.direction':'out','axes.spines.top':False,'axes.spines.right':False,
      'legend.frameon':False,'legend.handlelength':1.2,'legend.handletextpad':0.4,'legend.borderaxespad':0.2,'legend.labelspacing':0.3,
      'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300,'figure.dpi':150,'lines.linewidth':1.0,'axes.titlepad':0,
      'mathtext.default':'regular','axes.unicode_minus':True})
def wilson(k,n,z=1.96):
    k=np.asarray(k,float); n=np.asarray(n,float); p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return np.minimum(np.clip(c-h,0,1),p),np.maximum(np.clip(c+h,0,1),p)
SUPMAP=str.maketrans('⁻⁰¹²³⁴⁵⁶⁷⁸⁹','-0123456789')
def sup(t):
    import re
    return re.sub(r'10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)',lambda m:'10$^{'+m.group(1).translate(SUPMAP)+'}$',t)
def fmt_p(p,prefix='P'):
    if p is None or not np.isfinite(p): return ''
    if p<1e-3:
        e=int(np.floor(np.log10(p))); m=p/10**e
        if round(m,1)>=10: m/=10; e+=1
        return f'{prefix} = {m:.1f}×10$^{{{e}}}$'
    if p<0.01: return f'{prefix} = {p:.2g}'
    return f'{prefix} = {p:.2g}' if p<0.1 else f'{prefix} = {p:.2f}'
class Fig:
    """Figure in mm. Panels are regions (x,y_top,w,h) in mm from the top-left; the letter sits at the region's top-left,
    the title is left-aligned 3.2 mm to its right on the same top line, and the axes fill the region minus margins."""
    def __init__(s,w_mm,h_mm,name):
        s.w,s.h,s.name=w_mm,h_mm,name; s.fig=plt.figure(figsize=(w_mm*MM,h_mm*MM)); s.panels=[]; s.extra_ok=set()
    def ax_mm(s,x,y,w,h,**kw):  # x,y top-left in mm
        return s.fig.add_axes([x/s.w,1-(y+h)/s.h,w/s.w,h/s.h],**kw)
    def panel(s,letter,title,x,y,w,h,ml=10,mr=1,mt=6,mb=8,row=None,**kw):
        ax=s.ax_mm(x+ml,y+mt,w-ml-mr,h-mt-mb,**kw)
        L=s.fig.text(x/s.w,1-y/s.h,letter,fontsize=FS_L,fontweight='bold',va='top',ha='left') if letter else None
        T=s.fig.text((x+3.2)/s.w,1-(y+0.35)/s.h,title,fontsize=FS_T,va='top',ha='left',color=C['text']) if title else None
        s.panels.append(dict(letter=letter,title=title,L=L,T=T,ax=ax,row=row if row is not None else y,region=(x,y,w,h)))
        return ax
    def verify(s):
        f=s.fig; f.canvas.draw(); r=f.canvas.get_renderer(); H=f.bbox.height; W=f.bbox.width
        hidden=set()
        for ax in f.axes:
            for axis,lim in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
                lo,hi=min(lim),max(lim)
                for tk in axis.get_major_ticks()+axis.get_minor_ticks():
                    if not (lo-1e-9<=tk.get_loc()<=hi+1e-9): hidden.add(tk.label1); hidden.add(tk.label2)
        texts=[(t,t.get_window_extent(r)) for t in f.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible() and t not in hidden]
        ticks=set()
        for ax in f.axes: ticks|=set(ax.get_xticklabels(which='both')+ax.get_yticklabels(which='both'))
        ov=[]
        for i,(a,ba) in enumerate(texts):
            for b,bb in texts[i+1:]:
                if ba.overlaps(bb) and (id(a),id(b)) not in s.extra_ok: ov.append((a.get_text()[:28],b.get_text()[:28]))
        spines=[(sp,sp.get_window_extent(r)) for ax in f.axes if ax.axison for sp in ax.spines.values() if sp.get_visible()]
        for t,bt in texts:
            for sp,bs in spines:
                if bt.overlaps(bs) and t not in ticks and not getattr(t,'_ok_on_spine',False): ov.append((t.get_text()[:28],'spine'))
        out=[t.get_text()[:28] for t,b in texts if b.x0<-0.5 or b.x1>W+0.5 or b.y0<-0.5 or b.y1>H+0.5]
        px2mm=lambda v: v/f.dpi*25.4
        rows={}
        for p in s.panels:
            rec={}
            if p['L'] is not None: bb=p['L'].get_window_extent(r); rec['letter_top_mm']=round(px2mm(H-bb.y1),2); rec['letter_left_mm']=round(px2mm(bb.x0),2)
            if p['T'] is not None: bb=p['T'].get_window_extent(r); rec['title_top_mm']=round(px2mm(H-bb.y1),2); rec['title_left_mm']=round(px2mm(bb.x0),2)
            pos=p['ax'].get_position(); rec['axes_top_mm']=round((1-pos.y1)*s.h,2); rec['axes_bottom_mm']=round((1-pos.y0)*s.h,2)
            rec['region_left_mm']=p['region'][0]
            rows.setdefault(p['row'],[]).append((p['letter'],rec))
        align=[]
        for rk,items in rows.items():
            for key in ['letter_top_mm','title_top_mm']:
                vals=[it[1][key] for it in items if key in it[1]]
                if len(vals)>1: align.append((rk,key,max(vals)-min(vals)))
        axspread=[]
        for rk,items in rows.items():
            vals=[it[1]['axes_top_mm'] for it in items]
            if len(vals)>1: axspread.append((rk,round(max(vals)-min(vals),2)))
        lead=[ (L,round(rec.get('letter_left_mm',0)-rec['region_left_mm'],2)) for items in rows.values() for L,rec in items if L]
        return dict(figure=s.name,overlaps=ov,outside=out,rows={str(k):v for k,v in rows.items()},max_row_misalignment_mm=max([a[2] for a in align],default=0.0),axes_top_spread_mm=axspread,letter_offset_from_region_mm=lead)
    def save(s,base):
        s.fig.savefig(base+'.pdf'); s.fig.savefig(base+'.png',dpi=300)
        try:
            from PIL import Image; Image.open(base+'.png').convert('RGB').save(base+'.tif',compression='tiff_lzw',dpi=(300,300))
        except Exception as e: print('tif',e)
def forest(ax,rows,xlog=True,ref=1.0,ms=3.2,lab_x=None,fmt='{:.2f} ({:.2f}–{:.2f})',show_vals=True,val_x=None):
    """rows: list of dict(label,est,lo,hi,color,marker='o',open=False,header=False,thr=None)."""
    n=len(rows); ys=np.arange(n)[::-1]
    for y,r in zip(ys,rows):
        if r.get('header'):
            ax.text(-0.01,y,r['label'],transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=FS,fontweight='bold',color=r.get('color',C['text'])); continue
        ax.plot([r['lo'],r['hi']],[y,y],color=r['color'],lw=1.1,solid_capstyle='butt',zorder=2)
        fc='white' if r.get('open') else r['color']
        ax.plot(r['est'],y,r.get('marker','o'),ms=r.get('ms',ms),mfc=fc,mec=r['color'],mew=0.9,zorder=3)
        if r.get('thr') is not None: ax.plot(r['thr'],y,marker='|',ms=6,color='#333333',mew=0.8,zorder=1)
        if show_vals and val_x is not None and r.get('est') is not None:
            ax.text(val_x,y,fmt.format(r['est'],r['lo'],r['hi']),transform=ax.get_yaxis_transform(),fontsize=FS,va='center',ha='left')
    ax.axvline(ref,color=C['ref'],lw=0.6,ls=(0,(3,2)),zorder=0)
    if xlog:
        ax.set_xscale('log'); ax.minorticks_off(); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_yticks([y for y,r in zip(ys,rows) if not r.get('header')]); ax.set_yticklabels([r['label'] for r in rows if not r.get('header')])
    ax.set_ylim(-0.7,n-0.3); ax.tick_params(axis='y',length=0)
    ax.spines['left'].set_visible(False)
    return ys

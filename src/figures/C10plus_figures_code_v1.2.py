# SUPERSEDED (release v1.3.3): replaced by src/figures/v1.3.3/ (run_all.py). Kept for the record; the figures it produced were removed from figures/.
# Inputs (place under ./inputs/ relative to the working directory; original store-relative names kept):
#   87fcab6d-.../v944352bd_C10plus_fullscale_results_2026-11-03.xlsx  -> results workbook v3 (Additional file 5 v4 is a superset)
#   d50db924-.../v055bbf03_C10plus_cis_pqtl_robustness_grades_v2.tsv.gz -> resource table (Additional file 3; v2.2 contains all v2 columns)
#   bc9c00a1-.../v83e32977_Fig3abc_roc_source.csv -> ROC source (code/analysis/Fig3abc_roc_source.csv)
# Requires figcommon.py on the import path (sys.path insert 'work/figs' in the setup block).
# C10+ manuscript figure code v1.2 (v1.1 + addendum A2: Fig. 1a/b text and rows, new Fig. 3 with panels d–f) (layout fixes: Fig2c legend; FigS2 ticks/titles/letters; S4 x ticks; S5 title; S6 margin/title). Run setup_v1.py (loads workbook v3, resource v2, ROC source, figure-style helpers), then the v0 function library, then the v1 overrides below.

# ===== figcommon.py =====
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

# ===== setup_v1.py =====
import sys, json, re, os, glob, inspect, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
sys.path.insert(0,'work/figs')
from figcommon import *
from matplotlib.patches import Patch, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
A='inputs/'
WB3=A+'87fcab6d-b91b-4909-bf08-990431447993/v944352bd_C10plus_fullscale_results_2026-11-03.xlsx'
xl=pd.read_excel(WB3,sheet_name=None)
D={}
for k in ['H1_detail','H2_detail','H3_detail','H4_detail','H5_detail']: D[k[:2]]=dict(zip(xl[k]['key'],xl[k]['value']))
D['A3']=dict(zip(xl['A1_H3_ext_detail']['key'],xl['A1_H3_ext_detail']['value']))
D['A3s']=dict(zip(xl['A1_H3_ext_secondary']['key'],xl['A1_H3_ext_secondary']['value']))
D['A4']=dict(zip(xl['A1_H4_ext_detail']['key'],xl['A1_H4_ext_detail']['value']))
def g(h,k,js=False):
    v=D[h][k]; return json.loads(v) if js else v
H3,H4,H5,A3,A3s,A4=D['H3'],D['H4'],D['H5'],D['A3'],D['A3s'],D['A4']
sens=xl['Sensitivity']; cs_=xl['01_Confirmatory_summary'].set_index('hypothesis'); a1=xl['A1_Addendum_summary'].set_index('hypothesis')
R=pd.read_csv(A+'d50db924-d69e-45f3-a8f1-d52ba983e282/v055bbf03_C10plus_cis_pqtl_robustness_grades_v2.tsv.gz',sep='\t',low_memory=False)
roc=pd.read_csv(A+'bc9c00a1-a59e-4386-b3f3-ddf0f356c262/v83e32977_Fig3abc_roc_source.csv')
asum=xl['Application_summary'].copy(); ap=xl['Application_CS'].copy(); ap['absb']=ap.lead_beta.abs()
ap=ap.merge(R[['studyLocusId','ms_ext_primary']],on='studyLocusId',how='left')
lead=ap.sort_values(['symbol','lead_mlog10p','absb'],ascending=[True,False,False]).groupby('symbol').head(1).set_index('symbol')
best=pd.DataFrame({s:{'best_T1_mlog10p':(d[d.tier==1].lead_mlog10p.max() if (d.tier==1).any() else np.nan),
    'best_T1_variant':(d[d.tier==1].sort_values('lead_mlog10p').lead_variant_id.iloc[-1] if (d.tier==1).any() else '')} for s,d in ap.groupby('symbol')}).T
orr,lo,hi=woolf(2,8,42,39)
ct=pd.crosstab(R.platform_concordance,R.tier)
SRC={}
apply_figure_style(sizes=(8,7,6))

# ===== v0 function library =====
def box(ax,x,y,w,h,title,body,ec,fc_='white',tsize=7,bsize=6,lw=1.0,ls='-'):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.004,rounding_size=0.012',lw=lw,ec=ec,fc=fc_,ls=ls,transform=ax.transAxes)
    ax.add_patch(p)
    ax.text(x+0.008,y+h-0.018,title,ha='left',va='top',fontsize=tsize,fontweight='bold',color=ec,transform=ax.transAxes)
    ax.text(x+0.008,y+h-0.075,body,ha='left',va='top',fontsize=bsize,color='#222222',transform=ax.transAxes,linespacing=1.25)


def arrow(ax,x0,y0,x1,y1,c='#666666'):
    ax.annotate('',xy=(x1,y1),xytext=(x0,y0),xycoords='axes fraction',arrowprops=dict(arrowstyle='-|>',color=c,lw=0.9,shrinkA=0,shrinkB=0))


def fig1():
    fig=plt.figure(figsize=(7.2,5.6))
    ax=fig.add_axes([0.01,0.50,0.98,0.47]); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
    for x,t in [(0.0,'Discovery'),(0.175,'Evidence per credible set'),(0.425,'Evaluation layers (pre-registered)'),(0.785,'Resource')]:
        ax.text(x+0.004,1.0,t,ha='left',va='top',fontsize=6,color='#555555',style='italic',transform=ax.transAxes)
    H=0.16; gap=0.02; tops=[0.92-H-i*(H+gap) for i in range(5)]
    box(ax,0.0,tops[4],0.16,0.92-tops[4],'UKB-PPP Olink','European ancestry\n14,472 cis credible\nsets (CS)\n1,879 genes\n\nOpen Targets 26.09\nSuSiE-inf, 95% sets',L_DISC)
    fh=(0.92-tops[4]-2*gap)/3
    fy=[0.92-fh, 0.92-2*fh-gap, tops[4]]
    box(ax,0.175,fy[0],0.23,fh,'Transcript anchor','coloc PP.H4 ≥ 0.8 with\ncognate eQTL (GTEx v10;\neQTL Catalogue)',C_ANCH)
    box(ax,0.175,fy[1],0.23,fh,'Protein-altering variant','cognate missense/LoF/splice\nvariant in 95% CS (VEP 116)',C_PAVNA)
    box(ax,0.175,fy[2],0.23,fh,'Negative control','coloc with non-cognate\n(wrong-gene) GTEx eQTL',GREY)
    ev=[('SomaScan v3 (Sun 2018)','424 genes · 4,153 CS · H1 (internal)',L_AFF,'white','-'),
        ('SomaScan v4 (deCODE)','888 new genes · 6,711 CS · H2',L_AFF,'white','-'),
        ('Mass spectrometry: MSPA','Suhre 2025 · 169 evaluable CS · H3, H4',L_MS,'white','-'),
        ('Mass spectrometry: Seer full','addendum A1 · n CS (superseded v0 label)',L_MS,PH_FACE,'--'),
        ('Olink, non-EUR ancestries','UKB-PPP AFR/CSA/EAS/MID/AMR\n· 8,065 CS · H5',L_XA,'white','-')]
    for t,(ti,bo,c,f,ls) in zip(tops,ev): box(ax,0.425,t,0.34,H,ti,bo,c,fc_=f,ls=ls)
    box(ax,0.785,tops[4],0.21,0.92-tops[4],'PAV-aware grades','Tier 1–4 for all\n14,472 CS\n\nT1 anchored, no PAV\nT2 anchored + PAV,\n     or isoform-only\nT3 no anchor, no PAV\nT4 no anchor + PAV\n\nCaution flag:\nconcordant + unanchored\n+ PAV (339 CS)',C_PAVNA)
    mid=(0.92+tops[4])/2
    arrow(ax,0.16,mid,0.175,mid); arrow(ax,0.405,mid,0.425,mid); arrow(ax,0.765,mid,0.785,mid)
    # layer-type key
    for i,(c,l) in enumerate([(L_AFF,'affinity, other platform'),(L_MS,'antibody-independent (MS)'),(L_XA,'same platform, other ancestry')]):
        ax.add_patch(Rectangle((0.175+i*0.21,-0.045),0.012,0.03,color=c,transform=ax.transAxes,clip_on=False))
        ax.text(0.175+i*0.21+0.017,-0.03,l,fontsize=6,va='center',transform=ax.transAxes)
    panel_letter(ax,'a',dx=0.005,dy=1.03)
    # ---- b
    axb=fig.add_axes([0.385,0.085,0.29,0.32])
    yy=np.arange(len(rows))[::-1]
    for y,r in zip(yy,rows):
        h,lab,e,l,u,ph,thr,sgn,col,n=r
        axb.plot([l,u],[y,y],color=col,lw=1.4,solid_capstyle='butt')
        axb.plot(e,y,'o',ms=5.5,color=col,mec='white',mew=0.6,zorder=3)
        axb.plot(thr,y,marker='|',ms=9,color='#222222',mew=1.0,zorder=2)
    axb.axvline(1,color='#999999',lw=0.8,ls=(0,(3,2)),zorder=0)
    axb.set_xscale('log'); axb.set_xlim(0.035,6.5)
    axb.set_xticks([0.05,0.1,0.25,0.5,1,2,4]); axb.set_xticklabels(['0.05','0.1','0.25','0.5','1','2','4']); axb.minorticks_off()
    axb.set_yticks(yy); axb.set_yticklabels([f'{r[0]}  {r[1]}' for r in rows])
    axb.set_ylim(-0.7,len(rows)-0.4)
    axb.set_xlabel('Odds ratio (95% CI), log scale')
    axb.text(0.04,-0.55,'|  pre-specified threshold',fontsize=6,va='center')
    fig.text(0.02,0.445,'Anchoring predicts affinity and cross-ancestry replication but not MS support',fontsize=8,ha='left')
    tx=fig.add_axes([0.69,0.085,0.30,0.32]); tx.set_axis_off(); tx.set_ylim(axb.get_ylim()); tx.set_xlim(0,1)
    for x,t in [(0.0,'OR (95% CI)'),(0.40,'Holm P'),(0.64,'Verdict')]: tx.text(x,len(rows)-0.5,t,fontsize=6,va='bottom',color='#555555')
    vs={'H2a':'supported','H2b':'supported','H2a-ARIC':'supported','H2b-ARIC':'supported','H2a-Fenland':'supported','H2b-Fenland':'supported','H5':'supported','H3i':'not supported','H3ii':'partial'}
    for y,r in zip(yy,rows):
        tx.text(0.0,y,f'{r[2]:.2f} ({r[3]:.2f}–{r[4]:.2f})',fontsize=6,va='center')
        tx.text(0.40,y,fmt_p(r[5]),fontsize=6,va='center')
        tx.text(0.64,y,vs[r[0]],fontsize=6,va='center')
    fig.text(0.02,0.445+0.03,'b',fontsize=10,fontweight='bold',va='bottom')
    fig.canvas.draw(); save(fig,'Fig1_design_and_confirmatory_endpoints'); return fig


def fourgroup(h):
    fg=pd.DataFrame(g(h,'four_group',True)); fg['k']=(fg['mean']*fg['size']).round().astype(int)
    lab={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}
    col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    fg['group']=[lab[(a,p)] for a,p in zip(fg.A_gtex,fg.PAV)]; fg['color']=[col[(a,p)] for a,p in zip(fg.A_gtex,fg.PAV)]
    fg['lo'],fg['hi']=wilson(fg.k,fg['size']); return fg


def fig2():
    fig=plt.figure(figsize=(7.2,5.6))
    axa=fig.add_axes([0.085,0.60,0.25,0.33]); axb_=fig.add_axes([0.365,0.60,0.25,0.33])
    for ax,h,tit,key in zip([axa,axb_],['H1','H2'],['Sun 2018 SomaScan v3','deCODE SomaScan v4 (H2)'],['Fig2a','Fig2b']):
        fg=fourgroup(h); x=np.arange(4)
        ax.bar(x,fg['mean'],color=fg.color,width=0.72)
        ax.errorbar(x,fg['mean'],yerr=[fg['mean']-fg.lo,fg.hi-fg['mean']],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
        for xi,(hh,n) in enumerate(zip(fg.hi,fg['size'])): ax.text(xi,hh+0.012,f'n={n:,}',ha='center',va='bottom',fontsize=6)
        ax.set_xticks(x); ax.set_xticklabels(['Neither','PAV\nonly','Anchor\nonly','Anchor\n+PAV'])
        ax.set_ylim(0,0.6); ax.set_title(tit,loc='left',fontsize=8)
        SRC[key]=fg[['A_gtex','PAV','group','k','size','mean','lo','hi']].rename(columns={'k':'replicated','size':'n','mean':'rate','lo':'wilson_lo','hi':'wilson_hi','A_gtex':'GTEx_anchor'})
    axa.set_ylabel('Cross-affinity replication rate'); axb_.set_yticklabels([])
    panel_letter(axa,'a',dx=-0.30); panel_letter(axb_,'b',dx=-0.10)
    ax=fig.add_axes([0.80,0.60,0.185,0.33]); covs=['mlp','absb','maf','logn95']; nm={'mlp':'−log10 P (capped)','absb':'|lead β|','maf':'MAF','logn95':'log 95% CS size'}
    b=[D['H2'][f'pav_smd_before.{c}'] for c in covs]; a=[D['H2'][f'pav_smd_after.{c}'] for c in covs]; y=np.arange(4)[::-1]
    for yi,bb,aa in zip(y,b,a): ax.plot([bb,aa],[yi,yi],color='#CCCCCC',lw=0.8,zorder=1)
    ax.scatter(b,y,s=20,facecolor='white',edgecolor=C_PAVNA,lw=1.0,label='before',zorder=2)
    ax.scatter(a,y,s=20,color=C_PAVNA,label='after 1:3 matching',zorder=3)
    ax.axvspan(-0.1,0.1,color='#EEEEEE',zorder=0); ax.axvline(0,color='#999999',lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([nm[c] for c in covs]); ax.set_xlim(-0.4,0.45); ax.set_xticks([-0.3,0,0.3])
    ax.set_ylim(-0.6,4.9); ax.set_xlabel('SMD (PAV+ − PAV−)'); ax.set_title('Power covariates balanced',loc='left',fontsize=8,x=-0.55)
    ax.legend(loc='upper center',fontsize=6,handletextpad=0.2,borderaxespad=0.3,ncol=2,columnspacing=0.8,frameon=False)
    panel_letter(ax,'c',dx=-0.70)
    SRC['Fig2c']=pd.DataFrame({'covariate':covs,'label':[nm[c] for c in covs],'SMD_before':b,'SMD_after':a,'layer':'deCODE v4 (H2b)','n_PAV_cs_total':D['H2']['pav_n_treated_total'],'n_matched_sets':D['H2']['pav_matched_n_sets']})
    # d
    ax=fig.add_axes([0.085,0.09,0.20,0.36]); rows_=[]
    for i,(h,lab) in enumerate([('H1','Sun 2018'),('H2','deCODE')]):
        n=D[h]['pav_matched_n']; n1=D[h]['pav_matched_n_sets']; n0=n-n1
        r1=D[h]['pav_matched_rates.PAV1']; r0=D[h]['pav_matched_rates.PAV0']; k1=round(r1*n1); k0=round(r0*n0)
        for j,(r,k,nn,c,gl) in enumerate([(r0,k0,n0,C_PAV0,'matched PAV−'),(r1,k1,n1,C_PAV1,'PAV+')]):
            lo_,hi_=wilson(k,nn); xx=i+(j-0.5)*0.38
            ax.bar(xx,r,width=0.36,color=c); ax.errorbar(xx,r,yerr=[[r-lo_],[hi_-r]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            rows_.append(dict(layer=lab,group=gl,replicated=k,n=nn,rate=r,wilson_lo=lo_,wilson_hi=hi_))
        ax.text(i,0.40,f'OR {D[h]["pav_matched.PAV.OR"]:.2f}\nP = {fmt_p(D[h]["pav_matched.PAV.P"])}',ha='center',fontsize=6)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Sun 2018\n(internal)','deCODE\n(H2b)']); ax.set_ylim(0,0.6); ax.set_xlim(-0.6,1.6)
    ax.set_ylabel('Replication rate'); ax.set_title('PAV excess survives\npower matching',loc='left',fontsize=8)
    ax.legend(handles=[Patch(color=C_PAV0,label='matched PAV−'),Patch(color=C_PAV1,label='PAV+')],loc='upper left',fontsize=6,borderaxespad=0.0,handlelength=1,ncol=1)
    panel_letter(ax,'d',dx=-0.40)
    SRC['Fig2d']=pd.DataFrame(rows_); SRC['Fig2d']['note']='PAV+ n = matched sets; PAV− n = matched sample minus PAV+; counts back-calculated from rates; OR from matched-set clustered logistic'
    # e
    ax=fig.add_axes([0.60,0.09,0.24,0.36]); fr=[]
    for h,lab in [('H1','Sun 2018 (internal)'),('H2','deCODE (H2a, confirmatory)')]:
        for key,el,mk,fcol in [('main.A_gtex','GTEx eQTL anchor','o',L_AFF),('nongtex.A_nongtex','eQTL Catalogue anchor','s',L_AFF),('power_matched_subset.A_gtex','GTEx anchor, strong-signal subset','D',L_AFF),('negctrl_wrong.W_gtex','wrong-gene eQTL, alone','o',GREY),('joint.W_gtex','wrong-gene eQTL, joint model','o','white')]:
            fr.append(dict(layer=lab,estimate=el,OR=D[h][f'{key}.OR'],lo=D[h][f'{key}.lo'],hi=D[h][f'{key}.hi'],P=D[h][f'{key}.P'],n=D[h][f'{key}.n'],marker=mk,face=fcol))
    fr=pd.DataFrame(fr); ys=[]
    for i,r in fr.iterrows():
        y=-(i+1+(1 if i>=5 else 0)); ec=GREY if r.face in (GREY,'white') else r.face
        ax.plot([r.lo,r.hi],[y,y],color=ec,lw=1.2); ax.plot(r.OR,y,r.marker,ms=4.5,mfc=r.face,mec=ec,mew=0.9)
        ax.text(1.03,y,f'{r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})',va='center',fontsize=6,transform=ax.get_yaxis_transform())
        ys.append(y)
    ax.set_yticks(ys); ax.set_yticklabels(fr.estimate)
    for yh,t in [(0,'Sun 2018 (internal)'),(-6,'deCODE (H2a, confirmatory)')]:
        ax.text(-0.03,yh,t,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6,fontweight='bold',color=L_AFF)
    ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2)))
    ax.set_xscale('log'); ax.set_xlim(0.8,5.5); ax.set_xticks([1,2,3,5]); ax.set_xticklabels(['1','2','3','5']); ax.minorticks_off()
    ax.set_xlabel('OR for cross-affinity replication (95% CI)'); ax.set_ylim(-11.6,0.6)
    ax.set_title('Cognate, not wrong-gene, eQTL coloc predicts replication',loc='left',fontsize=8,x=-1.0)
    panel_letter(ax,'e',dx=-1.12)
    SRC['Fig2e']=fr.drop(columns=['marker','face'])
    fig.canvas.draw(); save(fig,'Fig2_cross_affinity_replication_PAV_matching'); return fig


def fig3():
    fig=plt.figure(figsize=(7.2,2.8))
    axs=[fig.add_axes([0.075,0.17,0.25,0.68]),fig.add_axes([0.40,0.17,0.25,0.68]),fig.add_axes([0.73,0.17,0.255,0.68])]
    rowsA=[];rowsB=[]
    for i,(h,lab) in enumerate([('H1','Sun 2018'),('H2','deCODE')]):
        zt=json.loads(D[h]['effect_corr.fisher_z_PAV1_lt_PAV0']) if isinstance(D[h]['effect_corr.fisher_z_PAV1_lt_PAV0'],str) else D[h]['effect_corr.fisher_z_PAV1_lt_PAV0']
        for j,(pv,c,gl) in enumerate([(0,C_PAV0,'PAV−'),(1,C_PAV1,'PAV+')]):
            n=D[h][f'effect_corr.PAV{pv}.n']; r=D[h][f'effect_corr.PAV{pv}.r']; rs=D[h][f'effect_corr.PAV{pv}.r_samedir']; ns=D[h][f'effect_corr.PAV{pv}.n_samedir']
            lo_,hi_=fisher_ci(r,n); xx=i+(j-0.5)*0.36
            axs[0].plot([xx,xx],[lo_,hi_],color=c,lw=1.4); axs[0].plot(xx,r,'o',ms=6,color=c,mec='white',mew=0.6,zorder=3)
            axs[0].plot(xx+0.10,rs,'o',ms=4.5,mfc='white',mec=c,mew=1.0,zorder=3)
            rowsA.append(dict(layer=lab,PAV=gl,n_pairs=n,pearson_r=r,fisherCI_lo=lo_,fisherCI_hi=hi_,r_same_direction_only=rs,n_same_direction=ns))
            no=D[h][f'effect_corr.PAV{pv}.n_opp']; fr_=no/n; lo2,hi2=wilson(no,n)
            axs[1].bar(xx,fr_,width=0.34,color=c); axs[1].errorbar(xx,fr_,yerr=[[fr_-lo2],[hi2-fr_]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            axs[1].text(xx,hi2+0.006,f'{no}/{n}',ha='center',va='bottom',fontsize=6)
            rowsB.append(dict(layer=lab,PAV=gl,n_opposite=no,n_pairs=n,frac=fr_,wilson_lo=lo2,wilson_hi=hi2))
        axs[0].text(i,0.97,f'P = {fmt_p(zt[1])}',ha='center',fontsize=6)
        po=D[h]['effect_corr.opp_fisher_PAV1_gt.P_one_sided']; axs[1].text(i,0.28,f'P = {fmt_p(po)}',ha='center',fontsize=6)
    axs[0].set_ylim(0.0,1.02); axs[0].set_ylabel('Platform-to-platform β correlation (r)')
    axs[0].set_title('PAV lowers β agreement',loc='left',fontsize=8)
    from matplotlib.lines import Line2D
    axs[0].legend(handles=[Line2D([],[],marker='o',ls='',color=C_PAV0,label='PAV−'),Line2D([],[],marker='o',ls='',color=C_PAV1,label='PAV+'),Line2D([],[],marker='o',ls='',mfc='white',mec='#555555',label='same-direction pairs only')],loc='lower left',fontsize=6,borderaxespad=0.1,handletextpad=0.1)
    axs[1].set_ylim(0,0.31); axs[1].set_ylabel('Fraction with opposite direction')
    axs[1].set_title('…and raises sign reversals',loc='left',fontsize=8)
    for a in axs[:2]: a.set_xticks([0,1]); a.set_xticklabels(['Sun 2018','deCODE']); a.set_xlim(-0.6,1.6)
    # c by anchor x PAV
    lab={(0,0):'Neither',(0,1):'PAV\nonly',(1,0):'Anchor\nonly',(1,1):'Anchor\n+PAV'}; col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    rowsC=[]
    for i,(h,ln) in enumerate([('H1','Sun 2018'),('H2','deCODE')]):
        bb=pd.DataFrame(json.loads(D[h]['effect_corr.by_anchor_PAV']))
        for j,r in bb.iterrows():
            k=int(round(r.opp_frac*r.n)); lo2,hi2=wilson(k,r.n); xx=i*5+j
            axs[2].bar(xx,r.opp_frac,width=0.75,color=col[(r.A_eqtl,r.PAV)]); axs[2].errorbar(xx,r.opp_frac,yerr=[[r.opp_frac-lo2],[hi2-r.opp_frac]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
            axs[2].text(xx,hi2+0.008,f'{k}/{int(r.n)}',ha='center',va='bottom',fontsize=5.5 if False else 6,rotation=90)
            rowsC.append(dict(layer=ln,A_eqtl=r.A_eqtl,PAV=r.PAV,n_opposite=k,n_pairs=int(r.n),frac=r.opp_frac,wilson_lo=lo2,wilson_hi=hi2))
    axs[2].set_xticks([1.5,6.5]); axs[2].set_xticklabels(['Sun 2018','deCODE']); axs[2].set_ylim(0,0.48)
    axs[2].set_ylabel('Fraction with opposite direction'); axs[2].set_title('Sign reversals by anchor and PAV status',loc='left',fontsize=8)
    axs[2].legend(handles=[Patch(color=col[k],label=lab[k].replace('\n',' ')) for k in [(0,0),(0,1),(1,0),(1,1)]],loc='upper right',fontsize=6,ncol=2,borderaxespad=0.1,handlelength=1,columnspacing=0.6)
    for a,l in zip(axs,'abc'): panel_letter(a,l,dx=-0.30)
    SRC['Fig3a']=pd.DataFrame(rowsA); SRC['Fig3b']=pd.DataFrame(rowsB); SRC['Fig3c']=pd.DataFrame(rowsC)
    SRC['Fig3a']['note']='pairs = replicated CS with harmonisable lead/shared variant; Fisher z one-sided P (PAV+ < PAV−): Sun 2.4e-06, deCODE 0.0055'
    fig.canvas.draw(); save(fig,'Fig3_effect_agreement_sign_reversal'); return fig


def fig4():  # v1.2 Seer-only panel; the published Fig. 4 is figures/Fig4_mass_spec_four_layers.png (all four MS layers)
    fig=plt.figure(figsize=(7.2,5.9))
    # a
    ax=fig.add_axes([0.085,0.62,0.36,0.31]); rows_=[]
    feats=[('concordant_fisher','Platform\nconcordant',L_AFF),('H3i_fisher','eQTL\nanchored',C_ANCH),('PAV_fisher','PAV in\n95% CS',C_PAVNA)]
    for i,(k,lab,c) in enumerate(feats):
        t=json.loads(H3[f'primary.{k}.table']) ; # [[sup_g1, notsup_g1],[sup_g0,notsup_g0]]
        for j,(sup,nsup,face,gl) in enumerate([(t[0][0],t[0][1],c,'yes'),(t[1][0],t[1][1],'#C8C8C8','no')]):
            n=sup+nsup; r=sup/n; lo_,hi_=wilson(sup,n); xx=i+(j-0.5)*0.38
            ax.bar(xx,r,width=0.36,color=face); ax.errorbar(xx,r,yerr=[[r-lo_],[hi_-r]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            ax.text(xx,0.015,f'{sup}/{n}',ha='center',va='bottom',fontsize=6,color='white' if gl=='yes' else '#222222',rotation=90)
            ax.text(xx,-0.035,gl,ha='center',va='top',fontsize=6)
            rows_.append(dict(feature=lab.replace('\n',' '),status=gl,supported=sup,n=n,rate=r,wilson_lo=lo_,wilson_hi=hi_))
        ax.text(i,0.80,f'OR {H3[f"primary.{k}.OR"]:.2f}\nP = {fmt_p(H3[f"primary.{k}.P"])}',ha='center',fontsize=6)
    ax.set_xticks(range(3)); ax.set_xticklabels([f[1] for f in feats]); ax.tick_params(axis='x',pad=10)
    ax.set_ylim(0,0.92); ax.set_ylabel('MS-supported fraction')
    ax.set_title('MS-supported fraction by concordance, anchoring and PAV status (Seer-study labels; superseded panel)',loc='left',fontsize=8)
    panel_letter(ax,'a',dx=-0.20)
    SRC['Fig4a']=pd.DataFrame(rows_); SRC['Fig4a']['test']='two-sided Fisher exact; n evaluable = 169 CS (concordance: 144 assayed on SomaScan)'
    # b
    ax=fig.add_axes([0.56,0.62,0.42,0.31]); rows_=[]
    gc={'CA':C_ANCH,'CN0':C_NONE,'CNP':C_PAVNA}; gl={'CA':'Anchored','CN0':'No anchor\nno PAV','CNP':'No anchor\n+ PAV'}
    for i,(lab,nm_) in enumerate([('primary','Bonferroni MS replication\n(primary label)'),('mspa','MSPA peptide score\n(secondary label)')]):
        cg=pd.DataFrame(json.loads(H3[f'{lab}.concordant_groups'])).set_index('grp')
        for j,gk in enumerate(['CA','CN0','CNP']):
            k=int(cg.loc[gk,'sum']); n=int(cg.loc[gk,'size']); r=k/n; lo_,hi_=wilson(k,n); xx=i*4.1+j*1.25
            ax.bar(xx,r,width=0.8,color=gc[gk]); ax.errorbar(xx,r,yerr=[[r-lo_],[hi_-r]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            ax.text(xx,hi_+0.02,f'{k}/{n}',ha='center',va='bottom',fontsize=6)
            rows_.append(dict(label=nm_.replace('\n',' '),group=gk,group_label=gl[gk].replace('\n',' '),supported=k,n=n,rate=r,wilson_lo=lo_,wilson_hi=hi_))
        ax.text(i*4.1+1.25,1.0,('Primary label' if i==0 else 'MSPA peptide label')+f'\nCNP vs CA: OR {H3[f"{lab}.H3ii_CNP_vs_CA.OR"]:.2f}, P = {fmt_p(H3[f"{lab}.H3ii_CNP_vs_CA.P"])}',ha='center',fontsize=6)
    ax.set_xticks([0,1.25,2.5,4.1,5.35,6.6]); ax.set_xticklabels([gl[k] for k in ['CA','CN0','CNP']]*2,fontsize=6)
    ax.set_ylim(0,1.15); ax.set_ylabel('MS-supported fraction')
    ax.set_title('Among platform-concordant signals',loc='left',fontsize=8)
    panel_letter(ax,'b',dx=-0.13)
    SRC['Fig4b']=pd.DataFrame(rows_); SRC['Fig4b']['test']='two-sided Fisher exact CNP vs CA; Holm-adjusted P (primary) = 0.28'
    # c forest
    ax=fig.add_axes([0.30,0.08,0.24,0.42]); fr=[]
    def add(lab,OR,lo,hi,P,n,grp): fr.append(dict(analysis=lab,OR=OR,lo=lo,hi=hi,P=P,n=n,group=grp))
    add('Primary (H3i)',H3['primary.H3i_logit.A_eqtl.OR'],H3['primary.H3i_logit.A_eqtl.lo'],H3['primary.H3i_logit.A_eqtl.hi'],H3['primary.H3i_logit.A_eqtl.P'],169,'label')
    add('Wrong-gene eQTL (neg. control)',H3['primary.H3i_negctrl_wrong.W_gtex.OR'],H3['primary.H3i_negctrl_wrong.W_gtex.lo'],H3['primary.H3i_negctrl_wrong.W_gtex.hi'],H3['primary.H3i_negctrl_wrong.W_gtex.P'],169,'control')
    for k,l in [('mspa','MSPA peptide label'),('primary_ST10_only','UKB-PPP rows only (ST10)'),('primary_excl_pep_pav','Excluding PAV-peptide rows'),('seer','Seer discovery overlap'),('seer_leadcs','Seer, lead CS per gene')]:
        add(l,H3[f'{k}.H3i_logit.A_eqtl.OR'],H3[f'{k}.H3i_logit.A_eqtl.lo'],H3[f'{k}.H3i_logit.A_eqtl.hi'],H3[f'{k}.H3i_logit.A_eqtl.P'],H3[f'{k}.H3i_logit.A_eqtl.n'],'label')
    snm={'h4_0.9':'PP.H4 ≥ 0.9','h4_0.95':'PP.H4 ≥ 0.95','prior_p12_5e-6':'p12 = 5×10⁻⁶','prior_p12_1e-6':'p12 = 1×10⁻⁶','coloc_intersection':'coloc on shared variants','window_lead_TSS_le_1Mb':'lead ≤ 1 Mb of TSS','window_lead_TSS_le_500kb':'lead ≤ 500 kb of TSS','lead_cs_per_gene':'lead CS per gene','exclude_trans_hotspot':'excluding trans hotspots'}
    for sc,l in snm.items():
        r=sens[(sens.scenario==sc)&(sens.estimate=='H3i_Aeqtl')].iloc[0]; add(l,r.OR,r.lo,r.hi,r.P,int(r.n),'sensitivity')
    fr=pd.DataFrame(fr); ys=[]
    for i,r in fr.iterrows():
        y=-(i+ (1 if i>=2 else 0) + (1 if i>=7 else 0)); c=GREY if r.group=='control' else L_MS
        ax.plot([r.lo,r.hi],[y,y],color=c,lw=1.2); ax.plot(r.OR,y,'o',ms=4,color=c,mec='white',mew=0.5,zorder=3)
        ax.text(1.02,y,f'{r.OR:.2f}  n={int(r.n):,}',transform=ax.get_yaxis_transform(),va='center',fontsize=6); ys.append(y)
    ax.set_yticks(ys); ax.set_yticklabels(fr.analysis)
    ax.text(-0.03,-2.0,'Alternative MS labels',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6,fontweight='bold',color=L_MS)
    ax.text(-0.03,-8.0,'Pre-specified sensitivity (primary label)',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6,fontweight='bold',color=L_MS)
    ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2))); ax.axvline(1.5,color='#222222',lw=0.6,ls=':')
    ax.set_xscale('log'); ax.set_xlim(0.2,10); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8']); ax.minorticks_off()
    ax.set_ylim(min(ys)-1.7,0.7); ax.set_xlabel('OR, eQTL anchor → MS support (95% CI)')
    ax.text(0.98,0.01,'dotted line: pre-specified OR 1.5',transform=ax.transAxes,fontsize=6,ha='right',va='bottom')
    ax.set_title('Anchor–MS association is null across labels and scenarios',loc='left',fontsize=8,x=-0.85)
    panel_letter(ax,'c',dx=-1.0)
    SRC['Fig4c']=fr; SRC['Fig4c']['note']='gene-clustered logistic ms ~ A_eqtl + mlp; Seer rows additionally power-adjusted; zero-event cells not interpreted'
    # d placeholder table
    ax=fig.add_axes([0.70,0.08,0.29,0.44]); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.add_patch(Rectangle((0,0),1,0.93,fc=PH_FACE,ec=L_MS,ls='--',lw=0.9))
    ax.text(0.0,0.98,'Seer full extension (addendum A1): pending',fontsize=7,va='bottom',ha='left')
    rr=[('H3i OR (anchor → MS)',f'{H3["primary.H3i_logit.A_eqtl.OR"]:.2f}'),('H3ii OR (CNP vs CA)',f'{H3["primary.H3ii_CNP_vs_CA.OR"]:.2f}'),('H4 ΔAUC (held-out)',f'{D["H4"]["H4_ms.delong_full_vs_B0.dauc"]:.2f}'),('n evaluable CS','169')]
    y=0.88; recs=[]
    for ep,v in rr:
        ax.text(0.04,y,ep,fontsize=6,fontweight='bold',va='top'); y-=0.055
        for src_,val in [('MSPA (pre-registered)',v),('Seer full (A1)','superseded (v0)'),('ext-new subset','superseded (v0)')]:
            ax.text(0.08,y,src_,fontsize=6,va='top'); ax.text(0.62,y,val,fontsize=6,va='top',color='#222222' if not val.startswith('[[') else '#7A7A7A'); y-=0.048
            recs.append(dict(endpoint=ep,source=src_,value=val))
        y-=0.012
    recs.append(dict(endpoint='label agreement kappa',source='MSPA vs Seer full',value='superseded (v0)'))
    panel_letter(ax,'d',dx=-0.03,dy=1.04)
    SRC['Fig4d']=pd.DataFrame(recs)
    fig.canvas.draw(); save(fig,'FigS_v1.2_Seer_MS_layers_superseded'); return fig


def fig5():
    fig=plt.figure(figsize=(7.2,5.6))
    mstyle={'Full grade model':dict(color='#222222',lw=1.4),'Concordance + power':dict(color='#8A8A8A',lw=1.1),'Concordance only':dict(color='#BBBBBB',lw=1.1,ls='--')}
    panels=sorted(roc.panel.unique(),key=lambda s: 0 if 'MS' in s else 1)
    aucs={panels[0]:{'Full grade model':H4['H4_ms.auc_full'],'Concordance + power':H4['H4_ms.auc_B1'],'Concordance only':H4['H4_ms.auc_B0']},
          panels[1]:{'Full grade model':H4['H4_xa.auc_full'],'Concordance + power':H4['H4_xa.auc_B1'],'Concordance only':H4['H4_xa.auc_B0']}}
    tit={panels[0]:'Held-out MS support\n(n = 70)',panels[1]:'Held-out cross-ancestry\n(n = 3,493)'}
    for i,p in enumerate(panels):
        ax=fig.add_axes([0.07+i*0.25,0.63,0.20,0.30])
        for m in ['Concordance only','Concordance + power','Full grade model']:
            d=roc[(roc.panel==p)&(roc.model==m)]; ax.plot(d.fpr,d.tpr,label=f'{m} ({aucs[p][m]:.2f})',**mstyle[m])
        ax.plot([0,1],[0,1],color='#DDDDDD',lw=0.6,zorder=0)
        ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.02); ax.set_aspect('equal')
        ax.set_xlabel('False-positive rate'); 
        if i==0: ax.set_ylabel('True-positive rate')
        ax.set_title(tit[p],loc='left',fontsize=8)
        ax.text(0.98,0.03+0.0,'\n'.join([f'{m.split()[0] if False else m}: {aucs[p][m]:.2f}' for m in ['Full grade model','Concordance + power','Concordance only']]),transform=ax.transAxes,ha='right',va='bottom',fontsize=6) if False else None
        for k_,m in enumerate(['Full grade model','Concordance + power','Concordance only']): ax.text(0.98,0.04+k_*0.085,f'{aucs[p][m]:.2f}',transform=ax.transAxes,ha='right',va='bottom',fontsize=6,color=mstyle[m]['color'] if m!='Concordance only' else '#999999')
        ax.text(0.98,0.04+3*0.085,'AUC',transform=ax.transAxes,ha='right',va='bottom',fontsize=6)
        panel_letter(ax,'ab'[i],dx=-0.30)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([],[],**mstyle[m],label=m) for m in ['Full grade model','Concordance + power','Concordance only']],loc='upper left',bbox_to_anchor=(0.07,0.535),ncol=3,fontsize=6,handlelength=1.8,columnspacing=1.0)
    SRC['Fig5ab']=roc.copy()
    SRC['Fig5ab_auc']=pd.DataFrame([dict(target=p,model=m,AUC=v) for p in aucs for m,v in aucs[p].items()]+[
        dict(target='MS (H4 primary)',model='DeLong one-sided P full vs concordance-only',AUC=H4['H4_ms.delong_full_vs_B0.P_one_sided']),
        dict(target='MS (H4 primary)',model='permutation P (1,000 label permutations)',AUC=H4['H4_ms.perm_P']),
        dict(target='MS (H4 primary)',model='calibration slope (full)',AUC=H4['H4_ms.calib_full.slope']),
        dict(target='transfer MS-trained model → held-out cross-ancestry',model='AUC',AUC=H4['transfer_ms_to_xa.auc1']),
        dict(target='transfer MS-trained model → MSPA peptide label',model='AUC',AUC=H4['transfer_ms_to_ms_mspa.auc1']),
        dict(target='transfer MS-trained model → Seer discovery',model='AUC',AUC=H4['transfer_ms_to_ms_seer.auc1'])])
    # c: tier vs outcomes
    ax=fig.add_axes([0.665,0.63,0.325,0.30]); rows_=[]
    for gi,(key,lab) in enumerate([('ms_primary','MS support\n(MSPA primary)'),('xa_rep','Cross-ancestry\nreplication')]):
        tr=pd.DataFrame(json.loads(H4[f'tiers_vs_truth.{key}.rates']))
        for j,r in tr.iterrows():
            k=int(r['sum']); n=int(r['size']); lo_,hi_=wilson(k,n); xx=gi*5+j
            ax.bar(xx,r['mean'],width=0.8,color=TIER_C[int(r.tier)]); ax.errorbar(xx,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
            ax.text(xx,0.02,f'{n:,}',ha='center',va='bottom',fontsize=6,rotation=90,color='white' if int(r.tier) in (1,4) else '#222222')
            rows_.append(dict(outcome=key,tier=int(r.tier),k=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
        ax.text(gi*5+1.5,0.93,f'trend P = {fmt_p(H4[f"tiers_vs_truth.{key}.trend.P_two_sided"])}',ha='center',fontsize=6)
    ax.set_xticks([0,1,2,3,5,6,7,8]); ax.set_xticklabels(['T1','T2','T3','T4']*2)
    for gi,lab in enumerate(['MS support (MSPA)','Cross-ancestry']): ax.text(gi*5+1.5,-0.13,lab,ha='center',va='top',fontsize=6,transform=ax.get_xaxis_transform())
    ax.set_ylim(0,1.0); ax.set_ylabel('Fraction positive'); ax.set_title('Tiers are not monotonic for either outcome',loc='left',fontsize=8)
    panel_letter(ax,'c',dx=-0.22)
    SRC['Fig5c']=pd.DataFrame(rows_); SRC['Fig5c']['note']='Cochran–Armitage trend P (two-sided); n inside bars'
    # d: H5 forest
    ax=fig.add_axes([0.255,0.07,0.27,0.38]); fr=[]
    def add(l,pfx,var='A_eqtl',grp='main'):
        fr.append(dict(analysis=l,OR=H5[f'{pfx}.main.{var}.OR'] if f'{pfx}.main.{var}.OR' in H5 else H5[f'{pfx}.{var}.OR'],
                       lo=H5.get(f'{pfx}.main.{var}.lo',H5.get(f'{pfx}.{var}.lo')),hi=H5.get(f'{pfx}.main.{var}.hi',H5.get(f'{pfx}.{var}.hi')),
                       P=H5.get(f'{pfx}.main.{var}.P',H5.get(f'{pfx}.{var}.P')),n=H5.get(f'{pfx}.main.{var}.n',H5.get(f'{pfx}.{var}.n')),group=grp))
    add('Pooled non-EUR (H5)','primary')
    for a,l in [('afr','African'),('csa','Central/South Asian'),('eas','East Asian'),('mid','Middle Eastern'),('amr','Admixed American')]: add(l,f'ancestry_{a}',grp='ancestry')
    add('Wrong-gene eQTL (neg. control)','negctrl_wrong','W_gtex','control')
    add('PAV− stratum','PAV0_stratum',grp='sens'); add('PAV+ stratum','PAV1_stratum',grp='sens')
    add('|ΔMAF| ≤ 0.1','excl_mafdiff_gt0.1',grp='sens'); add('EUR lead variant only','lead_only_variant',grp='sens')
    r=sens[(sens.scenario=='lead_cs_per_gene')&(sens.estimate=='H5_Aeqtl')].iloc[0]; fr.append(dict(analysis='Lead CS per gene',OR=r.OR,lo=r.lo,hi=r.hi,P=r.P,n=int(r.n),group='sens'))
    pv=H5['PAV_term.PAV.OR']; fr.append(dict(analysis='PAV itself (same model)',OR=pv,lo=H5['PAV_term.PAV.lo'],hi=H5['PAV_term.PAV.hi'],P=H5['PAV_term.PAV.P'],n=H5['PAV_term.PAV.n'],group='pav'))
    fr=pd.DataFrame(fr); ys=[]; gap={'main':0,'ancestry':1,'control':2,'sens':3,'pav':4}
    for i,r in fr.iterrows():
        y=-(i+gap[r.group]*0.6); c={'main':L_XA,'ancestry':L_XA,'control':GREY,'sens':L_XA,'pav':C_PAVNA}[r.group]
        mfc='white' if r.group=='sens' else c
        ax.plot([r.lo,r.hi],[y,y],color=c,lw=1.2); ax.plot(r.OR,y,'o',ms=4.5 if r.group!='main' else 6,mfc=mfc,mec=c,mew=1.0,zorder=3)
        ax.text(1.02,y,f'{r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})  n={int(r.n):,}',transform=ax.get_yaxis_transform(),va='center',fontsize=6); ys.append(y)
    ax.set_yticks(ys); ax.set_yticklabels(fr.analysis)
    ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2)))
    ax.set_xscale('log'); ax.set_xlim(0.45,3.0); ax.set_xticks([0.5,1,1.5,2,3]); ax.set_xticklabels(['0.5','1','1.5','2','3']); ax.minorticks_off()
    ax.set_xlabel('OR for cross-ancestry replication (95% CI)'); ax.set_ylim(min(ys)-0.6,0.6)
    ax.set_title('Anchoring predicts non-EUR replication, but so does PAV',loc='left',fontsize=8,x=-0.68)
    panel_letter(ax,'d',dx=-0.74)
    SRC['Fig5d']=fr
    ax=fig.add_axes([0.83,0.07,0.16,0.38]); rt=pd.DataFrame(json.loads(H5['primary.rates']))
    lab={(0,0):'Neither',(0,1):'PAV\nonly',(1,0):'Anch.\nonly',(1,1):'Anch.\n+PAV'}; col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}; rows_=[]
    for j,r in rt.iterrows():
        k=int(round(r['mean']*r['size'])); n=int(r['size']); lo_,hi_=wilson(k,n)
        ax.bar(j,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.75); ax.errorbar(j,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
        ax.text(j,0.02,f'{n:,}',rotation=90,ha='center',va='bottom',fontsize=6,color='white' if (r.A_eqtl,r.PAV)!=(0,0) else '#222222')
        rows_.append(dict(A_eqtl=r.A_eqtl,PAV=r.PAV,group=lab[(r.A_eqtl,r.PAV)].replace('\n',' '),k=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
    ax.set_xticks(range(4)); ax.set_xticklabels([lab[(a,p)] for a,p in zip(rt.A_eqtl,rt.PAV)],fontsize=6); ax.set_ylim(0,1)
    ax.set_ylabel('Cross-ancestry replication'); ax.set_title('PAV-only signals replicate\nas often as anchored ones',loc='left',fontsize=8,x=-0.35)
    panel_letter(ax,'e',dx=-0.45)
    SRC['Fig5e']=pd.DataFrame(rows_)
    fig.canvas.draw(); save(fig,'Fig5_grade_validation_cross_ancestry'); return fig


def fig6():
    fig=plt.figure(figsize=(7.2,6.2))
    ax=fig.add_axes([0.10,0.07,0.36,0.88])
    sm=asum.copy(); sm['tot']=sm[['n_T1','n_T2','n_T3','n_T4']].sum(1); sm=sm.sort_values('tot').reset_index(drop=True)
    y=np.arange(len(sm))
    left=np.zeros(len(sm))
    for t in [1,2,3,4]:
        ax.barh(y,sm[f'n_T{t}'],left=left,color=TIER_C[t],height=0.72,label=TIER_L[t]); left+=sm[f'n_T{t}']
    ax.scatter(np.full(len(sm),-2.0),y,marker='s',s=22,c=[TIER_C[int(t)] for t in sm.lead_cs_tier],clip_on=False,zorder=3)
    ax.text(-2.0,len(sm)-0.3,'lead\nCS',ha='center',va='bottom',fontsize=6)
    labs=[f'{s}*' if n>1 else s for s,n in zip(sm.symbol,sm.n_cs_tied_at_top_P)]
    ax.set_yticks(y); ax.set_yticklabels(labs,fontstyle='italic',fontsize=6); ax.tick_params(axis='y',pad=14,length=0)
    ax.set_xlim(-3.2,50); ax.set_xticks([0,10,20,30,40,50]); ax.set_ylim(-0.7,len(sm)+0.6); ax.set_xlabel('cis credible sets per target (UKB-PPP)')
    ax.spines['left'].set_visible(False)
    ax.legend(loc='center right',bbox_to_anchor=(1.0,0.30),fontsize=6,borderaxespad=0.2,handlelength=1)
    ax.text(49,4.0,'* ≥2 CS tied at the capped P;\nlead chosen by larger |β|',ha='right',fontsize=6)
    ax.set_title('Instrument tiers for 42 cardiometabolic targets',loc='left',fontsize=8,x=-0.2,pad=16)
    panel_letter(ax,'a',dx=-0.26,dy=1.0)
    SRC['Fig6a']=sm[['symbol','n_T1','n_T2','n_T3','n_T4','tot','lead_cs_tier','n_cs_tied_at_top_P','tiers_among_tied_top_cs']].rename(columns={'tot':'n_cs'})
    # b evidence matrix for T4-lead targets
    axb=fig.add_axes([0.60,0.43,0.38,0.50])
    tg=list(sm[sm.lead_cs_tier==4].symbol)[::-1]
    L=lead.loc[tg]; B=best.loc[tg]
    cols=['eQTL anchor','isoform anchor','PAV in 95% CS','lead is PAV','SomaScan\nconcordant','MS (MSPA)','non-EUR\nreplicated']
    def cell(r,c):
        if c=='eQTL anchor': return ('no',C_NONE) if r.anchor_eqtl==0 else ('yes',C_ANCH)
        if c=='isoform anchor': return ('yes',C_ANCH) if r.anchor_isoform==1 else ('no',C_NONE)
        if c=='PAV in 95% CS': return ('yes',C_PAVNA) if r.pav_cognate_cs95_vep116_alltx==1 else ('no',C_NONE)
        if c=='lead is PAV': return ('yes',C_PAVNA) if r.pav_lead_variant==1 else ('no',C_NONE)
        if c.startswith('SomaScan'): return {'concordant':('yes',L_AFF),'discordant':('no',C_NONE),'na':('n.a.','white')}[r.platform_concordance]
        if c.startswith('MS'): return {'supported':('supp.',L_MS),'refuted':('ref.','#7A7A7A'),'not_assessed':('n.a.','white'),'indeterminate_underpowered':('indet.','white')}[r.ms_primary]
        if c.startswith('non-EUR'): 
            v=r.xa_replicated; return ('n.a.','white') if pd.isna(v) else (('yes',L_XA) if v==1 else ('no',C_NONE))
    rec=[]
    for i,s_ in enumerate(tg):
        r=L.loc[s_]
        for j,c in enumerate(cols):
            t,fcol=cell(r,c); axb.add_patch(Rectangle((j,i),0.94,0.9,fc=fcol,ec='#BBBBBB',lw=0.5))
            axb.text(j+0.47,i+0.45,t,ha='center',va='center',fontsize=6,color='white' if fcol in (C_ANCH,C_PAVNA,L_AFF,L_MS,L_XA,'#7A7A7A') else '#222222')
        bt=B.loc[s_,'best_T1_mlog10p']; txt='none' if pd.isna(bt) else f'{bt:.0f}'
        axb.text(len(cols)+0.15,i+0.45,txt,ha='left',va='center',fontsize=6)
        rec.append(dict(symbol=s_,lead_variant=r.lead_variant_id,lead_mlog10p=r.lead_mlog10p,lead_beta=r.lead_beta,tier=int(r.tier),anchor_eqtl=r.anchor_eqtl,anchor_isoform=r.anchor_isoform,pav_cs95=r.pav_cognate_cs95_vep116_alltx,pav_lead=r.pav_lead_variant,platform_concordance=r.platform_concordance,ms_primary=r.ms_primary,xa_replicated=r.xa_replicated,flag_concordant_unanchored_PAV=r.flag_concordant_unanchored_PAV,best_T1_cs_mlog10p=bt,best_T1_lead_variant=B.loc[s_,'best_T1_variant']))
    axb.text(len(cols)+0.15,len(tg)+0.1,'best T1 CS\n−log10 P',ha='left',va='bottom',fontsize=6)
    axb.set_xlim(0,len(cols)+1.6); axb.set_ylim(0,len(tg)); axb.set_axis_off()
    for j,c in enumerate(cols): axb.text(j+0.47,len(tg)+0.1,c,ha='center',va='bottom',fontsize=6,rotation=90)
    for i,s_ in enumerate(tg): axb.text(-0.1,i+0.45,s_,ha='right',va='center',fontsize=6,fontstyle='italic')
    axb.set_title('Lead instruments graded T4: evidence available',loc='left',fontsize=8,pad=58)
    panel_letter(axb,'b',dx=-0.12,dy=1.22)
    SRC['Fig6b']=pd.DataFrame(rec)
    # c note box
    axc=fig.add_axes([0.60,0.07,0.38,0.27]); axc.set_axis_off()
    nT1=(asum.lead_cs_tier==1).sum(); nT4=(asum.lead_cs_tier==4).sum(); nalt=int(B.best_T1_mlog10p.notna().sum())
    axc.text(0,1,f'Lead CS tier across 42 targets: T1 {nT1}, T2 {(asum.lead_cs_tier==2).sum()}, T3 {(asum.lead_cs_tier==3).sum()}, T4 {nT4}.\n'
             f'{nalt} of {nT4} T4-lead targets also carry a T1 credible set\nthat could serve as an alternative instrument.\n\n'
             'T4 does not mean artefact: functional PAVs (e.g. in PCSK9,\nAPOE, LPL) are graded T4 because affinity data alone cannot\nseparate abundance from epitope effects.',fontsize=6,va='top',ha='left',linespacing=1.3)
    panel_letter(axc,'c',dx=-0.12,dy=1.0)
    SRC['Fig6c']=pd.DataFrame([dict(lead_T1=nT1,lead_T2=(asum.lead_cs_tier==2).sum(),lead_T3=(asum.lead_cs_tier==3).sum(),lead_T4=nT4,T4_lead_with_alternative_T1=nalt)])
    fig.canvas.draw(); save(fig,'Fig6_drug_target_instrument_grades'); return fig


def figS1():
    eps=[('H1_G1_Agtex','H1 Sun 2018:\nGTEx anchor',L_AFF),('H2a_Agtex','H2a deCODE:\nGTEx anchor',L_AFF),('H2b_PAV_matched','H2b deCODE:\nPAV, matched',L_AFF),('H3i_Aeqtl','H3i MS:\neQTL anchor',L_MS),('H5_Aeqtl','H5 non-EUR:\neQTL anchor',L_XA)]
    scn=['primary','h4_0.9','h4_0.95','prior_p12_5e-6','prior_p12_1e-6','coloc_intersection','window_lead_TSS_le_1Mb','window_lead_TSS_le_500kb','lead_cs_per_gene','exclude_HLA','exclude_trans_hotspot','PAV_canonical','PAV_ppsum_ge_0.2','PAV_lead_variant']
    snl={'primary':'Primary analysis','h4_0.9':'PP.H4 ≥ 0.9','h4_0.95':'PP.H4 ≥ 0.95','prior_p12_5e-6':'p12 = 5×10⁻⁶','prior_p12_1e-6':'p12 = 1×10⁻⁶','coloc_intersection':'coloc on shared variants','window_lead_TSS_le_1Mb':'lead ≤ 1 Mb of TSS','window_lead_TSS_le_500kb':'lead ≤ 500 kb of TSS','lead_cs_per_gene':'lead CS per gene','exclude_HLA':'excluding HLA (vacuous)','exclude_trans_hotspot':'excluding trans hotspots','PAV_canonical':'PAV: canonical transcript','PAV_ppsum_ge_0.2':'PAV: summed PIP ≥ 0.2','PAV_lead_variant':'PAV: lead variant'}
    fig,axs=plt.subplots(1,5,figsize=(7.2,3.6),sharey=True,gridspec_kw=dict(wspace=0.08,left=0.20,right=0.99,top=0.86,bottom=0.13))
    y=np.arange(len(scn))[::-1]; recs=[]
    for ax,(e,lab,c) in zip(axs,eps):
        for yi,sc in zip(y,scn):
            r=sens[(sens.scenario==sc)&(sens.estimate==e)].iloc[0]
            ax.plot([r.lo,r.hi],[yi,yi],color=c,lw=1.1); ax.plot(r.OR,yi,'o',ms=3.5 if sc!='primary' else 5,color=c,mec='white',mew=0.4)
            recs.append(dict(endpoint=e,scenario=sc,scenario_label=snl[sc],OR=r.OR,lo=r.lo,hi=r.hi,P=r.P,n=r.n,n_cs_scenario=r.n_cs_scenario))
        ax.axvline(1,color='#999999',lw=0.6,ls=(0,(3,2))); ax.set_xscale('log'); ax.minorticks_off()
        ax.set_xticks([0.5,1,2,4]); ax.set_xticklabels(['0.5','1','2','4']); ax.set_xlim(0.3,6.5)
        ax.set_title(lab,loc='left',fontsize=7)
    axs[0].set_yticks(y); axs[0].set_yticklabels([snl[s] for s in scn])
    fig.supxlabel('Odds ratio (95% CI), log scale',fontsize=8,y=0.01)
    fig.text(0.01,0.95,'Estimates are stable across pre-specified sensitivity analyses except the one-CS-per-gene restriction',fontsize=8)
    SRC['FigS1']=pd.DataFrame(recs)
    fig.canvas.draw(); save(fig,'FigS1_sensitivity_analyses'); plt.close(fig)


def figS2():
    fig=plt.figure(figsize=(7.2,2.9)); recs={}
    ax=fig.add_axes([0.07,0.17,0.22,0.68]); tc=R.tier.value_counts().sort_index()
    ax.bar(tc.index,tc.values,color=[TIER_C[t] for t in tc.index],width=0.7)
    for t,v in tc.items(): ax.text(t,v+150,f'{v:,}',ha='center',fontsize=6)
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.set_ylabel('Credible sets'); ax.set_ylim(0,11200)
    ax.set_title('Tier distribution (n = 14,472)',loc='left',fontsize=8); panel_letter(ax,'a',dx=-0.35)
    recs['FigS2a']=tc.rename('n_cs').reset_index().rename(columns={'index':'tier'})
    ax=fig.add_axes([0.40,0.17,0.25,0.68])
    ct=pd.crosstab(R.platform_concordance,R.tier); ct=ct.loc[['concordant','discordant','na']]; pct=ct.div(ct.sum(1),axis=0)
    left=np.zeros(3)
    for t in [1,2,3,4]:
        ax.barh(range(3),pct[t],left=left,color=TIER_C[t],height=0.65); left+=pct[t]
    ax.set_yticks(range(3)); ax.set_yticklabels([f'{i}\n(n={ct.loc[i].sum():,})' for i in ['concordant','discordant','not assayed']] if False else [f'Concordant\n(n={ct.loc["concordant"].sum():,})',f'Discordant\n(n={ct.loc["discordant"].sum():,})',f'Not assayed\n(n={ct.loc["na"].sum():,})'])
    ax.invert_yaxis(); ax.set_xlim(0,1); ax.set_xlabel('Fraction of credible sets')
    t4c=ct.loc['concordant',4]; ax.text(1-pct.loc['concordant',4]/2,0,f'{t4c}',ha='center',va='center',fontsize=6,color='white')
    ax.set_title('Tier by SomaScan concordance',loc='left',fontsize=8); panel_letter(ax,'b',dx=-0.42)
    recs['FigS2b']=ct.reset_index().merge(pct.add_prefix('frac_T').reset_index(),on='platform_concordance')
    ax=fig.add_axes([0.75,0.17,0.23,0.68])
    data=[R.loc[R.tier==t,'lead_mlog10p'].clip(upper=300).values for t in [1,2,3,4]]
    bp=ax.boxplot(data,positions=[1,2,3,4],widths=0.6,patch_artist=True,showfliers=False,medianprops=dict(color='#222222',lw=1))
    for p_,t in zip(bp['boxes'],[1,2,3,4]): p_.set_facecolor(TIER_C[t]); p_.set_edgecolor('#444444')
    ax.set_yscale('log'); ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.set_ylabel('Lead −log10 P (capped at 300)')
    ax.set_yticks([10,30,100,300]); ax.set_yticklabels(['10','30','100','300']); ax.minorticks_off()
    ax.set_title('Signal strength differs by tier',loc='left',fontsize=8); panel_letter(ax,'c',dx=-0.38)
    recs['FigS2c']=R.groupby('tier').lead_mlog10p.describe().reset_index()
    SRC.update(recs); fig.canvas.draw(); save(fig,'FigS2_resource_overview'); plt.close(fig)


def figS3():
    feats=['A_eqtl','log1p_nanchor','A_iso','PAV','pav_pp','PAVxNoAnchor','conc_concordant','conc_discordant','mlp','absb','maf','logn95','purityMinR2','logdist']
    fl={'A_eqtl':'eQTL anchor','log1p_nanchor':'log(1 + anchor datasets)','A_iso':'isoform anchor','PAV':'PAV in 95% CS','pav_pp':'summed PAV PIP','PAVxNoAnchor':'PAV × no anchor','conc_concordant':'SomaScan concordant','conc_discordant':'SomaScan discordant','mlp':'−log10 P (capped)/10','absb':'|lead β|','maf':'MAF','logn95':'log 95% CS size','purityMinR2':'CS purity (min r²)','logdist':'log distance to TSS'}
    fig=plt.figure(figsize=(7.2,3.4))
    ax=fig.add_axes([0.20,0.13,0.38,0.75]); y=np.arange(len(feats))[::-1]; rec=[]
    for yi,f_ in zip(y,feats):
        a=H4[f'H4_ms.coef_full_std.{f_}']; b=H4[f'H4_xa.coef_full_std.{f_}']
        ax.plot([a,b],[yi,yi],color='#DDDDDD',lw=0.8,zorder=1)
        ax.plot(a,yi,'o',ms=4.5,color=L_MS,zorder=3); ax.plot(b,yi,'s',ms=4,color=L_XA,zorder=3)
        rec.append(dict(feature=f_,label=fl[f_],coef_std_MS_target=a,coef_std_crossancestry_target=b))
    ax.axvline(0,color='#999999',lw=0.6); ax.set_yticks(y); ax.set_yticklabels([fl[f] for f in feats])
    ax.set_xlabel('Standardized L2-logistic coefficient (training chromosomes)'); ax.set_xlim(-1.15,1.55)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=L_MS,label='MS support (n = 99)'),Line2D([],[],marker='s',ls='',color=L_XA,label='cross-ancestry (n = 4,572)')],loc='lower right',fontsize=6,borderaxespad=0.1,title='training target',title_fontsize=6)
    ax.set_title('Signal strength has opposite weights for the two targets',loc='left',fontsize=8,x=-0.45)
    panel_letter(ax,'a',dx=-0.50)
    SRC['FigS3a']=pd.DataFrame(rec)
    ax=fig.add_axes([0.72,0.13,0.26,0.75]); rows_=[]
    items=[('MS target, held-out MS',H4['H4_ms.auc_full'],L_MS),('→ held-out cross-ancestry',H4['transfer_ms_to_xa.auc1'],L_XA),('→ MSPA peptide label',H4['transfer_ms_to_ms_mspa.auc1'],L_MS),('→ Seer discovery overlap',H4['transfer_ms_to_ms_seer.auc1'],L_MS)]
    for i,(l,v,c) in enumerate(items):
        yy=len(items)-1-i; ax.plot([0.5,v],[yy,yy],color='#BBBBBB',lw=0.9); ax.plot(v,yy,'o',ms=5,color=c); ax.text(v+(0.03 if v>0.5 else -0.03),yy,f'{v:.2f}',va='center',ha='left' if v>0.5 else 'right',fontsize=6)
        rows_.append(dict(evaluation=l,AUC=v))
    ax.axvline(0.5,color='#999999',lw=0.6,ls=(0,(3,2))); ax.set_xlim(0,1); ax.set_yticks(range(len(items))[::-1]); ax.set_yticklabels([i[0] for i in items])
    ax.set_xlabel('AUC (0.5 = chance)'); ax.set_title('MS-trained model does not transfer',loc='left',fontsize=8,x=-0.6)
    panel_letter(ax,'b',dx=-0.68)
    SRC['FigS3b']=pd.DataFrame(rows_)
    fig.canvas.draw(); save(fig,'FigS3_grade_model_coefficients_transfer'); plt.close(fig)


def figS4():
    fig=plt.figure(figsize=(7.2,3.0)); rec=[]
    anc=[('ancestry_afr','African'),('ancestry_csa','Central/South Asian'),('ancestry_eas','East Asian'),('ancestry_mid','Middle Eastern'),('ancestry_amr','Admixed American')]
    col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    for i,(k,l) in enumerate(anc):
        ax=fig.add_axes([0.07+i*0.108,0.20,0.092,0.62]); rt=pd.DataFrame(json.loads(H5[f'{k}.rates']))
        for j,r in rt.iterrows():
            n=int(r['size']); kk=int(round(r['mean']*n)); lo_,hi_=wilson(kk,n)
            ax.bar(j,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.8); ax.errorbar(j,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.6,capsize=1)
            rec.append(dict(ancestry=l,A_eqtl=r.A_eqtl,PAV=r.PAV,k=kk,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
        ax.set_ylim(0,0.8); ax.set_xticks([]); ax.set_title(l.replace('/','/\n') if '/' in l else l.replace(' ','\n'),loc='left',fontsize=7)
        if i==0: ax.set_ylabel('Replication rate (same sign, P < 0.05)'); panel_letter(ax,'a',dx=-0.75)
        else: ax.set_yticklabels([])
    fig.legend(handles=[Patch(color=col[k],label=n_) for k,n_ in [((0,0),'Neither'),((0,1),'PAV only'),((1,0),'Anchor only'),((1,1),'Anchor + PAV')]],loc='lower left',bbox_to_anchor=(0.06,0.0),ncol=4,fontsize=6,handlelength=1)
    SRC['FigS4a']=pd.DataFrame(rec)
    ph=xl['PostHoc_exploratory']; ph=ph[ph.analysis.str.startswith('Cross-ancestry slope')].copy(); ph['label']=ph.analysis.str.replace('Cross-ancestry slope: ','')
    lbl={'WLS_all(primary secondary)':'All (pre-specified)','OLS_unweighted':'Unweighted OLS','WLS_|eurb|<=1':'|EUR β| ≤ 1','WLS_noPAV':'PAV− only','WLS_PAV':'PAV+ only','WLS_African':'African only','WLS_Central_South Asian':'C/S Asian only','WLS_excl_beta_inconsistent':'excl. 204 β-inconsistent CS'}
    ax=fig.add_axes([0.765,0.20,0.215,0.62]); y=np.arange(len(ph))[::-1]
    for yi,(_,r) in zip(y,ph.iterrows()):
        ax.plot([r.slope_unanchored,r.slope_anchored],[yi,yi],color='#DDDDDD',lw=0.8)
        ax.plot(r.slope_unanchored,yi,'o',ms=4,mfc='white',mec=L_XA,mew=1); ax.plot(r.slope_anchored,yi,'o',ms=4,color=L_XA)
    ax.axvline(1,color='#999999',lw=0.6,ls=(0,(3,2))); ax.set_yticks(y); ax.set_yticklabels([lbl.get(l,l) for l in ph.label]); ax.set_xlabel('non-EUR β on EUR β slope',labelpad=1)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',mec=L_XA,label='unanchored'),Line2D([],[],marker='o',ls='',color=L_XA,label='anchored')],loc='upper left',bbox_to_anchor=(0.0,-0.20),ncol=2,fontsize=6,borderaxespad=0.1)
    ax.set_title('Slope gap after β-QC (post hoc)',loc='left',fontsize=8,x=-0.6); panel_letter(ax,'b',dx=-0.68)
    SRC['FigS4b']=ph.assign(label=[lbl.get(l,l) for l in ph.label])[['label','slope_unanchored','slope_anchored','P','n']].rename(columns={'P':'P_interaction','n':'n_rows'})
    fig.canvas.draw(); save(fig,'FigS4_cross_ancestry_by_ancestry_slopes'); plt.close(fig)


def figS5():
    fig=plt.figure(figsize=(7.2,2.8)); rec=[]
    for gi,(lab,nm_) in enumerate([('primary','Primary label (Bonferroni MS replication)'),('mspa','MSPA peptide label')]):
        ax=fig.add_axes([0.07+gi*0.49,0.20,0.43,0.62]); gd=pd.DataFrame(json.loads(H3[f'{lab}.grid']))
        col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}; x=0; xt=[];xl_=[]
        for cc,cn in [('concordant','Concordant'),('discordant','Discordant'),('na','Not assayed')]:
            sub=gd[gd.conc_cat==cc]
            for _,r in sub.iterrows():
                k=int(r['sum']); n=int(r['size']); lo_,hi_=wilson(k,n)
                ax.bar(x,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.8); ax.errorbar(x,r['mean'],yerr=[[max(r['mean']-lo_,0)],[max(hi_-r['mean'],0)]],fmt='none',ecolor='#222222',lw=0.6,capsize=1)
                (ax.plot(x,0.008,marker='_',ms=8,mew=2,color=col[(r.A_eqtl,r.PAV)],clip_on=False) if k==0 else None)
                ax.text(x,hi_+0.02,f'{k}/{n}',ha='center',va='bottom',fontsize=6,rotation=90)
                rec.append(dict(label=lab,concordance=cc,A_eqtl=r.A_eqtl,PAV=r.PAV,supported=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_)); x+=1
            xt.append(x-2.5); xl_.append(cn); x+=0.8
        ax.set_xticks(xt); ax.set_xticklabels(xl_); ax.set_ylim(0,1.25); ax.set_yticks([0,0.25,0.5,0.75,1.0])
        ax.set_title(nm_,loc='left',fontsize=8)
        if gi==0: ax.set_ylabel('MS-supported fraction')
        panel_letter(ax,'ab'[gi],dx=-0.12)
    fig.legend(handles=[Patch(color=col[k],label=n_) for k,n_ in [((0,0),'Neither'),((0,1),'PAV only'),((1,0),'Anchor only'),((1,1),'Anchor + PAV')]],loc='lower center',ncol=4,fontsize=6,handlelength=1,bbox_to_anchor=(0.5,0.0))
    SRC['FigS5']=pd.DataFrame(rec)
    fig.canvas.draw(); save(fig,'FigS5_MS_support_full_grid'); plt.close(fig)

# ===== v1 overrides =====
# --- fig1
def fig1():
    fig=plt.figure(figsize=(7.2,6.0))
    ax=fig.add_axes([0.01,0.53,0.98,0.44]); ax.set_axis_off(); ax.set_xlim(0,1); ax.set_ylim(0,1)
    for x,t in [(0.0,'Discovery'),(0.175,'Evidence per credible set'),(0.425,'Evaluation layers (pre-registered)'),(0.785,'Resource')]:
        ax.text(x+0.004,1.0,t,ha='left',va='top',fontsize=6,color='#555555',style='italic',transform=ax.transAxes)
    H=0.16; gap=0.02; tops=[0.92-H-i*(H+gap) for i in range(5)]
    box(ax,0.0,tops[4],0.16,0.92-tops[4],'UKB-PPP Olink','European ancestry\n14,472 cis credible\nsets (CS)\n1,879 genes\n\nOpen Targets 26.09\nSuSiE-inf, 95% sets',L_DISC)
    fh=(0.92-tops[4]-2*gap)/3
    fy=[0.92-fh, 0.92-2*fh-gap, tops[4]]
    box(ax,0.175,fy[0],0.23,fh,'Transcript anchor','coloc PP.H4 ≥ 0.8 with\ncognate eQTL (GTEx v10;\neQTL Catalogue)',C_ANCH)
    box(ax,0.175,fy[1],0.23,fh,'Protein-altering variant','cognate missense/LoF/splice\nvariant in 95% CS (VEP 116)',C_PAVNA)
    box(ax,0.175,fy[2],0.23,fh,'Negative control','coloc with non-cognate\n(wrong-gene) GTEx eQTL',GREY)
    ev=[('SomaScan v3 (Sun 2018)','424 genes · 4,153 CS · H1 (internal)',L_AFF,'white','-'),
        ('SomaScan v4: deCODE, ARIC, Fenland','deCODE 888 new genes · 6,711 CS · H2\nARIC EA/AA, Fenland · A2 (same platform)',L_AFF,'white','-'),
        ('Mass spectrometry: MSPA','Suhre 2025 · 169 evaluable CS · H3, H4',L_MS,'white','-'),
        ('Mass spectrometry: Seer full','addendum A1 · 564 studies · 367 evaluable CS',L_MS,'white','-'),
        ('Olink, non-EUR ancestries','UKB-PPP AFR/CSA/EAS/MID/AMR\n· 8,065 CS · H5',L_XA,'white','-')]
    for t,(ti,bo,c,f,ls) in zip(tops,ev): box(ax,0.425,t,0.34,H,ti,bo,c,fc_=f,ls=ls)
    box(ax,0.785,tops[4],0.21,0.92-tops[4],'PAV-aware grades','Tier 1–4 for all\n14,472 CS\n\nT1 anchored, no PAV\nT2 anchored + PAV,\n     or isoform-only\nT3 no anchor, no PAV\nT4 no anchor + PAV\n\nCaution flag:\nconcordant + unanchored\n+ PAV (339 CS)',C_PAVNA)
    mid=(0.92+tops[4])/2
    arrow(ax,0.16,mid,0.175,mid); arrow(ax,0.405,mid,0.425,mid); arrow(ax,0.765,mid,0.785,mid)
    # layer-type key
    for i,(c,l) in enumerate([(L_AFF,'affinity, other platform'),(L_MS,'antibody-independent (MS)'),(L_XA,'same platform, other ancestry')]):
        ax.add_patch(Rectangle((0.175+i*0.21,-0.045),0.012,0.03,color=c,transform=ax.transAxes,clip_on=False))
        ax.text(0.175+i*0.21+0.017,-0.03,l,fontsize=6,va='center',transform=ax.transAxes)
    fig.text(0.005,0.974,'a',fontsize=10,fontweight='bold',va='bottom')
    # ---- b
    axb=fig.add_axes([0.385,0.075,0.29,0.36])
    yy=np.arange(len(rows))[::-1]
    for y,r in zip(yy,rows):
        h,lab,e,l,u,ph,thr,sgn,col,n=r
        axb.plot([l,u],[y,y],color=col,lw=1.4,solid_capstyle='butt')
        axb.plot(e,y,'o',ms=5.5,color=col,mec='white',mew=0.6,zorder=3)
        axb.plot(thr,y,marker='|',ms=9,color='#222222',mew=1.0,zorder=2)
    axb.axvline(1,color='#999999',lw=0.8,ls=(0,(3,2)),zorder=0)
    axb.set_xscale('log'); axb.set_xlim(0.035,6.5)
    axb.set_xticks([0.05,0.1,0.25,0.5,1,2,4]); axb.set_xticklabels(['0.05','0.1','0.25','0.5','1','2','4']); axb.minorticks_off()
    axb.set_yticks(yy); axb.set_yticklabels([f'{r[0]}  {r[1]}' for r in rows])
    axb.set_ylim(-1.2,len(rows)-0.4)
    axb.set_xlabel('Odds ratio (95% CI), log scale')
    axb.text(0.04,-0.95,'|  pre-specified threshold',fontsize=6,va='center')
    fig.text(0.02,0.475,'PAVs inflate affinity replication in three SomaScan v4 cohorts; anchoring predicts affinity and cross-ancestry, not MS, support',fontsize=8,ha='left')
    tx=fig.add_axes([0.69,0.075,0.30,0.36]); tx.set_axis_off(); tx.set_ylim(axb.get_ylim()); tx.set_xlim(0,1)
    for x,t in [(0.0,'OR (95% CI)'),(0.40,'Holm P'),(0.64,'Verdict')]: tx.text(x,len(rows)-0.5,t,fontsize=6,va='bottom',color='#555555')
    vs={'H2a':'supported','H2b':'supported','H2a-ARIC':'supported','H2b-ARIC':'supported','H2a-Fenland':'supported','H2b-Fenland':'supported','H5':'supported','H3i':'not supported','H3i-ext':'not supported','H3ii':'partial','H3ii-ext':'not supported'}
    for y,r in zip(yy,rows):
        tx.text(0.0,y,f'{r[2]:.2f} ({r[3]:.2f}–{r[4]:.2f})',fontsize=6,va='center')
        tx.text(0.40,y,fmt_p(r[5]),fontsize=6,va='center')
        tx.text(0.64,y,vs[r[0]],fontsize=6,va='center')
    fig.text(0.02,0.475+0.025,'b',fontsize=10,fontweight='bold',va='bottom')
    fig.canvas.draw(); save(fig,'Fig1_design_and_confirmatory_endpoints'); return fig




# --- fig4
def fig4():  # v1.2 Seer-only panel; the published Fig. 4 is figures/Fig4_mass_spec_four_layers.png (all four MS layers)
    fig=plt.figure(figsize=(7.2,6.6))
    # a: features x layers
    rows_=[]
    for li,(dct,pre,lab,nn) in enumerate([(H3,'primary.','MSPA primary label (n = 169)',169),(A3,'','Seer full, A1 (n = 367)',367)]):
        ax=fig.add_axes([0.075+li*0.235,0.64,0.20,0.29])
        feats=[('concordant_fisher','Platform\nconcordant',L_AFF),('H3i_fisher','eQTL\nanchored',C_ANCH),('PAV_fisher','PAV in\n95% CS',C_PAVNA)]
        for i,(k,fl,c) in enumerate(feats):
            t=tab(dct,f'{pre}{k}.table')
            for j,(sup,ns_,face,gl) in enumerate([(t[0][0],t[0][1],c,'yes'),(t[1][0],t[1][1],'#C8C8C8','no')]):
                n=sup+ns_; r=sup/n; lo_,hi_=wilson(sup,n); xx=i+(j-0.5)*0.38
                ax.bar(xx,r,width=0.36,color=face); ax.errorbar(xx,r,yerr=[[r-lo_],[hi_-r]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
                ax.text(xx,0.015,f'{sup}/{n}',ha='center',va='bottom',fontsize=6,color='white' if gl=='yes' else '#222222',rotation=90)
                rows_.append(dict(layer=lab,feature=fl.replace('\n',' '),status=gl,supported=sup,n=n,rate=r,wilson_lo=lo_,wilson_hi=hi_))
            ax.text(i,0.84,f'P = {fmt_p(dct[f"{pre}{k}.P"])}',ha='center',fontsize=6)
        ax.set_xticks(range(3)); ax.set_xticklabels([f[1] for f in feats],fontsize=6); ax.set_ylim(0,0.92)
        ax.set_title(lab,loc='left',fontsize=7)
        if li==0: ax.set_ylabel('MS-supported fraction'); panel_letter(ax,'a',dx=-0.30)
        else: ax.set_yticklabels([])
    fig.text(0.075,0.965,'MS-supported fraction by feature, Seer-study labels only (MSPA, Seer full); superseded by the four-layer Fig. 4 (figures/Fig4_mass_spec_four_layers.png)',fontsize=8)
    fig.legend(handles=[Patch(color='#555555',label='feature present (coloured)'),Patch(color='#C8C8C8',label='absent')],loc='upper left',bbox_to_anchor=(0.07,0.585),ncol=2,fontsize=6,handlelength=1) if False else None
    SRC['Fig4a']=pd.DataFrame(rows_); SRC['Fig4a']['test']='two-sided Fisher exact; coloured = feature present, grey = absent'
    # b: H3ii groups
    ax=fig.add_axes([0.62,0.64,0.37,0.29]); rows_=[]
    gc={'CA':C_ANCH,'CN0':C_NONE,'CNP':C_PAVNA}; gl={'CA':'Anchored','CN0':'No anchor\nno PAV','CNP':'No anchor\n+ PAV'}
    for i,(dct,key,nm_,pk) in enumerate([(H3,'primary.concordant_groups','MSPA (H3ii)','primary.H3ii_CNP_vs_CA'),(A3,'concordant_groups','Seer full (H3ii-ext)','H3ii_CNP_vs_CA')]):
        cg=pd.DataFrame(json.loads(dct[key])).set_index('grp')
        for j,gk in enumerate(['CA','CN0','CNP']):
            k=int(cg.loc[gk,'sum']); n=int(cg.loc[gk,'size']); r=k/n; lo_,hi_=wilson(k,n); xx=i*4.6+j*1.4
            ax.bar(xx,r,width=1.0,color=gc[gk]); ax.errorbar(xx,r,yerr=[[r-lo_],[hi_-r]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            ax.text(xx,hi_+0.02,f'{k}/{n}',ha='center',va='bottom',fontsize=6)
            rows_.append(dict(label=nm_,group=gk,group_label=gl[gk].replace('\n',' '),supported=k,n=n,rate=r,wilson_lo=lo_,wilson_hi=hi_))
        ax.text(i*4.6+1.4,1.02,f'{nm_}\nCNP vs CA: OR {float(dct[pk+".OR"]):.2f}, P = {fmt_p(float(dct[pk+".P"]))}',ha='center',fontsize=6)
    ax.set_xticks([0,1.4,2.8,4.6,6.0,7.4]); ax.set_xticklabels(['CA','CN0','CNP']*2,fontsize=6); ax.text(0.5,-0.13,'CA anchored · CN0 no anchor, no PAV · CNP no anchor + PAV',transform=ax.transAxes,ha='center',va='top',fontsize=6)
    ax.set_ylim(0,1.18); ax.set_ylabel('MS-supported fraction'); ax.set_title('Concordant signals: the MSPA hint is not replicated',loc='left',fontsize=7)
    panel_letter(ax,'b',dx=-0.14)
    SRC['Fig4b']=pd.DataFrame(rows_); SRC['Fig4b']['test']='two-sided Fisher exact CNP vs CA; Holm P: MSPA 0.28 (main family), Seer full 1.0 (addendum family)'
    # c: anchor forest
    ax=fig.add_axes([0.27,0.07,0.20,0.48]); fr=[]
    def add(l,OR,lo_,hi_,P,n,grp): fr.append(dict(analysis=l,OR=float(OR),lo=float(lo_),hi=float(hi_),P=float(P),n=int(n),group=grp))
    def addk(dct,pre,var,l,grp):
        add(l,dct[f'{pre}{var}.OR'],dct[f'{pre}{var}.lo'],dct[f'{pre}{var}.hi'],dct[f'{pre}{var}.P'],dct[f'{pre}{var}.n'],grp)
    addk(H3,'primary.','H3i_logit.A_eqtl','MSPA primary (H3i)','mspa')
    for k,l in [('mspa','MSPA peptide label'),('primary_ST10_only','UKB-PPP rows only (ST10)'),('primary_excl_pep_pav','excl. PAV-peptide rows'),('seer','Seer discovery overlap')]: addk(H3,k+'.','H3i_logit.A_eqtl',l,'mspa')
    addk(H3,'primary.','H3i_negctrl_wrong.W_gtex','wrong-gene eQTL (MSPA)','control')
    addk(A3,'','H3i_logit.A_eqtl','Seer full (H3i-ext)','seer')
    for k,l in [('H3_ext_new','ext-new (no MSPA overlap)'),('H3_ext_nominal','nominal label (P < 0.05)'),('H3_ext_cswide','CS-wide label')]: addk(A3s,k+'.','H3i_logit.A_eqtl',l,'seer')
    addk(A3,'','H3i_negctrl_wrong.W_gtex','wrong-gene eQTL (Seer full)','control')
    snm={'h4_0.9':'PP.H4 ≥ 0.9','h4_0.95':'PP.H4 ≥ 0.95','prior_p12_5e-6':'p12 = 5×10⁻⁶','prior_p12_1e-6':'p12 = 1×10⁻⁶','coloc_intersection':'coloc on shared variants','window_lead_TSS_le_1Mb':'lead ≤ 1 Mb of TSS','window_lead_TSS_le_500kb':'lead ≤ 500 kb of TSS','lead_cs_per_gene':'lead CS per gene','exclude_trans_hotspot':'excl. trans hotspots'}
    for sc,l in snm.items():
        r=sens[(sens.scenario==sc)&(sens.estimate=='H3i_Aeqtl')].iloc[0]; add(l,r.OR,r.lo,r.hi,r.P,r.n,'sens')
    fr=pd.DataFrame(fr); ys=[]; gap={'mspa':0,'control':0,'seer':1,'sens':2}
    y=0; prev=None; headers=[]
    for i,r in fr.iterrows():
        grpblk={'mspa':'A','control':None,'seer':'B','sens':'C'}[r.group]
        if grpblk and grpblk!=prev: y-=1.0; headers.append((y,{'A':'MSPA labels (Suhre 2025)','B':'Seer full summary statistics (A1)','C':'Sensitivity, MSPA primary'}[grpblk])); y-=1.0; prev=grpblk
        else: y-=1.0
        c=GREY if r.group=='control' else (L_MS)
        mfc='white' if r.group=='sens' else c
        ax.plot([r.lo,r.hi],[y,y],color=c,lw=1.1); ax.plot(r.OR,y,'o',ms=3.8,mfc=mfc,mec=c,mew=0.9,zorder=3)
        ax.text(1.02,y,f'{r.OR:.2f}  n={r.n:,}',transform=ax.get_yaxis_transform(),va='center',fontsize=6); ys.append(y)
    for yh,t in headers: ax.text(-0.03,yh,t,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6,fontweight='bold',color=L_MS)
    ax.set_yticks(ys); ax.set_yticklabels(fr.analysis)
    ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2))); ax.axvline(1.5,color='#222222',lw=0.6,ls=':')
    ax.set_xscale('log'); ax.set_xlim(0.2,10); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8']); ax.minorticks_off()
    ax.set_ylim(min(ys)-1.8,0.2); ax.set_xlabel('OR, eQTL anchor → MS support')
    ax.text(0.98,0.005,'dotted: threshold OR 1.5',transform=ax.transAxes,fontsize=6,ha='right',va='bottom')
    ax.set_title('Anchoring does not predict MS support',loc='left',fontsize=8,x=-1.05)
    panel_letter(ax,'c',dx=-1.12,dy=1.0)
    SRC['Fig4c']=fr; SRC['Fig4c']['note']='gene-clustered logistic ms ~ A_eqtl + mlp (Seer discovery additionally power-adjusted); Seer full = addendum A1'
    # d: PAV and concordance ORs across labels (Woolf CI from 2x2)
    ax=fig.add_axes([0.745,0.07,0.17,0.44]); rec=[]
    labs=[('MSPA primary',H3,'primary.'),('MSPA peptide',H3,'mspa.'),('MSPA, ST10 only',H3,'primary_ST10_only.'),('Seer full (A1)',A3,''),('Seer, ext-new',A3s,'H3_ext_new.'),('Seer, nominal',A3s,'H3_ext_nominal.'),('Seer, CS-wide',A3s,'H3_ext_cswide.')]
    yy=np.arange(len(labs))[::-1]*1.0
    for y_,(l,dct,pre) in zip(yy,labs):
        for off,k,c,nm_ in [(0.17,'PAV_fisher',C_PAVNA,'PAV'),(-0.17,'concordant_fisher',L_AFF,'concordant')]:
            t=tab(dct,f'{pre}{k}.table'); a_,b_,c_,d_=t[0][0],t[0][1],t[1][0],t[1][1]
            o,l_,h_=woolf(a_+0.5*(min(t[0]+t[1])==0),b_+0.5*(min(t[0]+t[1])==0),c_+0.5*(min(t[0]+t[1])==0),d_+0.5*(min(t[0]+t[1])==0))
            ax.plot([l_,h_],[y_+off]*2,color=c,lw=1.1); ax.plot(o,y_+off,'o',ms=3.8,color=c,zorder=3)
            rec.append(dict(label=l,feature=nm_,sup_with=a_,notsup_with=b_,sup_without=c_,notsup_without=d_,OR=o,woolf_lo=l_,woolf_hi=h_,fisher_P=float(dct[f'{pre}{k}.P'])))
    ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2))); ax.set_xscale('log'); ax.set_xlim(0.12,6); ax.set_xticks([0.25,0.5,1,2,4]); ax.set_xticklabels(['0.25','0.5','1','2','4']); ax.minorticks_off()
    ax.set_yticks(yy); ax.set_yticklabels([l for l,_,_ in labs]); ax.set_ylim(-0.7,len(labs)-0.3)
    ax.set_xlabel('OR for MS support (95% CI)')
    ax.legend(handles=[Line2D([],[],marker='o',color=C_PAVNA,label='PAV in 95% CS'),Line2D([],[],marker='o',color=L_AFF,label='platform concordant')],loc='upper center',bbox_to_anchor=(0.4,-0.12),fontsize=6,ncol=1,borderaxespad=0)
    ax.set_title('PAV signals are less often MS-supported\n(exploratory; both layers from one study)',loc='left',fontsize=7,x=-0.9)
    panel_letter(ax,'d',dx=-0.95,dy=1.13)
    SRC['Fig4d']=pd.DataFrame(rec); SRC['Fig4d']['note']='sample OR with Woolf 95% CI computed for display; P two-sided Fisher exact; MSPA and Seer full derive from the same MS study (label agreement kappa 0.89, n=148)'
    fig.canvas.draw(); save(fig,'FigS_v1.2_Seer_MS_layers_superseded'); return fig


# --- fig5
def fig5():
    fig=plt.figure(figsize=(7.2,6.4))
    # a: H5 forest
    ax=fig.add_axes([0.25,0.58,0.30,0.36]); fr=[]
    def add(l,pfx,var='A_eqtl',grp='main'):
        k=lambda s: H5.get(f'{pfx}.main.{var}.{s}',H5.get(f'{pfx}.{var}.{s}'))
        fr.append(dict(analysis=l,OR=k('OR'),lo=k('lo'),hi=k('hi'),P=k('P'),n=k('n'),group=grp))
    add('Pooled non-EUR (H5)','primary')
    for a_,l in [('afr','African'),('csa','Central/South Asian'),('eas','East Asian'),('mid','Middle Eastern'),('amr','Admixed American')]: add(l,f'ancestry_{a_}',grp='ancestry')
    add('Wrong-gene eQTL (neg. control)','negctrl_wrong','W_gtex','control')
    add('PAV− stratum','PAV0_stratum',grp='sens'); add('PAV+ stratum','PAV1_stratum',grp='sens')
    add('|ΔMAF| ≤ 0.1','excl_mafdiff_gt0.1',grp='sens'); add('EUR lead variant only','lead_only_variant',grp='sens')
    r=sens[(sens.scenario=='lead_cs_per_gene')&(sens.estimate=='H5_Aeqtl')].iloc[0]; fr.append(dict(analysis='Lead CS per gene',OR=r.OR,lo=r.lo,hi=r.hi,P=r.P,n=int(r.n),group='sens'))
    fr.append(dict(analysis='PAV itself (same model)',OR=H5['PAV_term.PAV.OR'],lo=H5['PAV_term.PAV.lo'],hi=H5['PAV_term.PAV.hi'],P=H5['PAV_term.PAV.P'],n=H5['PAV_term.PAV.n'],group='pav'))
    fr=pd.DataFrame(fr); ys=[]; gap={'main':0,'ancestry':1,'control':2,'sens':3,'pav':4}
    for i,r in fr.iterrows():
        y=-(i+gap[r.group]*0.6); c={'main':L_XA,'ancestry':L_XA,'control':GREY,'sens':L_XA,'pav':C_PAVNA}[r.group]
        ax.plot([r.lo,r.hi],[y,y],color=c,lw=1.2); ax.plot(r.OR,y,'o',ms=4.5 if r.group!='main' else 6,mfc='white' if r.group=='sens' else c,mec=c,mew=1.0,zorder=3)
        ax.text(1.02,y,f'{r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})',transform=ax.get_yaxis_transform(),va='center',fontsize=6); ys.append(y)
    ax.set_yticks(ys); ax.set_yticklabels(fr.analysis); ax.axvline(1,color='#999999',lw=0.7,ls=(0,(3,2)))
    ax.set_xscale('log'); ax.set_xlim(0.45,3.0); ax.set_xticks([0.5,1,1.5,2,3]); ax.set_xticklabels(['0.5','1','1.5','2','3']); ax.minorticks_off()
    ax.set_xlabel('OR, cross-ancestry replication'); ax.set_ylim(min(ys)-0.6,0.6)
    ax.set_title('Anchoring predicts non-EUR replication, but so does PAV',loc='left',fontsize=8,x=-0.8)
    panel_letter(ax,'a',dx=-0.86)
    SRC['Fig5a']=fr
    ax=fig.add_axes([0.83,0.58,0.16,0.36]); rt=pd.DataFrame(json.loads(H5['primary.rates'])); rows_=[]
    lab={(0,0):'Neither',(0,1):'PAV\nonly',(1,0):'Anch.\nonly',(1,1):'Anch.\n+PAV'}; col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    for j,r in rt.iterrows():
        n=int(r['size']); k=int(round(r['mean']*n)); lo_,hi_=wilson(k,n)
        ax.bar(j,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.75); ax.errorbar(j,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
        ax.text(j,0.02,f'{n:,}',rotation=90,ha='center',va='bottom',fontsize=6,color='white' if (r.A_eqtl,r.PAV)!=(0,0) else '#222222')
        rows_.append(dict(A_eqtl=r.A_eqtl,PAV=r.PAV,group=lab[(r.A_eqtl,r.PAV)].replace('\n',' '),k=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
    ax.set_xticks(range(4)); ax.set_xticklabels([lab[(a_,p)] for a_,p in zip(rt.A_eqtl,rt.PAV)],fontsize=6); ax.set_ylim(0,1)
    ax.set_ylabel('Cross-ancestry replication'); ax.set_title('PAV-only signals replicate\nas often as anchored ones',loc='left',fontsize=8,x=-0.35)
    panel_letter(ax,'b',dx=-0.45); SRC['Fig5b']=pd.DataFrame(rows_)
    # c-e ROC
    mstyle={'Full grade model':dict(color='#222222',lw=1.4),'Concordance + power':dict(color='#8A8A8A',lw=1.1),'Concordance only':dict(color='#BBBBBB',lw=1.1,ls='--')}
    pans=['Held-out MSPA support (n=70)','Held-out Seer-full support (n=139)','Held-out cross-ancestry (n=3,493)']
    aucs={pans[0]:{'Full grade model':H4['H4_ms.auc_full'],'Concordance + power':H4['H4_ms.auc_B1'],'Concordance only':H4['H4_ms.auc_B0']},
          pans[1]:{'Full grade model':A4['H4_ext.auc_full'],'Concordance + power':A4['H4_ext.auc_B1'],'Concordance only':A4['H4_ext.auc_B0']},
          pans[2]:{'Full grade model':H4['H4_xa.auc_full'],'Concordance + power':H4['H4_xa.auc_B1'],'Concordance only':H4['H4_xa.auc_B0']}}
    tit=['MS support, MSPA\n(H4; n = 70)','MS support, Seer full\n(H4-ext; n = 139)','Cross-ancestry\n(n = 3,493)']
    for i,p in enumerate(pans):
        ax=fig.add_axes([0.07+i*0.175,0.09,0.14,0.30])
        for m in ['Concordance only','Concordance + power','Full grade model']:
            d_=roc[(roc.panel==p)&(roc.model==m)]; ax.plot(d_.fpr,d_.tpr,**mstyle[m])
        ax.plot([0,1],[0,1],color='#DDDDDD',lw=0.6,zorder=0); ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.02); ax.set_aspect('equal')
        ax.set_xticks([0,0.5,1]); ax.set_yticks([0,0.5,1]); ax.set_xlabel('False-positive rate')
        if i==0: ax.set_ylabel('True-positive rate')
        else: ax.set_yticklabels([])
        ax.set_title(tit[i],loc='left',fontsize=7)
        for k_,m in enumerate(['Full grade model','Concordance + power','Concordance only']):
            ax.text(0.98,0.04+k_*0.10,f'{aucs[p][m]:.2f}',transform=ax.transAxes,ha='right',va='bottom',fontsize=6,color=mstyle[m]['color'] if m!='Concordance only' else '#999999')
        ax.text(0.98,0.04+3*0.10,'AUC',transform=ax.transAxes,ha='right',va='bottom',fontsize=6)
        panel_letter(ax,'cde'[i],dx=-0.25 if i==0 else -0.1)
    fig.legend(handles=[Line2D([],[],**mstyle[m],label=m) for m in ['Full grade model','Concordance + power','Concordance only']],loc='upper left',bbox_to_anchor=(0.06,0.465),ncol=3,fontsize=6,handlelength=1.8,columnspacing=1.0)
    SRC['Fig5cde']=roc.copy()
    SRC['Fig5cde_auc']=pd.DataFrame([dict(target=p,model=m,AUC=v) for p in aucs for m,v in aucs[p].items()]+[
        dict(target='MSPA (H4)',model='DeLong one-sided P full vs B0',AUC=H4['H4_ms.delong_full_vs_B0.P_one_sided']),
        dict(target='Seer full (H4-ext)',model='DeLong one-sided P full vs B0',AUC=A4['H4_ext.delong_full_vs_B0.P_one_sided']),
        dict(target='Seer full (H4-ext)',model='DeLong one-sided P full vs B1',AUC=A4['H4_ext.delong_full_vs_B1.P_one_sided']),
        dict(target='Seer full ext-new',model='dAUC full vs B0',AUC=A4['H4_ext_new.dauc']),
        dict(target='Seer full ext-new',model='DeLong one-sided P',AUC=A4['H4_ext_new.P_one_sided'])])
    # f tiers
    ax=fig.add_axes([0.665,0.09,0.325,0.30]); rows_=[]
    outs=[('MSPA',json.loads(H4['tiers_vs_truth.ms_primary.rates']),H4['tiers_vs_truth.ms_primary.trend.P_two_sided']),('Seer full',json.loads(A4['tiers_ext.rates']),A4['tiers_ext.trend.P_two_sided']),('Cross-ancestry',json.loads(H4['tiers_vs_truth.xa_rep.rates']),H4['tiers_vs_truth.xa_rep.trend.P_two_sided'])]
    for gi,(lab_,tr,pt) in enumerate(outs):
        for j,r in enumerate(tr):
            k=int(r['sum']); n=int(r['size']); lo_,hi_=wilson(k,n); xx=gi*5+j
            ax.bar(xx,r['mean'],width=0.8,color=TIER_C[int(r['tier'])]); ax.errorbar(xx,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
            ax.text(xx,0.02,f'{n:,}',ha='center',va='bottom',fontsize=6,rotation=90,color='white' if int(r['tier']) in (1,4) else '#222222')
            rows_.append(dict(outcome=lab_,tier=int(r['tier']),k=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
        ax.text(gi*5+1.5,0.93,f'trend P = {fmt_p(pt)}',ha='center',fontsize=6)
        ax.text(gi*5+1.5,-0.13,lab_,ha='center',va='top',fontsize=6,transform=ax.get_xaxis_transform())
    ax.set_xticks([0,1,2,3,5,6,7,8,10,11,12,13]); ax.set_xticklabels(['T1','T2','T3','T4']*3,fontsize=6)
    ax.set_ylim(0,1.0); ax.set_ylabel('Fraction positive'); ax.set_title('Tiers are not monotonic for any outcome',loc='left',fontsize=8)
    panel_letter(ax,'f',dx=-0.17)
    SRC['Fig5f']=pd.DataFrame(rows_); SRC['Fig5f']['note']='Cochran–Armitage trend P (two-sided); n inside bars'
    fig.canvas.draw(); save(fig,'Fig5_anchoring_ancestry_grade_validation'); return fig


# --- fig6
def fig6():
    fig=plt.figure(figsize=(7.2,6.2))
    ax=fig.add_axes([0.10,0.07,0.36,0.88])
    sm=asum.copy(); sm['tot']=sm[['n_T1','n_T2','n_T3','n_T4']].sum(1); sm=sm.sort_values('tot').reset_index(drop=True)
    y=np.arange(len(sm))
    left=np.zeros(len(sm))
    for t in [1,2,3,4]:
        ax.barh(y,sm[f'n_T{t}'],left=left,color=TIER_C[t],height=0.72,label=TIER_L[t]); left+=sm[f'n_T{t}']
    ax.scatter(np.full(len(sm),-2.0),y,marker='s',s=22,c=[TIER_C[int(t)] for t in sm.lead_cs_tier],clip_on=False,zorder=3)
    ax.text(-2.0,len(sm)-0.3,'lead\nCS',ha='center',va='bottom',fontsize=6)
    labs=[f'{s}*' if n>1 else s for s,n in zip(sm.symbol,sm.n_cs_tied_at_top_P)]
    ax.set_yticks(y); ax.set_yticklabels(labs,fontstyle='italic',fontsize=6); ax.tick_params(axis='y',pad=14,length=0)
    ax.set_xlim(-3.2,50); ax.set_xticks([0,10,20,30,40,50]); ax.set_ylim(-0.7,len(sm)+0.6); ax.set_xlabel('cis credible sets per target (UKB-PPP)')
    ax.spines['left'].set_visible(False)
    ax.legend(loc='center right',bbox_to_anchor=(1.0,0.30),fontsize=6,borderaxespad=0.2,handlelength=1)
    ax.text(49,4.0,'* ≥2 CS tied at the capped P;\nlead chosen by larger |β|',ha='right',fontsize=6)
    ax.set_title('Instrument tiers for 42 cardiometabolic targets',loc='left',fontsize=8,x=-0.2,pad=16)
    panel_letter(ax,'a',dx=-0.26,dy=1.0)
    SRC['Fig6a']=sm[['symbol','n_T1','n_T2','n_T3','n_T4','tot','lead_cs_tier','n_cs_tied_at_top_P','tiers_among_tied_top_cs']].rename(columns={'tot':'n_cs'})
    # b evidence matrix for T4-lead targets
    axb=fig.add_axes([0.58,0.43,0.40,0.50])
    tg=list(sm[sm.lead_cs_tier==4].symbol)[::-1]
    L=lead.loc[tg]; B=best.loc[tg]
    cols=['eQTL anchor','isoform anchor','PAV in 95% CS','lead is PAV','SomaScan\nconcordant','MS (MSPA)','MS (Seer full)','non-EUR\nreplicated']
    def cell(r,c):
        if c=='eQTL anchor': return ('no',C_NONE) if r.anchor_eqtl==0 else ('yes',C_ANCH)
        if c=='isoform anchor': return ('yes',C_ANCH) if r.anchor_isoform==1 else ('no',C_NONE)
        if c=='PAV in 95% CS': return ('yes',C_PAVNA) if r.pav_cognate_cs95_vep116_alltx==1 else ('no',C_NONE)
        if c=='lead is PAV': return ('yes',C_PAVNA) if r.pav_lead_variant==1 else ('no',C_NONE)
        if c.startswith('SomaScan'): return {'concordant':('yes',L_AFF),'discordant':('no',C_NONE),'na':('n.a.','white')}[r.platform_concordance]
        if c=='MS (Seer full)':
            v=r.ms_ext_primary; return {'supported':('supp.',L_MS),'refuted':('ref.','#7A7A7A'),'indeterminate':('indet.','white')}.get(v,('n.a.','white'))
        if c.startswith('MS'): return {'supported':('supp.',L_MS),'refuted':('ref.','#7A7A7A'),'not_assessed':('n.a.','white'),'indeterminate_underpowered':('indet.','white')}[r.ms_primary]
        if c.startswith('non-EUR'): 
            v=r.xa_replicated; return ('n.a.','white') if pd.isna(v) else (('yes',L_XA) if v==1 else ('no',C_NONE))
    rec=[]
    for i,s_ in enumerate(tg):
        r=L.loc[s_]
        for j,c in enumerate(cols):
            t,fcol=cell(r,c); axb.add_patch(Rectangle((j,i),0.94,0.9,fc=fcol,ec='#BBBBBB',lw=0.5))
            axb.text(j+0.47,i+0.45,t,ha='center',va='center',fontsize=6,color='white' if fcol in (C_ANCH,C_PAVNA,L_AFF,L_MS,L_XA,'#7A7A7A') else '#222222')
        bt=B.loc[s_,'best_T1_mlog10p']; txt='none' if pd.isna(bt) else f'{bt:.0f}'
        axb.text(len(cols)+0.15,i+0.45,txt,ha='left',va='center',fontsize=6)
        rec.append(dict(symbol=s_,lead_variant=r.lead_variant_id,lead_mlog10p=r.lead_mlog10p,lead_beta=r.lead_beta,tier=int(r.tier),anchor_eqtl=r.anchor_eqtl,anchor_isoform=r.anchor_isoform,pav_cs95=r.pav_cognate_cs95_vep116_alltx,pav_lead=r.pav_lead_variant,platform_concordance=r.platform_concordance,ms_primary=r.ms_primary,ms_ext_primary=r.ms_ext_primary,xa_replicated=r.xa_replicated,flag_concordant_unanchored_PAV=r.flag_concordant_unanchored_PAV,best_T1_cs_mlog10p=bt,best_T1_lead_variant=B.loc[s_,'best_T1_variant']))
    axb.text(len(cols)+0.15,len(tg)+0.1,'best T1 CS\n−log10 P',ha='left',va='bottom',fontsize=6)
    axb.set_xlim(0,len(cols)+1.6); axb.set_ylim(0,len(tg)); axb.set_axis_off()
    for j,c in enumerate(cols): axb.text(j+0.47,len(tg)+0.1,c,ha='center',va='bottom',fontsize=6,rotation=90)
    for i,s_ in enumerate(tg): axb.text(-0.1,i+0.45,s_,ha='right',va='center',fontsize=6,fontstyle='italic')
    axb.set_title('Lead instruments graded T4: evidence available',loc='left',fontsize=8,pad=58)
    panel_letter(axb,'b',dx=-0.12,dy=1.22)
    SRC['Fig6b']=pd.DataFrame(rec)
    # c note box
    axc=fig.add_axes([0.60,0.07,0.38,0.27]); axc.set_axis_off()
    nT1=(asum.lead_cs_tier==1).sum(); nT4=(asum.lead_cs_tier==4).sum(); nalt=int(B.best_T1_mlog10p.notna().sum())
    axc.text(0,1,f'Lead CS tier across 42 targets: T1 {nT1}, T2 {(asum.lead_cs_tier==2).sum()}, T3 {(asum.lead_cs_tier==3).sum()}, T4 {nT4}.\n'
             f'{nalt} of {nT4} T4-lead targets also carry a T1 credible set\nthat could serve as an alternative instrument.\n\n'
             'T4 does not mean artefact: functional PAVs (e.g. in PCSK9,\nAPOE, LPL) are graded T4 because affinity data alone cannot\nseparate abundance from epitope effects.',fontsize=6,va='top',ha='left',linespacing=1.3)
    panel_letter(axc,'c',dx=-0.12,dy=1.0)
    SRC['Fig6c']=pd.DataFrame([dict(lead_T1=nT1,lead_T2=(asum.lead_cs_tier==2).sum(),lead_T3=(asum.lead_cs_tier==3).sum(),lead_T4=nT4,T4_lead_with_alternative_T1=nalt)])
    fig.canvas.draw(); save(fig,'Fig6_drug_target_instrument_grades'); return fig




# --- figS2
def figS2grid():
    fig=plt.figure(figsize=(7.2,2.9)); rec=[]
    col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    for gi,(dct,key,nm_) in enumerate([(H3,'primary.grid','MSPA primary label'),(H3,'mspa.grid','MSPA peptide label'),(A3,'grid','Seer full, A1 primary label')]):
        ax=fig.add_axes([0.07+gi*0.315,0.22,0.28,0.60]); gd=pd.DataFrame(json.loads(dct[key])); x=0; xt=[]; xl_=[]
        for cc,cn in [('concordant','Concordant'),('discordant','Discordant'),('na','Not assayed')]:
            for _,r in gd[gd.conc_cat==cc].iterrows():
                k=int(r['sum']); n=int(r['size']); lo_,hi_=wilson(k,n)
                ax.bar(x,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.8); ax.errorbar(x,r['mean'],yerr=[[max(r['mean']-lo_,0)],[max(hi_-r['mean'],0)]],fmt='none',ecolor='#222222',lw=0.6,capsize=1)
                if k==0: ax.plot(x,0.008,marker='_',ms=6,mew=2,color=col[(r.A_eqtl,r.PAV)],clip_on=False)
                ax.text(x,hi_+0.02,f'{k}/{n}',ha='center',va='bottom',fontsize=6,rotation=90)
                rec.append(dict(label=nm_,concordance=cc,A_eqtl=r.A_eqtl,PAV=r.PAV,supported=k,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_)); x+=1
            xt.append(x-2.5); xl_.append(cn); x+=0.8
        ax.set_xticks(xt); ax.set_xticklabels(xl_,fontsize=6); ax.set_ylim(0,1.3); ax.set_yticks([0,0.5,1.0]); ax.set_title(nm_,loc='left',fontsize=7)
        if gi==0: ax.set_ylabel('MS-supported fraction')
        else: ax.set_yticklabels([])
        panel_letter(ax,'abc'[gi],dx=-0.12 if gi else -0.2)
    fig.legend(handles=[Patch(color=col[k],label=n_) for k,n_ in [((0,0),'Neither'),((0,1),'PAV only'),((1,0),'Anchor only'),((1,1),'Anchor + PAV')]],loc='lower center',ncol=4,fontsize=6,handlelength=1,bbox_to_anchor=(0.5,0.0))
    SRC['FigS2']=pd.DataFrame(rec); fig.canvas.draw(); save(fig,'FigS2_MS_support_full_grid'); plt.close(fig)


# --- figS3
def figS3x():
    fig=plt.figure(figsize=(7.2,3.0)); rec=[]
    anc=[('ancestry_afr','African'),('ancestry_csa','Central/South Asian'),('ancestry_eas','East Asian'),('ancestry_mid','Middle Eastern'),('ancestry_amr','Admixed American')]
    col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    for i,(k,l) in enumerate(anc):
        ax=fig.add_axes([0.07+i*0.108,0.20,0.092,0.62]); rt=pd.DataFrame(json.loads(H5[f'{k}.rates']))
        for j,r in rt.iterrows():
            n=int(r['size']); kk=int(round(r['mean']*n)); lo_,hi_=wilson(kk,n)
            ax.bar(j,r['mean'],color=col[(r.A_eqtl,r.PAV)],width=0.8); ax.errorbar(j,r['mean'],yerr=[[r['mean']-lo_],[hi_-r['mean']]],fmt='none',ecolor='#222222',lw=0.6,capsize=1)
            rec.append(dict(ancestry=l,A_eqtl=r.A_eqtl,PAV=r.PAV,k=kk,n=n,rate=r['mean'],wilson_lo=lo_,wilson_hi=hi_))
        ax.set_ylim(0,0.8); ax.set_xticks([]); ax.set_title(l.replace('/','/\n') if '/' in l else l.replace(' ','\n'),loc='left',fontsize=7)
        if i==0: ax.set_ylabel('Replication rate (same sign, P < 0.05)'); panel_letter(ax,'a',dx=-0.75)
        else: ax.set_yticklabels([])
    fig.legend(handles=[Patch(color=col[k],label=n_) for k,n_ in [((0,0),'Neither'),((0,1),'PAV only'),((1,0),'Anchor only'),((1,1),'Anchor + PAV')]],loc='lower left',bbox_to_anchor=(0.06,0.0),ncol=4,fontsize=6,handlelength=1)
    SRC['FigS3a']=pd.DataFrame(rec)
    ph=xl['PostHoc_exploratory']; ph=ph[ph.analysis.str.startswith('Cross-ancestry slope')].copy(); ph['label']=ph.analysis.str.replace('Cross-ancestry slope: ','')
    lbl={'WLS_all(primary secondary)':'All (pre-specified)','OLS_unweighted':'Unweighted OLS','WLS_|eurb|<=1':'|EUR β| ≤ 1','WLS_noPAV':'PAV− only','WLS_PAV':'PAV+ only','WLS_African':'African only','WLS_Central_South Asian':'C/S Asian only','WLS_excl_beta_inconsistent':'excl. 204 β-inconsistent CS'}
    ax=fig.add_axes([0.765,0.20,0.215,0.62]); y=np.arange(len(ph))[::-1]
    for yi,(_,r) in zip(y,ph.iterrows()):
        ax.plot([r.slope_unanchored,r.slope_anchored],[yi,yi],color='#DDDDDD',lw=0.8)
        ax.plot(r.slope_unanchored,yi,'o',ms=4,mfc='white',mec=L_XA,mew=1); ax.plot(r.slope_anchored,yi,'o',ms=4,color=L_XA)
    ax.axvline(1,color='#999999',lw=0.6,ls=(0,(3,2))); ax.set_yticks(y); ax.set_yticklabels([lbl.get(l,l) for l in ph.label]); ax.set_xlabel('non-EUR β on EUR β slope',labelpad=1); ax.set_xticks([0.5,1.0,1.5])
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',mec=L_XA,label='unanchored'),Line2D([],[],marker='o',ls='',color=L_XA,label='anchored')],loc='upper left',bbox_to_anchor=(0.0,-0.20),ncol=2,fontsize=6,borderaxespad=0.1)
    ax.set_title('Slope gap after β-QC (post hoc)',loc='left',fontsize=8,x=-0.6); panel_letter(ax,'b',dx=-0.68)
    SRC['FigS3b']=ph.assign(label=[lbl.get(l,l) for l in ph.label])[['label','slope_unanchored','slope_anchored','P','n']].rename(columns={'P':'P_interaction','n':'n_rows'})
    fig.canvas.draw(); save(fig,'FigS3_cross_ancestry_by_ancestry_slopes'); plt.close(fig)




# --- figS4
def figS4x():
    feats=['A_eqtl','log1p_nanchor','A_iso','PAV','pav_pp','PAVxNoAnchor','conc_concordant','conc_discordant','mlp','absb','maf','logn95','purityMinR2','logdist']
    fl={'A_eqtl':'eQTL anchor','log1p_nanchor':'log(1 + anchor datasets)','A_iso':'isoform anchor','PAV':'PAV in 95% CS','pav_pp':'summed PAV PIP','PAVxNoAnchor':'PAV × no anchor','conc_concordant':'SomaScan concordant','conc_discordant':'SomaScan discordant','mlp':'−log10 P (capped)/10','absb':'|lead β|','maf':'MAF','logn95':'log 95% CS size','purityMinR2':'CS purity (min r²)','logdist':'log distance to TSS'}
    fig=plt.figure(figsize=(7.2,3.8))
    ax=fig.add_axes([0.20,0.22,0.38,0.68]); y=np.arange(len(feats))[::-1]; rec=[]
    for yi,f_ in zip(y,feats):
        a=H4[f'H4_ms.coef_full_std.{f_}']; b=H4[f'H4_xa.coef_full_std.{f_}']
        ax.plot([a,b],[yi,yi],color='#DDDDDD',lw=0.8,zorder=1)
        c_=A4[f'H4_ext.coef_full_std.{f_}']; ax.plot(a,yi,'o',ms=4.5,color=L_MS,zorder=3); ax.plot(c_,yi,'D',ms=3.6,mfc='white',mec=L_MS,mew=1,zorder=3); ax.plot(b,yi,'s',ms=4,color=L_XA,zorder=3)
        rec.append(dict(feature=f_,label=fl[f_],coef_std_MSPA_target=a,coef_std_SeerFull_target=c_,coef_std_crossancestry_target=b))
    ax.axvline(0,color='#999999',lw=0.6); ax.set_yticks(y); ax.set_yticklabels([fl[f] for f in feats])
    ax.set_xlabel('Standardized L2-logistic coefficient (training chromosomes)'); ax.set_xlim(-1.5,1.6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=L_MS,label='MSPA support (n = 99)'),Line2D([],[],marker='D',ls='',mfc='white',mec=L_MS,label='Seer-full support (n = 228)'),Line2D([],[],marker='s',ls='',color=L_XA,label='cross-ancestry (n = 4,572)')],loc='upper left',bbox_to_anchor=(-0.05,-0.16),ncol=3,fontsize=6,borderaxespad=0,title='training target',title_fontsize=6)
    ax.set_title('Signal strength: opposite weights for MS vs cross-ancestry',loc='left',fontsize=8,x=-0.45)
    panel_letter(ax,'a',dx=-0.50)
    SRC['FigS4a']=pd.DataFrame(rec)
    ax=fig.add_axes([0.74,0.22,0.24,0.68]); rows_=[]
    items=[('MSPA target, held-out',H4['H4_ms.auc_full'],L_MS),('Seer-full target, held-out',A4['H4_ext.auc_full'],L_MS),('→ held-out cross-ancestry',H4['transfer_ms_to_xa.auc1'],L_XA),('→ MSPA peptide label',H4['transfer_ms_to_ms_mspa.auc1'],L_MS),('→ Seer discovery overlap',H4['transfer_ms_to_ms_seer.auc1'],L_MS)]
    for i,(l,v,c) in enumerate(items):
        yy=len(items)-1-i; ax.plot([0.5,v],[yy,yy],color='#BBBBBB',lw=0.9); ax.plot(v,yy,'o',ms=5,color=c); ax.text(v+(0.03 if v>0.5 else -0.03),yy,f'{v:.2f}',va='center',ha='left' if v>0.5 else 'right',fontsize=6)
        rows_.append(dict(evaluation=l,AUC=v))
    ax.axvline(0.5,color='#999999',lw=0.6,ls=(0,(3,2))); ax.set_xlim(0,1); ax.set_yticks(range(len(items))[::-1]); ax.set_yticklabels([i[0] for i in items])
    ax.set_xlabel('AUC (0.5 = chance)'); ax.set_title('MS-trained model does not transfer',loc='left',fontsize=8,x=-0.6)
    panel_letter(ax,'b',dx=-0.68)
    SRC['FigS4b']=pd.DataFrame(rows_)
    fig.canvas.draw(); save(fig,'FigS4_grade_model_coefficients_transfer'); plt.close(fig)




# --- figS5
def figS5x():
    fig=plt.figure(figsize=(7.2,2.9)); recs={}
    ax=fig.add_axes([0.095,0.17,0.20,0.68]); tc=R.tier.value_counts().sort_index()
    ax.bar(tc.index,tc.values,color=[TIER_C[t] for t in tc.index],width=0.7)
    for t,v in tc.items(): ax.text(t,v+150,f'{v:,}',ha='center',fontsize=6)
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.set_ylabel('Credible sets'); ax.set_ylim(0,11200)
    ax.set_title('Tiers (n = 14,472)',loc='left',fontsize=8); panel_letter(ax,'a',dx=-0.40)
    recs['FigS5a']=tc.rename('n_cs').reset_index().rename(columns={'index':'tier'})
    ax=fig.add_axes([0.40,0.17,0.25,0.68])
    ct=pd.crosstab(R.platform_concordance,R.tier); ct=ct.loc[['concordant','discordant','na']]; pct=ct.div(ct.sum(1),axis=0)
    left=np.zeros(3)
    for t in [1,2,3,4]:
        ax.barh(range(3),pct[t],left=left,color=TIER_C[t],height=0.65); left+=pct[t]
    ax.set_yticks(range(3)); ax.set_yticklabels([f'{i}\n(n={ct.loc[i].sum():,})' for i in ['concordant','discordant','not assayed']] if False else [f'Concordant\n(n={ct.loc["concordant"].sum():,})',f'Discordant\n(n={ct.loc["discordant"].sum():,})',f'Not assayed\n(n={ct.loc["na"].sum():,})'])
    ax.invert_yaxis(); ax.set_xlim(0,1); ax.set_xlabel('Fraction of credible sets')
    t4c=ct.loc['concordant',4]; ax.text(1-pct.loc['concordant',4]/2,0,f'{t4c}',ha='center',va='center',fontsize=6,color='white')
    ax.set_title('Tier by SomaScan concordance',loc='left',fontsize=8); panel_letter(ax,'b',dx=-0.42)
    recs['FigS5b']=ct.reset_index().merge(pct.add_prefix('frac_T').reset_index(),on='platform_concordance')
    ax=fig.add_axes([0.75,0.17,0.23,0.68])
    data=[R.loc[R.tier==t,'lead_mlog10p'].clip(upper=300).values for t in [1,2,3,4]]
    bp=ax.boxplot(data,positions=[1,2,3,4],widths=0.6,patch_artist=True,showfliers=False,medianprops=dict(color='#222222',lw=1))
    for p_,t in zip(bp['boxes'],[1,2,3,4]): p_.set_facecolor(TIER_C[t]); p_.set_edgecolor('#444444')
    ax.set_yscale('log'); ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.set_ylabel('Lead −log10 P (capped at 300)')
    ax.set_yticks([10,30,100,300]); ax.set_yticklabels(['10','30','100','300']); ax.minorticks_off()
    ax.set_title('Signal strength differs by tier',loc='left',fontsize=8); panel_letter(ax,'c',dx=-0.38)
    recs['FigS5c']=R.groupby('tier').lead_mlog10p.describe().reset_index()
    SRC.update(recs); fig.canvas.draw(); save(fig,'FigS5_resource_overview'); plt.close(fig)




# --- figS6
def figS6():
    fig=plt.figure(figsize=(7.2,2.7))
    ax=fig.add_axes([0.215,0.18,0.20,0.62])
    steps=[('Seer studies in range',1980),('mapped to a gene',1970),('mapped to UKB-PPP cis genes',564),('autosomal windows validated',555),('full read, MD5 verified',218),('early stop, MD5 not verifiable',340)]
    y=np.arange(len(steps))[::-1]
    ax.barh(y,[s[1] for s in steps],color=[L_MS,L_MS,L_MS,L_MS,'#7FBFB8','#CFE5E2'],height=0.65)
    for yi,(l,v) in zip(y,steps): ax.text(v+30,yi,f'{v:,}',va='center',fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([s[0] for s in steps]); ax.set_xlim(0,2350); ax.set_xlabel('Seer GWAS studies (files)')
    ax.set_title('Extraction (addendum A1)',loc='left',fontsize=8); fig.text(0.015,0.88,'a',fontsize=10,fontweight='bold')
    ax=fig.add_axes([0.50,0.20,0.20,0.58])
    st=['supported','refuted','indeterminate']; vals=[167,200,1344]
    ax.bar(range(3),vals,color=[L_MS,'#7A7A7A','#DDDDDD'],width=0.65)
    for i,v in enumerate(vals): ax.text(i,v+25,f'{v:,}',ha='center',fontsize=6)
    ax.set_xticks(range(3)); ax.set_xticklabels(['supported','refuted','indeter-\nminate']); ax.set_ylim(0,1500); ax.set_yticks([0,400,800,1200])
    ax.set_ylabel('UKB-PPP cis CS'); ax.set_title('Seer-full labels (1,711 CS)\nrefuted: power ≥ 0.99',loc='left',fontsize=6.5); fig.text(0.455,0.88,'b',fontsize=10,fontweight='bold')
    ax=fig.add_axes([0.81,0.24,0.15,0.50]); t=np.array(json.loads(A4['agreement_mspa_vs_ext.table']))
    ax.imshow(t,cmap='Greens',vmin=0,vmax=t.max()*1.3)
    for i in range(2):
        for j in range(2): ax.text(j,i,str(t[i,j]),ha='center',va='center',fontsize=7)
    ax.set_xticks([0,1]); ax.set_yticks([0,1])
    lab_=['refuted','supported']
    ax.set_xticklabels(lab_); ax.set_yticklabels(lab_); ax.set_xlabel('Seer full (A1)'); ax.set_ylabel('MSPA primary')
    ax.set_title(f'Label agreement {A4["agreement_mspa_vs_ext.agree"]*100:.1f}%,\nκ = {A4["agreement_mspa_vs_ext.kappa"]:.2f} (n = 148)',loc='left',fontsize=7)
    for s_ in ax.spines.values(): s_.set_visible(False)
    fig.text(0.745,0.88,'c',fontsize=10,fontweight='bold')
    SRC['FigS6a']=pd.DataFrame(steps,columns=['step','n']); SRC['FigS6b']=pd.DataFrame({'label':st,'n_cs':vals})
    SRC['FigS6c']=pd.DataFrame(t,index=['MSPA_row0','MSPA_row1'],columns=['A1_col0','A1_col1']).reset_index()
    SRC['FigS6c']['note']='table as stored in workbook A1_H4_ext_detail agreement_mspa_vs_ext.table; agreement 0.946, kappa 0.892'
    fig.canvas.draw(); save(fig,'FigS6_A1_extraction_and_label_agreement'); plt.close(fig)

# NOTE: after drawing, supplementary files were renumbered by first citation: S6->S2 (A1 QC), S2->S3 (MS grid), S3->S4 (cross-ancestry), S4->S5 (model coefficients), S5->S6 (resource overview).

# ===== v1.2 overrides (addendum A2) =====

def fig3():
    fig=plt.figure(figsize=(7.2,5.6))
    top=0.60; hgt=0.33
    axs=[fig.add_axes([0.075,top,0.25,hgt]),fig.add_axes([0.40,top,0.25,hgt]),fig.add_axes([0.73,top,0.255,hgt])]
    rowsA=[];rowsB=[]
    for i,(h,lab) in enumerate([('H1','Sun 2018'),('H2','deCODE')]):
        zt=json.loads(D[h]['effect_corr.fisher_z_PAV1_lt_PAV0']) if isinstance(D[h]['effect_corr.fisher_z_PAV1_lt_PAV0'],str) else D[h]['effect_corr.fisher_z_PAV1_lt_PAV0']
        for j,(pv,c,gl) in enumerate([(0,C_PAV0,'PAV−'),(1,C_PAV1,'PAV+')]):
            n=D[h][f'effect_corr.PAV{pv}.n']; r=D[h][f'effect_corr.PAV{pv}.r']; rs=D[h][f'effect_corr.PAV{pv}.r_samedir']; ns=D[h][f'effect_corr.PAV{pv}.n_samedir']
            lo_,hi_=fisher_ci(r,n); xx=i+(j-0.5)*0.36
            axs[0].plot([xx,xx],[lo_,hi_],color=c,lw=1.4); axs[0].plot(xx,r,'o',ms=6,color=c,mec='white',mew=0.6,zorder=3)
            axs[0].plot(xx+0.10,rs,'o',ms=4.5,mfc='white',mec=c,mew=1.0,zorder=3)
            rowsA.append(dict(layer=lab,PAV=gl,n_pairs=n,pearson_r=r,fisherCI_lo=lo_,fisherCI_hi=hi_,r_same_direction_only=rs,n_same_direction=ns))
            no=D[h][f'effect_corr.PAV{pv}.n_opp']; fr_=no/n; lo2,hi2=wilson(no,n)
            axs[1].bar(xx,fr_,width=0.34,color=c); axs[1].errorbar(xx,fr_,yerr=[[fr_-lo2],[hi2-fr_]],fmt='none',ecolor='#222222',lw=0.8,capsize=2)
            axs[1].text(xx,hi2+0.006,f'{no}/{n}',ha='center',va='bottom',fontsize=6)
            rowsB.append(dict(layer=lab,PAV=gl,n_opposite=no,n_pairs=n,frac=fr_,wilson_lo=lo2,wilson_hi=hi2))
        axs[0].text(i,0.97,f'P = {fmt_p(zt[1])}',ha='center',fontsize=6)
        po=D[h]['effect_corr.opp_fisher_PAV1_gt.P_one_sided']; axs[1].text(i,0.285,f'P = {fmt_p(po)}',ha='center',fontsize=6)
    axs[0].set_ylim(0.0,1.04); axs[0].set_ylabel('Cross-platform β correlation (r)')
    axs[0].set_title('PAV lowers β agreement',loc='left',fontsize=8)
    axs[0].legend(handles=[Line2D([],[],marker='o',ls='',color=C_PAV0,label='PAV−'),Line2D([],[],marker='o',ls='',color=C_PAV1,label='PAV+'),Line2D([],[],marker='o',ls='',mfc='white',mec='#555555',label='same-direction pairs only')],loc='lower left',fontsize=6,borderaxespad=0.1,handletextpad=0.1)
    axs[1].set_ylim(0,0.32); axs[1].set_ylabel('Fraction with opposite sign')
    axs[1].set_title('…and raises sign reversals',loc='left',fontsize=8)
    for a in axs[:2]: a.set_xticks([0,1]); a.set_xticklabels(['Sun 2018','deCODE']); a.set_xlim(-0.6,1.6)
    lab={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}; col={(0,0):C_NONE,(0,1):C_PAVNA,(1,0):C_ANCH,(1,1):C_ANCHPAV}
    rowsC=[]
    for i,(h,ln) in enumerate([('H1','Sun 2018'),('H2','deCODE')]):
        bb=pd.DataFrame(json.loads(D[h]['effect_corr.by_anchor_PAV']))
        for j,r in bb.iterrows():
            k=int(round(r.opp_frac*r.n)); lo2,hi2=wilson(k,r.n); xx=i*5+j
            axs[2].bar(xx,r.opp_frac,width=0.75,color=col[(r.A_eqtl,r.PAV)]); axs[2].errorbar(xx,r.opp_frac,yerr=[[r.opp_frac-lo2],[hi2-r.opp_frac]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
            axs[2].text(xx,hi2+0.008,f'{k}/{int(r.n)}',ha='center',va='bottom',fontsize=6,rotation=90)
            rowsC.append(dict(layer=ln,A_eqtl=r.A_eqtl,PAV=r.PAV,n_opposite=k,n_pairs=int(r.n),frac=r.opp_frac,wilson_lo=lo2,wilson_hi=hi2))
    axs[2].set_xticks([1.5,6.5]); axs[2].set_xticklabels(['Sun 2018','deCODE']); axs[2].set_ylim(0,0.52)
    axs[2].set_ylabel('Fraction with opposite sign'); axs[2].set_title('Sign reversals by anchor and PAV',loc='left',fontsize=8)
    axs[2].legend(handles=[Patch(color=col[k],label=lab[k]) for k in [(0,0),(0,1),(1,0),(1,1)]],loc='upper right',fontsize=6,ncol=2,borderaxespad=0.1,handlelength=1,columnspacing=0.6)
    for a,l in zip(axs,'abc'): panel_letter(a,l,dx=-0.30)
    # ---- row 2: addendum A2 (independent SomaScan v4 cohorts)
    bot=0.12; h2=0.33
    ax=fig.add_axes([0.20,bot,0.21,h2]); rowsD=[]
    order=['deCODE (H2)','ARIC EA','Fenland','ARIC AA (secondary)','Random-effects meta (deCODE+ARIC EA+Fenland)']
    short={'deCODE (H2)':'deCODE','ARIC EA':'ARIC EA','Fenland':'Fenland','ARIC AA (secondary)':'ARIC AA (sec.)','Random-effects meta (deCODE+ARIC EA+Fenland)':'Meta (3 EA cohorts)'}
    yy=[]; lbl=[]; y0=11.5
    for ep,cc_,mk in [('H2a: eQTL anchor',L_AFF,'o'),('H2b: PAV (power-matched)',C_PAV1,'s')]:
        ax.text(-0.02,y0+0.15,'H2a eQTL anchor' if ep.startswith('H2a') else 'H2b PAV, power-matched',transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=6.5,fontweight='bold',color=cc_)
        for i,co in enumerate(order):
            r=FOREST[(FOREST.endpoint==ep)&(FOREST.cohort==co)].iloc[0]; y=y0-1-i
            meta=co.startswith('Random'); sec='AA' in co
            ax.plot([r.lo,r.hi],[y,y],color=cc_,lw=1.2)
            ax.plot(r.OR,y,'D' if meta else mk,ms=5.5 if meta else 4.5,color=cc_ if not sec else 'white',mec=cc_,mew=1.0,zorder=3)
            ax.text(1.02,y,f'{r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f})',transform=ax.get_yaxis_transform(),va='center',fontsize=6)
            yy.append(y); lbl.append(short[co]); rowsD.append(dict(endpoint=ep,cohort=co,OR=r.OR,lo=r.lo,hi=r.hi))
        y0-=6.5
    ax.set_yticks(yy); ax.set_yticklabels(lbl); ax.set_ylim(-1.6,12.2); ax.set_xscale('log'); ax.set_xlim(1.5,5)
    ax.set_xticks([2,3,4]); ax.set_xticklabels(['2','3','4']); ax.minorticks_off(); ax.set_xlabel('OR for replication (95% CI)')
    ax.set_title('Reproduced in independent SomaScan v4 cohorts',loc='left',fontsize=8,x=-0.70)
    ax.text(1.02,-1.25,'meta I²: H2a 0%, H2b 21%',transform=ax.get_yaxis_transform(),fontsize=6,color='#444444')
    panel_letter(ax,'d',dx=-0.78)
    # e/f: effect agreement by cohort
    coh=[('deCODE',D['H2'],'effect_corr'),('ARIC EA',A2D['A2_aric_ea_detail'],'effect_corr'),('Fenland',A2D['A2_fenland_detail'],'effect_corr'),('ARIC AA',A2D['A2_aric_aa_detail'],'effect_corr')]
    axe=fig.add_axes([0.625,bot,0.155,h2]); axf=fig.add_axes([0.845,bot,0.14,h2]); rowsE=[]
    for i,(cn,dd,pre) in enumerate(coh):
        for j,(pv,c) in enumerate([(0,C_PAV0),(1,C_PAV1)]):
            n=int(dd[f'{pre}.PAV{pv}.n']); r=float(dd[f'{pre}.PAV{pv}.r']); no=int(dd[f'{pre}.PAV{pv}.n_opp'])
            lo_,hi_=fisher_ci(r,n); xx=i+(j-0.5)*0.36
            axe.plot([xx,xx],[lo_,hi_],color=c,lw=1.2); axe.plot(xx,r,'o',ms=4.5,color=c,mec='white',mew=0.5,zorder=3)
            fr_=no/n; lo2,hi2=wilson(no,n)
            axf.bar(xx,fr_,width=0.34,color=c); axf.errorbar(xx,fr_,yerr=[[fr_-lo2],[hi2-fr_]],fmt='none',ecolor='#222222',lw=0.7,capsize=1.5)
            rowsE.append(dict(cohort=cn,PAV=pv,n_pairs=n,pearson_r=r,fisherCI_lo=lo_,fisherCI_hi=hi_,n_opposite=no,frac_opposite=fr_,wilson_lo=lo2,wilson_hi=hi2))
    for a in (axe,axf):
        a.set_xticks(range(4)); a.set_xticklabels([c[0] for c in coh],rotation=90); a.set_xlim(-0.6,3.6)
    axe.set_ylim(0.0,1.0); axe.set_ylabel('β correlation (r)'); axf.set_ylim(0,0.32); axf.set_ylabel('Fraction opposite sign')
    axe.set_title('β agreement',loc='left',fontsize=8); axf.set_title('Sign reversals',loc='left',fontsize=8)
    axe.legend(handles=[Line2D([],[],marker='o',ls='',color=C_PAV0,label='PAV−'),Line2D([],[],marker='o',ls='',color=C_PAV1,label='PAV+')],loc='lower left',fontsize=6,borderaxespad=0.1,handletextpad=0.1)
    panel_letter(axe,'e',dx=-0.50); panel_letter(axf,'f',dx=-0.52)
    SRC['Fig3a']=pd.DataFrame(rowsA); SRC['Fig3b']=pd.DataFrame(rowsB); SRC['Fig3c']=pd.DataFrame(rowsC); SRC['Fig3d']=pd.DataFrame(rowsD); SRC['Fig3ef']=pd.DataFrame(rowsE)
    SRC['Fig3a']['note']='pairs = replicated CS with harmonisable lead/shared variant; Fisher z one-sided P (PAV+ < PAV−): Sun 2.4e-06, deCODE 0.0055'
    SRC['Fig3d']['note']='addendum A2 (workbook v4 sheets A2_summary, A2_meta_SomaScan_v4); meta = DerSimonian–Laird; ARIC AA secondary'
    fig.canvas.draw(); save(fig,'Fig3_effect_agreement_sign_reversal'); return fig

# Fig. 1b endpoint rows (v1.2): built from Additional file 2 sheet Fig1b (in earlier versions this list was a session variable)
def build_fig1_rows(af2_path):
    s=pd.read_excel(af2_path,sheet_name='Fig1b'); out=[]
    for _,r in s.iterrows():
        col=L_AFF if r.hypothesis.startswith(('H2a','H2b')) else (L_XA if r.hypothesis=='H5' else L_MS)
        out.append((r.hypothesis,r.endpoint,r.OR,r.lo95,r.hi95,r.P_holm,float(r.threshold),None,col,r.n))
    return out
# A2 inputs (v1.2): FOREST = FigS_A2_forest_source.csv; A2D = workbook v4 sheets A2_aric_ea_detail / A2_aric_aa_detail / A2_fenland_detail as key->value dicts

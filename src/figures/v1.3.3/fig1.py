from common import *
def make(WB):
    cs=WB['01_Confirmatory_summary'].set_index('hypothesis'); a2=WB['A2_summary'].set_index('hypothesis'); a1=WB['A1_Addendum_summary'].set_index('hypothesis')
    F=Fig(170,178,'Fig1')
    ax=F.panel('a','Each UKB-PPP cis credible set is tested against four kinds of evidence',0,0,170,74,ml=0,mr=0,mt=6,mb=0,row=0)
    ax.set_xlim(0,170); ax.set_ylim(0,68); ax.axis('off')
    for x,t in [(1,'Discovery'),(33,'Features of each credible set'),(73,'Evaluation layers'),(132,'Annotation resource')]:
        ax.text(x,66,t,fontsize=FS,color=C['sub'],style='italic',va='top')
    box(ax,1,2,29,59,'UKB-PPP Olink','European ancestry\n14,472 cis credible\nsets (CS) in\n1,879 genes\n\nOpen Targets 26.09\nSuSiE-inf 95% sets',C['text'])
    fh=18.3
    box(ax,33,61-fh,36,fh,'eQTL anchor','colocalisation with a\ncognate eQTL (GTEx,\neQTL Catalogue)',C['anchor'])
    box(ax,33,61-2*fh-2,36,fh,'Protein-altering variant','missense, LoF or splice\nvariant in the 95% CS',C['pav'])
    box(ax,33,2,36,fh,'Negative control','colocalisation with a\nnon-cognate eQTL',C['ctrl'])
    layers=[('SomaScan v3: Sun 2018','4,153 CS; internal reproduction',C['conc']),('SomaScan v4: deCODE, ARIC, Fenland','6,711, 7,783 and 7,296 CS',C['conc']),
            ('MS, Seer study: MSPA, Seer','169 and 367 evaluable CS (shared samples)',MSL['MSPA']),('MS, independent: Gen. Scotland, HOLBAEK','450 and 403 evaluable CS',MSL['Seer']),
            ('Olink, non-European ancestries','8,065 evaluable CS',C['xa']),('Published drug-target MR','20,027 protein–outcome pairs, 11 studies',C['mr'])]
    lh=8.6; gap=1.5
    for k,(t,b,c) in enumerate(layers):
        y0=61-(k+1)*lh-k*gap
        ax.add_patch(FancyBboxPatch((73,y0),56,lh,boxstyle='round,pad=0,rounding_size=1.2',lw=0.8,ec=c,fc='white'))
        ax.text(74.5,y0+lh-1.3,t,fontsize=FS,fontweight='bold',color=c,va='top'); ax.text(74.5,y0+1.3,b,fontsize=FS,color=C['text'],va='bottom')
    ax.add_patch(FancyBboxPatch((132,2),37,59,boxstyle='round,pad=0,rounding_size=1.2',lw=0.8,ec=C['text'],fc='white'))
    ax.text(133.5,59.6,'PAV-aware annotations',fontsize=FS,fontweight='bold',va='top')
    ax.text(133.5,55.0,'Tier 1–4 for all 14,472 CS;\nevidence patterns, not a\nreliability ranking',fontsize=FS,va='top',linespacing=1.25)
    tiers=[(1,'anchored, no PAV'),(2,'anchored + PAV, or\nisoform-anchored only'),(3,'no anchor, no PAV'),(4,'no anchor + PAV')]
    yy=42.5
    for t,lab in tiers:
        ax.add_patch(Rectangle((133.5,yy-2.6),2.6,2.6,fc=TIER[t],ec='none'))
        ax.text(137.3,yy,f'T{t}  {lab}',fontsize=FS,va='top',linespacing=1.25); yy-=(6.4 if '\n' in lab else 3.6)
    ax.text(133.5,yy-1.6,'Flag: concordant, unanchored,\nPAV-containing (339 CS)',fontsize=FS,va='top',linespacing=1.25)
    for xa_,xb_ in [(30,33),(69,73),(129,132)]:
        ax.annotate('',xy=(xb_,31.5),xytext=(xa_,31.5),arrowprops=dict(arrowstyle='-|>',lw=0.8,color='#666666',shrinkA=0,shrinkB=0))
    def r_cs(h,lab,col,thr,mk='o'):
        c=cs.loc[h]; return dict(label=lab,est=c.estimate,lo=c.lo,hi=c.hi,P=c.P_holm,thr=thr,color=col,marker=mk,v=c.verdict,source=f'01_Confirmatory_summary:{h}')
    def r_a2(h,lab,col,thr,mk):
        c=a2.loc[h]; return dict(label=lab,est=c.OR,lo=c.lo,hi=c.hi,P=c.P_holm_A2,thr=thr,color=col,marker=mk,v=c.verdict,source=f'A2_summary:{h}')
    def wrow(lab,k,P,v,src):
        o,l,h=woolf(*k); return dict(label=lab,est=o,lo=l,hi=h,P=P,thr=0.5,color=C['pav'],marker='D',ms=2.8,v=v,source=src+' (Woolf CI from counts '+'/'.join(map(str,k))+')')
    rows=[dict(header=True,label='Replication on SomaScan',color=C['conc']),
          r_cs('H2a','eQTL anchor → replication, deCODE',C['anchor'],1.5),r_a2('H2a-ARIC','eQTL anchor → replication, ARIC',C['anchor'],1.5,'o'),r_a2('H2a-Fenland','eQTL anchor → replication, Fenland',C['anchor'],1.5,'o'),
          r_cs('H2b','PAV → replication, deCODE (matched)',C['pav'],1.3,'s'),r_a2('H2b-ARIC','PAV → replication, ARIC (matched)',C['pav'],1.3,'s'),r_a2('H2b-Fenland','PAV → replication, Fenland (matched)',C['pav'],1.3,'s'),
          dict(header=True,label='Replication in non-European ancestries',color=C['xa']),
          r_cs('H5','eQTL anchor → replication, pooled',C['anchor'],1.3),
          dict(header=True,label='Mass-spectrometry support',color=C['ms']),
          r_cs('H3i','eQTL anchor → MS support, MSPA',C['anchor'],1.5),
          dict(label='eQTL anchor → MS support, Seer',est=a1.loc['H3i-ext','estimate'],lo=a1.loc['H3i-ext','lo'],hi=a1.loc['H3i-ext','hi'],P=a1.loc['H3i-ext','P_holm_addendum'],thr=1.5,color=C['anchor'],marker='o',v=a1.loc['H3i-ext','verdict'],source='A1_Addendum_summary:H3i-ext'),
          r_cs('H3i-A4','eQTL anchor → MS support, Generation Scotland',C['anchor'],1.5),r_cs('H3i-HB','eQTL anchor → MS support, HOLBAEK',C['anchor'],1.5),
          wrow('Unanchored PAV vs anchored*, MSPA',(2,8,42,39),cs.loc['H3ii','P_holm'],cs.loc['H3ii','verdict'],'01_Confirmatory_summary:H3ii'),
          wrow('Unanchored PAV vs anchored*, Seer',(5,6,57,68),a1.loc['H3ii-ext','P_holm_addendum'],a1.loc['H3ii-ext','verdict'],'A1_Addendum_summary:H3ii-ext'),
          wrow('Unanchored PAV vs anchored*, Gen. Scotland',(17,5,43,14),cs.loc['H3ii-A4','P_holm'],cs.loc['H3ii-A4','verdict'],'01_Confirmatory_summary:H3ii-A4'),
          dict(header=True,label='Published MR instruments',color=C['mr']),
          r_cs('H-A3b (A3 P2)','Tier 4 instrument → significant result',C['pav'],1.0,'s')]
    for r in rows:
        if not r.get('header'):
            v=str(r['v']).lower(); r['v']='partial' if ('部分' in v or 'partial' in v) else ('not supported' if ('不成立' in v or 'not' in v) else 'supported')
    F.rows_b=rows
    axb=F.panel('b','PAVs and eQTL anchoring predict SomaScan replication; anchoring does not predict MS support',0,77,170,101,ml=63,mr=59,mt=7,mb=9,row=1)
    ys=forest(axb,rows)
    axb.set_xlim(0.04,6); axb.set_xticks([0.1,0.25,0.5,1,2,4]); axb.set_xticklabels(['0.1','0.25','0.5','1','2','4'])
    axb.set_xlabel('Odds ratio (95% CI), log scale')
    axb.plot(0.055,ys[0],marker='|',ms=6,color='#333333',mew=0.8); axb.text(0.068,ys[0],'pre-specified threshold',fontsize=FS,va='center',color='#333333')
    
    tx=F.ax_mm(114,77+7,56,101-7-9); tx.axis('off'); tx.set_ylim(axb.get_ylim()); tx.set_xlim(0,1)
    for xx,t in [(0.0,'OR (95% CI)'),(0.42,'P (Holm)'),(0.70,'Verdict')]: tx.text(xx,len(rows)-0.35,t,fontsize=FS,color=C['sub'],va='bottom')
    tx.text(0.0,ys[9],'* concordant CS only',fontsize=FS,va='center',color=C['sub'])
    for y,r in zip(ys,rows):
        if r.get('header'): continue
        tx.text(0.0,y,f"{r['est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})",fontsize=FS,va='center')
        tx.text(0.42,y,fmt_p(r['P'],'').replace('= ','').strip(),fontsize=FS,va='center')
        tx.text(0.70,y,r['v'],fontsize=FS,va='center')
    rep=F.verify(); F.save(OUT+'/Fig1')
    src=pd.DataFrame([dict(panel='Fig1b',row=r['label'],OR=r['est'],lo95=r['lo'],hi95=r['hi'],P_holm=r['P'],threshold=r['thr'],verdict=r['v'],source=r['source']) for r in rows if not r.get('header')])
    return F,rep,src

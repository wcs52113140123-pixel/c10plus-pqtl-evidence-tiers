from common import *
from fig1_icons import *
TCOUNT={1:2539,2:682,3:10046,4:1205}
def panel_a(F):
    ax=F.panel('a','A protein-altering variant can make both affinity platforms agree for the wrong reason',0,0,170,52,ml=0,mr=0,mt=6,mb=0,row=0)
    ax.set_xlim(0,170); ax.set_ylim(0,46); ax.axis('off')
    cards=[(0,'#EEF3F9',C['anchor'],'Abundance effect','regulatory\nvariant','more protein',('up','up','up'),'All readouts agree'),
           (87.5,'#FCEFE6',C['pav'],'Binding artefact','missense\nvariant (PAV)','same amount,\naltered epitope',('down','down','eq'),'Affinity platforms agree; MS need not')]
    for x0,fc,ec,head_,vlab,plab,res,foot in cards:
        card(ax,x0,0.5,82.5,45,fc)
        ax.text(x0+3,43.2,head_,fontsize=FS_T,fontweight='bold',color=ec,va='top')
        cy=24.5
        dna(ax,x0+4,cy,17,h=4.4,mark=True,mark_c=ec,mark_x=0.55)
        ax.text(x0+12.5,cy-5.2,vlab,fontsize=FS,ha='center',va='top',color=C['text'],linespacing=1.15)
        arrow(ax,x0+22.8,cy,x0+27.6,cy)
        if res[0]=='up':
            pts=[(32.5,31),(38.5,33),(44.8,30.5),(35.5,24.3),(41.8,24),(32.8,17.7),(39,17.4),(45.2,18.2)]; epi='normal'
        else:
            pts=[(33.5,29.5),(41,30.5),(37.3,21.5)]; epi='altered'
        for (px,py) in pts: protein(ax,x0+px,py,r=2.25,c='#8FA9C9',epi=epi,epi_c=C['pav'])
        ax.text(x0+39,12.3,plab,fontsize=FS,ha='center',va='top',color=C['text'],linespacing=1.15)
        arrow(ax,x0+50.2,cy,x0+55,cy)
        rows=[('Olink',antibody,'#3a3a3a'),('SomaScan',aptamer,C['conc']),('MS',None,MSL['Seer'])]
        for k,(nm,fn,col) in enumerate(rows):
            yy=34.5-k*10
            if nm=='MS': spectrum(ax,x0+56.3,yy-2.3,w=5.8,h=4.6,c=col)
            else: fn(ax,x0+59.2,yy-2.5,s=1.05,c=col)
            ax.text(x0+63.6,yy,nm,fontsize=FS,va='center',color=C['text'])
            badge(ax,x0+78.6,yy,res[k],r=1.6)
        ax.text(x0+41.25,4.3,foot,fontsize=FS,ha='center',va='center',color=ec,style='italic',fontweight='bold')
def panel_b(F,y,h):
    ax=F.panel('b','Each UKB-PPP cis credible set is annotated and tested against four independent kinds of evidence',0,y,170,h,ml=0,mr=0,mt=6,mb=0,row=1)
    H=h-6; ax.set_xlim(0,170); ax.set_ylim(0,H); ax.axis('off')
    W=[34,34,62,32]; g=(170-sum(W))/3; xs=[0]; 
    for w in W[:-1]: xs.append(xs[-1]+w+g)
    names=['1  Discover','2  Annotate','3  Test','4  Release']
    for x,w,n in zip(xs,W,names):
        card(ax,x,1,w,H-2,'#F3F4F6')
        ax.text(x+2.5,H-3.2,n,fontsize=FS,color=C['sub'],fontweight='bold',va='top')
    for k in range(3):
        xa=xs[k]+W[k]+0.3; arrow(ax,xa,H/2,xa+g-0.6,H/2,c='#9a9a9a',lw=0.9,hw=1.4)
    # 1 discover
    x=xs[0]
    for k in range(3): person(ax,x+6+k*3.6,H-15.5,s=0.9,c='#7d7d7d')
    tube(ax,x+20,H-16,s=1.0,c='#E9C46A')
    ax.text(x+17,H-21.5,'14,472',fontsize=11,fontweight='bold',ha='center',va='top',color=C['text'])
    ax.text(x+17,H-26.6,'cis credible sets',fontsize=FS,ha='center',va='top')
    ax.text(x+17,H-30.4,'UKB-PPP Olink · 1,879 genes',fontsize=6.2,ha='center',va='top',color=C['sub'])
    # 2 annotate
    x=xs[1]
    for k,(lab,col) in enumerate([('PAV in 95% CS',C['pav']),('eQTL anchor',C['anchor']),('control: wrong-gene eQTL',C['ctrl'])]):
        yy=H-13.5-k*8.6
        card(ax,x+3,yy-2.9,28,5.8,col if k<2 else 'white',ec=col,lw=0.8,r=2.9)
        ax.text(x+17,yy,lab,fontsize=FS if k<2 else 6.3,ha='center',va='center',color='white' if k<2 else col,fontweight='bold' if k<2 else 'normal')
    # 3 test: 2x2 tiles
    x=xs[2]; tw=28.5; th=14.6
    tiles=[('3','SomaScan\ncohorts','apt',C['conc']),('4','MS\ndatasets','ms',MSL['Seer']),('5','non-European\nancestries','globe',C['xa']),('20,027','published MR\npairs','doc',C['mr'])]
    for k,(num,lab,ic,col) in enumerate(tiles):
        cx=x+2+(k%2)*(tw+1.5); cy=H-8-(k//2)*(th+1.6)-th
        card(ax,cx,cy,tw,th,'white')
        if ic=='apt': aptamer(ax,cx+4.5,cy+th/2-2.4,s=1.0,c=col)
        elif ic=='ms': spectrum(ax,cx+1.8,cy+th/2-2.3,w=5.4,h=4.6,c=col)
        elif ic=='globe': globe(ax,cx+4.6,cy+th/2,r=2.4,c=col)
        else: doc(ax,cx+2.8,cy+th/2-2.3,s=1.0,c=col)
        ax.text(cx+9.5,cy+th/2+0.4,num,fontsize=10,fontweight='bold',va='bottom',color=col)
        ax.text(cx+9.5,cy+th/2-0.2,lab,fontsize=6.4,va='top',color=C['text'],linespacing=1.1)
    # 4 release
    x=xs[3]; tot=sum(TCOUNT.values()); bx=x+3; bw=26; by=H-15.5; acc=bx
    for t in (1,2,3,4):
        ww=bw*TCOUNT[t]/tot; ax.add_patch(Rectangle((acc,by),ww,4.2,fc=TIER[t],ec='white',lw=0.4,zorder=5)); acc+=ww
    ax.text(x+16,by-1.6,'Tier 1–4 annotation',fontsize=FS,ha='center',va='top')
    ax.text(x+16,by-5.2,'for every credible set',fontsize=6.2,ha='center',va='top',color=C['sub'])
    magnifier(ax,x+8.5,H-31.5,s=1.0)
    ax.text(x+13,H-31.5,'online\nportal',fontsize=FS,va='center',linespacing=1.1)
def make(WB):
    cs=WB['01_Confirmatory_summary'].set_index('hypothesis'); a2=WB['A2_summary'].set_index('hypothesis'); a1=WB['A1_Addendum_summary'].set_index('hypothesis')
    global Y_C,H_C
    Y_C=104; H_C=101
    F=Fig(170,Y_C+H_C,'Fig1')
    panel_a(F); panel_b(F,55,46)
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
    axb=F.panel('c','PAVs and eQTL anchoring predict SomaScan replication; anchoring does not predict MS support',0,Y_C,170,H_C,ml=63,mr=59,mt=7,mb=9,row=2)
    ys=forest(axb,rows)
    axb.set_xlim(0.04,6); axb.set_xticks([0.1,0.25,0.5,1,2,4]); axb.set_xticklabels(['0.1','0.25','0.5','1','2','4'])
    axb.set_xlabel('Odds ratio (95% CI), log scale')
    axb.plot(0.055,ys[0],marker='|',ms=6,color='#333333',mew=0.8); axb.text(0.068,ys[0],'pre-specified threshold',fontsize=FS,va='center',color='#333333')
    
    tx=F.ax_mm(114,Y_C+7,56,H_C-7-9); tx.axis('off'); tx.set_ylim(axb.get_ylim()); tx.set_xlim(0,1)
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

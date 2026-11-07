from common import *
from fig1_icons import *
from fig1_icons2 import *
TCOUNT={1:2539,2:682,3:10046,4:1205}
def mini_plot(F,ax,H,top_offset,x0,y0,w,h,trend,col):
    """genotype mini-plot placed at panel-axes data box (x0,y0,w,h) in mm"""
    ia=F.ax_mm(x0,top_offset+(H-(y0+h)),w,h)
    means={'up':[0.28,0.5,0.72],'down':[0.72,0.5,0.28],'flat':[0.5,0.5,0.5]}[trend]
    off=[(-0.16,-0.10),(0.10,-0.04),(-0.06,0.06),(0.15,0.11)]
    for g,m in enumerate(means):
        for dx,dy in off: ia.plot(g+dx,m+dy,'o',ms=1.5,mfc=col,mec='none',zorder=3)
        ia.plot([g-0.25,g+0.25],[m,m],color=col,lw=0.9,zorder=4,solid_capstyle='butt')
    ia.plot([0,1,2],means,color=col,lw=0.7,zorder=2)
    ia.set_xlim(-0.55,2.55); ia.set_ylim(0.05,0.98); ia.set_yticks([])
    ia.set_xticks([0,1,2]); ia.set_xticklabels(['0','1','2'],fontsize=5.3); ia.tick_params(axis='x',length=1.4,width=0.4,pad=0.8,colors='#555555')
    for sp in ('left','bottom'): ia.spines[sp].set_linewidth(0.45); ia.spines[sp].set_color('#8C8C8C')
    ia.patch.set_alpha(0)
    return ia
def panel_a(F):
    PH=56; MT=5.5; H=PH-MT
    ax=F.panel('a','A protein-altering variant can make both affinity platforms agree for the wrong reason',0,0,170,PH,ml=0,mr=0,mt=MT,mb=0,row=0)
    ax.set_xlim(0,170); ax.set_ylim(0,H); ax.axis('off')
    cards=[dict(x0=0,fc='#EEF3F9',ec=C['anchor'],title='Abundance effect',vx=4.4,vin=None,vlab='regulatory\nvariant',
                carrier=6,epi='normal',plab='more protein',trends=('up','up','up'),cols=(C['anchor'],)*3,foot='All readouts agree',art=False),
           dict(x0=87.5,fc='#FCEFE6',ec=C['pav'],title='Binding artefact',vx=13.6,vin=1,vlab='missense variant\n(PAV)',
                carrier=3,epi='altered',plab='same amount,\naltered epitope',trends=('down','down','flat'),cols=(C['pav'],C['pav'],C['pavneg']),foot='Affinity platforms agree; MS need not',art=True)]
    for c in cards:
        x0=c['x0']; cy=27.5
        card(ax,x0+0.3,0.4,81.9,H-0.8,c['fc'],r=2.4)
        ax.text(x0+3,H-2.6,c['title'],fontsize=7.6,fontweight='bold',color=c['ec'],va='top')
        # 1 genome
        gene_model(ax,x0+2.6,x0+20.6,cy+2.4,var_x=x0+c['vx'],var_c=c['ec'],exons=(x0+9.2,x0+13.6,x0+18.0),var_in_exon=c['vin'],tss_x=x0+6.9)
        helix(ax,x0+3.0,cy-3.3,17.4,h=2.6)
        ax.text(x0+11.6,cy-6.6,c['vlab'],fontsize=6.6,ha='center',va='top',linespacing=1.1)
        arrow(ax,x0+21.4,cy,x0+24.0,cy)
        # 2 plasma
        ccx,ccy,R=x0+36.6,cy-0.5,9.3
        tube(ax,x0+26.6,cy+5.6,s=0.85)
        ax.plot([x0+27.4,ccx-R*0.62],[cy+9.2,ccy+R*0.79],color=GREY,lw=0.5,ls=(0,(1.6,1.2)),zorder=5)
        ax.add_patch(Circle((ccx,ccy),R,fc=CREAM,ec=GREY,lw=0.55,zorder=4))
        ax.plot([ccx,ccx],[ccy-R+0.4,ccy+R-0.4],color=GREY,lw=0.5,ls=(0,(1.6,1.2)),zorder=5)
        for bx,by in [(-2.3,2.6),(-2.0,-4.4),(2.4,6.7),(4.3,-7.0) if c['carrier']==3 else (1.9,-6.9)]: bgprotein(ax,ccx+bx,ccy+by,r=0.75)
        left=[(-5.4,4.3),(-4.6,-0.9),(-5.0,-5.8)]
        right6=[(2.6,4.3),(6.3,2.9),(2.5,-0.6),(6.4,-1.9),(1.9,-5.9),(5.6,-5.0)][:6]
        right3=[(5.0,4.3),(4.4,-0.9),(5.2,-5.6)]
        for dx,dy in left: protein(ax,ccx+dx,ccy+dy,r=1.45,epi='normal')
        for dx,dy in (right6 if c['carrier']==6 else right3): protein(ax,ccx+dx,ccy+dy,r=1.45,epi=c['epi'])
        ax.text(ccx-4.6,ccy-R-1.0,'non-carrier',fontsize=5.6,color='#6b6b6b',ha='center',va='top')
        ax.text(ccx+4.6,ccy-R-1.0,'carrier',fontsize=5.6,color='#6b6b6b',ha='center',va='top')
        ax.text(ccx,ccy-R-3.9,c['plab'],fontsize=6.6,ha='center',va='top',linespacing=1.1,zorder=9)
        arrow(ax,x0+46.6,cy,x0+49.2,cy)
        # 3 read-outs
        tx=x0+49.9; tw=31.0; th=11.4; gap=1.1; tops=[H-5.6-k*(th+gap) for k in range(3)]
        ax.text(tx+tw-0.6,H-2.6,'level by allele copies',fontsize=5.4,color='#6b6b6b',ha='right',va='top')
        names=['Olink','SomaScan','MS']
        for k in range(3):
            ty=tops[k]-th
            card(ax,tx,ty,tw,th,'white',ec=LG,lw=0.45,r=1.0,z=3)
            ax.text(tx+10.9,ty+th/2,names[k],fontsize=6.1,fontweight='bold',va='center',color=C['text'],zorder=8)
            mini_plot(F,ax,H,MT,tx+23.3,ty+2.6,7.4,th-3.6,c['trends'][k],c['cols'][k])
            if k==0:  # Olink PEA
                px,py=tx+5.5,ty+2.7; top=protein(ax,px,py,r=1.35,epi=c['epi'])
                if not c['art']:
                    s1,_=igg(ax,(px-1.25,py+0.9),128,s=0.82); s2,_=igg(ax,(px+1.25,py+0.9),52,s=0.82)
                    oligo_pair(ax,s1,s2,(px,ty+7.9))
                else:
                    s2,_=igg(ax,(px+1.25,py+0.9),52,s=0.82); s1,_=igg(ax,(px-2.35,py+2.0),128,s=0.82)
                    cross(ax,px-1.45,py+1.25,s=0.42); dangling(ax,s1,-1); dangling(ax,s2,1)
            elif k==1:  # SomaScan aptamer
                px,py=tx+5.5,ty+2.5; top=protein(ax,px,py,r=1.35,epi=c['epi'])
                if not c['art']: aptamer(ax,px,top+0.45,s=0.82)
                else: aptamer(ax,px,top+1.75,s=0.82); cross(ax,px+0.05,top+1.05,s=0.4)
            else:  # MS
                py=ty+7.3; protein(ax,tx+2.2,py,r=1.1,epi=c['epi'])
                scissors(ax,tx+4.4,py,s=0.75)
                for m,xx in enumerate((tx+6.2,tx+7.7,tx+9.2)): peptide(ax,xx,py,w=1.25,h=0.6,dot=(c['art'] and m==1))
                spectrum(ax,tx+5.4,ty+1.6,w=4.5,h=2.6)
        ax.text(x0+41.25,2.7,c['foot'],fontsize=6.8,ha='center',va='center',color=c['ec'],style='italic',fontweight='bold')
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
    tube(ax,x+20,H-16.5,s=1.0)
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
    Y_C=106; H_C=101
    F=Fig(170,Y_C+H_C,'Fig1')
    panel_a(F); panel_b(F,58,46)
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

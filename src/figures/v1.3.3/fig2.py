from common import *
def make(R,WB,S24):
    v=lambda sh,key: float(WB[sh].set_index('key').loc[key,'value'])
    out=[]
    F=Fig(170,132,'Fig2')
    dec=R[(R.somascan_v4_decode_assayed==1)&(R.somascan_v3_sun2018_assayed!=1)]
    sun=R[R.somascan_v3_sun2018_assayed==1]
    def fourgrp(d,rep):
        g=d.groupby(['anchor_eqtl_gtex','pav'])[rep].agg(['sum','count']).reset_index()
        nm={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}
        g['group']=[nm[(int(a),int(b))] for a,b in zip(g.anchor_eqtl_gtex,g.pav)]; return g.rename(columns={'sum':'k','count':'n'})
    gd=fourgrp(dec,'rep_decode_v4'); gs=fourgrp(sun,'rep_sun2018')
    ax=F.panel('a','deCODE: replication by PAV and anchor',0,0,54,62,ml=10,mr=2,row=0); grp_bars(ax,gd,ymax=58)
    ax.text(0.02,0.98,f'{len(dec):,} CS, {dec.gene_id.nunique()} genes',transform=ax.transAxes,fontsize=FS,va='top',color=C['sub'])
    ax=F.panel('b','Sun 2018: the same ordering',56,0,52,62,ml=8,mr=2,row=0); grp_bars(ax,gs,ylab='',ymax=58)
    ax.text(0.02,0.98,f'{len(sun):,} CS, {sun.gene_id.nunique()} genes',transform=ax.transAxes,fontsize=FS,va='top',color=C['sub'])
    for nm,g in [('Fig2a deCODE',gd),('Fig2b Sun 2018',gs)]:
        for _,r in g.iterrows(): lo,hi=wilson(r.k,r.n); out.append(dict(panel=nm.split()[0],layer=nm.split(' ',1)[1],group=r.group,replicated=int(r.k),n=int(r.n),rate=r.k/r.n,wilson_lo=lo,wilson_hi=hi))
    # c SMD
    sm=S24[S24.panel_sheet=='Fig2c']
    ax=F.panel('c','Matching balances signal strength',110,0,60,62,ml=21,mr=2,row=0)
    ys=np.arange(len(sm))[::-1]
    ax.axvspan(-0.1,0.1,color='#EEEEEE',zorder=0,lw=0)
    for y,(_,r) in zip(ys,sm.iterrows()):
        ax.plot([r.SMD_before,r.SMD_after],[y,y],color='#BBBBBB',lw=0.7,zorder=1)
        ax.plot(r.SMD_before,y,'o',mfc='white',mec=C['pav'],ms=4,mew=0.9,zorder=2); ax.plot(r.SMD_after,y,'o',color=C['pav'],ms=4,zorder=3)
        out.append(dict(panel='Fig2c',covariate=r.label,SMD_before=r.SMD_before,SMD_after=r.SMD_after,n_PAV_cs_total=r.n_PAV_cs_total,n_matched_sets=r.n_matched_sets))
    ax.axvline(0,color=C['ref'],lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels(sm.label.str.replace('−log10 P (capped)','Capped −log$_{10}$P').str.replace('log 95% CS size','log CS size')); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.set_xlim(-0.35,0.45); ax.set_xlabel('Standardised mean difference\n(PAV+ minus PAV−)'); ax.set_ylim(-0.6,len(sm)-0.2)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',mec=C['pav'],ms=4,label='before'),Line2D([],[],marker='o',ls='',color=C['pav'],ms=4,label='after matching')],loc='lower right',fontsize=FS,handletextpad=0.2,borderaxespad=0)
    ax.set_ylim(-1.6,len(sm)-0.4)
    # d matched
    md=S24[S24.panel_sheet=='Fig2d']
    ax=F.panel('d','The PAV excess persists after matching',0,67,62,65,ml=10,mr=2,mb=12,row=1)
    Ps={'deCODE':(v('H2_detail','pav_matched.PAV.OR'),cs_holm(WB,'H2b'),'Holm P'),'Sun 2018':(v('H1_detail','pav_matched.PAV.OR'),v('H1_detail','pav_matched.PAV.P'),'P')}
    for j,lay in enumerate(['deCODE','Sun 2018']):
        for k,(grp,col) in enumerate([('matched PAV−',C['pavneg']),('PAV+',C['pav'])]):
            r=md[(md.layer==lay)&(md.group==grp)].iloc[0]; p=r.replicated/r.n; lo,hi=wilson(r.replicated,r.n); x=j*2.6+k
            ax.bar(x,100*p,0.8,color=col,ec='none'); ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.5,capthick=0.7)
            ax.text(x,100*hi+1,f'{100*p:.1f}',ha='center',va='bottom',fontsize=FS)
            out.append(dict(panel='Fig2d',layer=lay,group=grp,replicated=int(r.replicated),n=int(r.n),rate=p,wilson_lo=lo,wilson_hi=hi))
        o,pv,pl=Ps[lay]; ax.text(j*2.6+0.5,43,f'OR {o:.2f}\n{fmt_p(pv,pl)}',ha='center',va='bottom',fontsize=FS)
        out[-1].update(OR=o,P=pv,P_type='Holm-adjusted' if pl=='Holm P' else 'unadjusted (internal layer)')
    ax.set_xticks([0,1,2.6,3.6]); ax.set_xticklabels(['PAV−','PAV+','PAV−','PAV+']); ax.tick_params(axis='x',length=0)
    for xx,t in [(0.5,'deCODE\n701 matched sets'),(3.1,'Sun 2018\n386 matched sets')]: ax.text(xx,-0.11,t,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=FS)
    ax.set_ylim(0,58); ax.set_yticks(range(0,51,10)); ax.set_ylabel('Replicated (%)'); ax.set_xlim(-0.7,4.3)
    # e anchoring forest
    me=S24[S24.panel_sheet=='Fig2e']
    rows=[]
    style={'GTEx eQTL anchor':(C['anchor'],'o',False),'eQTL Catalogue anchor':(C['anchor'],'s',False),'GTEx anchor, strong-signal subset':(C['anchor'],'o',True),
           'wrong-gene eQTL, alone':(C['ctrl'],'o',False),'wrong-gene eQTL, joint model':(C['ctrl'],'o',True)}
    lab={'GTEx eQTL anchor':'GTEx anchor','eQTL Catalogue anchor':'eQTL Catalogue anchor','GTEx anchor, strong-signal subset':'GTEx anchor, strong signals','wrong-gene eQTL, alone':'Wrong-gene eQTL, alone','wrong-gene eQTL, joint model':'Wrong-gene eQTL, joint model'}
    for layer,hdr in [('deCODE (H2a, confirmatory)','deCODE (confirmatory)'),('Sun 2018 (internal)','Sun 2018 (internal)')]:
        rows.append(dict(header=True,label=hdr,color=C['conc']))
        for _,r in me[me.layer==layer].iterrows():
            c,m,op=style[r.estimate]; rows.append(dict(label=lab[r.estimate],est=r.OR,lo=r.lo,hi=r.hi,color=c,marker=m,open=op,ms=3.4))
            out.append(dict(panel='Fig2e',layer=hdr,row=lab[r.estimate],n=int(r.n),OR=r.OR,lo95=r.lo,hi95=r.hi,P=r.P,P_type='unadjusted Wald (gene-clustered)'))
    ax=F.panel('e','Anchoring predicts replication; the wrong-gene control is null in a joint model',64,67,106,65,ml=40,mr=30,mb=12,row=1)
    ys=forest(ax,rows); ax.set_xlim(0.7,6.5); ax.set_xticks([1,2,4]); ax.set_xticklabels(['1','2','4']); ax.set_xlabel('Odds ratio (95% CI), log scale')
    tx=F.ax_mm(64+106-29,67+6,29,65-6-12); tx.axis('off'); tx.set_ylim(ax.get_ylim()); tx.set_xlim(0,1)
    tx.text(0.05,len(rows)-0.35,'OR (95% CI)',fontsize=FS,color=C['sub'],va='bottom')
    for y,r in zip(ys,rows):
        if not r.get('header'): tx.text(0.05,y,f"{r['est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})",fontsize=FS,va='center')
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/Fig2')
    return F,rep,pd.DataFrame(out)
def cs_holm(WB,h): return float(WB['01_Confirmatory_summary'].set_index('hypothesis').loc[h,'P_holm'])

from common import *
TL={'MSPA':'ms_primary','Seer':'ms_ext_primary','Generation Scotland':'ms_gs_primary','HOLBAEK':'ms_hb_primary'}
def ca_trend(k,n,sc=(1,2,3,4)):
    from scipy.stats import norm
    k=np.array(k,float); n=np.array(n,float); s=np.array(sc,float); N=n.sum(); p=k.sum()/N
    T=np.sum(s*(k-n*p)); V=p*(1-p)*(np.sum(n*s*s)-np.sum(n*s)**2/N); z=T/np.sqrt(V); return z,2*norm.sf(abs(z))
def make(R,WB,S27):
    out=[]; F=Fig(170,140,'Fig5')
    h4=WB['H4_detail'].set_index('key').value; a1h=WB['A1_H4_ext_detail'].set_index('key').value; a4=WB['A4_H4_detail']
    # a forest
    fa=S27[S27.panel_sheet=='Fig5a']; rows=[]
    lab={'Pooled non-EUR (H5)':'Pooled, five ancestries','Wrong-gene eQTL (neg. control)':'Wrong-gene eQTL (control)','PAV− stratum':'Within PAV− CS','PAV+ stratum':'Within PAV+ CS',
         '|ΔMAF| ≤ 0.1':'|ΔMAF| ≤ 0.1','EUR lead variant only':'European lead variant only','Lead CS per gene':'Lead CS per gene','PAV itself (same model)':'PAV itself (same model)'}
    hdr={'ancestry':'By ancestry','sens':'Anchoring, sensitivity analyses','control':None,'pav':None}
    seen=set()
    for _,r in fa.iterrows():
        g=r.group
        if g in ('ancestry','sens') and g not in seen: rows.append(dict(header=True,label=hdr[g],color=C['sub'])); seen.add(g)
        col={'main':C['anchor'],'ancestry':C['anchor'],'control':C['ctrl'],'sens':C['anchor'],'pav':C['pav']}[g]
        rows.append(dict(label=lab.get(r.analysis,r.analysis),est=r.OR,lo=r.lo,hi=r.hi,color=col,marker='s' if g=='pav' else 'o',open=g=='sens',ms=4.2 if g=='main' else 3.4))
        out.append(dict(panel='Fig5a',row=lab.get(r.analysis,r.analysis),group=g,n=int(r.n),OR=r.OR,lo95=r.lo,hi95=r.hi,P=r.P,P_type='unadjusted Wald (gene-clustered); pooled row Holm-adjusted P = 8.5e-11 in Table 1'))
    ax=F.panel('a','Anchoring predicts replication in non-European ancestries',0,0,100,70,ml=42,mr=26,mb=10,row=0)
    ys=forest(ax,rows); ax.set_xlim(0.45,3.2); ax.set_xticks([0.5,1,2,3]); ax.set_xticklabels(['0.5','1','2','3']); ax.set_xlabel('Odds ratio for replication (95% CI)')
    tx=F.ax_mm(100-25,6,25,70-6-10); tx.axis('off'); tx.set_ylim(ax.get_ylim()); tx.set_xlim(0,1)
    tx.text(0.04,len(rows)-0.35,'OR (95% CI)',fontsize=FS,color=C['sub'],va='bottom')
    for y,r in zip(ys,rows):
        if not r.get('header'): tx.text(0.04,y,f"{r['est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})",fontsize=FS,va='center')
    # b groups
    fb=S27[S27.panel_sheet=='Fig5b'].copy(); fb['group']=fb.group.replace({'Anch. only':'Anchor only','Anch. +PAV':'Anchor + PAV'})
    d=R[R.xa_replicated.notna()]; chk=d.groupby(['anchor_eqtl','pav']).xa_replicated.agg(['sum','count'])
    for _,r in fb.iterrows(): assert chk.loc[(r.A_eqtl,r.PAV)].tolist()==[r.k,r.n]
    ax=F.panel('b','PAV-only CS replicate as often as anchored CS',102,0,68,70,ml=10,mr=1,mb=10,row=0)
    grp_bars(ax,fb.rename(columns={'k':'k','n':'n'}),ylab='Replicated in non-European\nancestries (%)',ymax=95)
    for _,r in fb.iterrows():
        lo,hi=wilson(r.k,r.n); out.append(dict(panel='Fig5b',group=r.group,replicated=int(r.k),n=int(r.n),rate=r.k/r.n,wilson_lo=lo,wilson_hi=hi))
    ax.text(0.02,0.98,f'{len(d):,} evaluable CS',transform=ax.transAxes,fontsize=FS,va='top',color=C['sub'])
    # c AUC dot plot
    g4=a4.set_index(['layer','set'])
    T=[('MSPA',h4['H4_ms.n_test'],h4['H4_ms.auc_full'],h4['H4_ms.auc_B1'],h4['H4_ms.auc_B0'],h4['H4_ms.delong_full_vs_B1.P_one_sided']),
       ('Seer',a1h['H4_ext.n_test'],a1h['H4_ext.auc_full'],a1h['H4_ext.auc_B1'],a1h['H4_ext.auc_B0'],a1h['H4_ext.delong_full_vs_B1.P_one_sided'])]
    for L_,k in [('Generation Scotland','GS'),('HOLBAEK','HOLBAEK')]:
        r=g4.loc[(k,'even-chr held-out (primary)')]; T.append((L_,r.n,r.auc_full,r.auc_B1,r.auc_B0,r.P_one_sided_vs_B1))
    T.append(('Non-European replication',h4['H4_xa.n_test'],h4['H4_xa.auc_full'],h4['H4_xa.auc_B1'],h4['H4_xa.auc_B0'],h4['H4_xa.delong_full_vs_B1.P_one_sided']))
    ax=F.panel('c','The full model adds little beyond concordance + signal strength',0,74,100,66,ml=42,mr=26,mt=9,mb=16,row=1)
    n=len(T); ys=np.arange(n)[::-1]
    for y,(L_,nt,af,ab1,ab0,p1) in zip(ys,T):
        ax.plot([min(af,ab1,ab0),max(af,ab1,ab0)],[y,y],color='#DDDDDD',lw=2.4,solid_capstyle='round',zorder=1)
        ax.plot(ab0,y,'o',mfc='white',mec='#8C8C8C',ms=4.4,mew=0.9,zorder=2); ax.plot(ab1,y,'D',color=C['b1'],ms=3.8,zorder=3); ax.plot(af,y,'o',color=C['full'],ms=4.4,zorder=4)
        out.append(dict(panel='Fig5c',target=L_,n_heldout=int(nt),AUC_full=float(af),AUC_B1=float(ab1),AUC_B0=float(ab0),dAUC_full_vs_B1=float(af-ab1),P_full_vs_B1=float(p1),P_type='one-sided DeLong, unadjusted'))
    ax.axvline(0.5,color=C['ref'],lw=0.6,ls=(0,(3,2)))
    ax.set_yticks(ys); ax.set_yticklabels([f'{t[0]} (n={int(t[1]):,})' for t in T]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.set_ylim(-0.6,n-0.4); ax.set_xlim(0.25,0.85); ax.set_xticks([0.3,0.4,0.5,0.6,0.7,0.8]); ax.set_xlabel('Held-out AUC')
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=C['full'],ms=4.4,label='full model'),Line2D([],[],marker='D',ls='',color=C['b1'],ms=3.8,label='concordance + signal strength'),
                       Line2D([],[],marker='o',ls='',mfc='white',mec='#8C8C8C',ms=4.4,label='concordance only')],loc='upper left',bbox_to_anchor=(-0.62,-0.20),ncol=3,fontsize=FS,borderaxespad=0,columnspacing=0.9)
    tx=F.ax_mm(100-25,74+9,25,66-9-16); tx.axis('off'); tx.set_ylim(ax.get_ylim()); tx.set_xlim(0,1)
    tx.text(0.04,n-0.45,'Full minus conc.\n+ signal (P)',fontsize=FS,color=C['sub'],va='bottom')
    for y,t in zip(ys,T): tx.text(0.04,y,f'{t[2]-t[3]:+.3f} ({fmt_p(t[5],"").replace("= ","").strip()})',fontsize=FS,va='center')
    # d tiers
    ax=F.panel('d','Support does not follow tier order',102,74,68,66,ml=10,mr=1,mt=9,mb=16,row=1)
    layers=[('MSPA',TL['MSPA']),('Seer',TL['Seer']),('Gen.\nScotland',TL['Generation Scotland']),('HOLBAEK',TL['HOLBAEK']),('Non-Eur.\nrepl.',None)]
    for i,(L_,col) in enumerate(layers):
        if col: dd=R[R[col].isin(['supported','refuted'])]; yv=(dd[col]=='supported')
        else: dd=R[R.xa_replicated.notna()]; yv=dd.xa_replicated==1
        k=[int(yv[dd.tier==t].sum()) for t in (1,2,3,4)]; nn=[int((dd.tier==t).sum()) for t in (1,2,3,4)]; z,p=ca_trend(k,nn)
        for t in range(4):
            x=i*5+t; pr=k[t]/nn[t]; lo,hi=wilson(k[t],nn[t])
            ax.bar(x,100*pr,0.86,color=TIER[t+1]); ax.errorbar(x,100*pr,yerr=[[100*(pr-lo)],[100*(hi-pr)]],fmt='none',ecolor='#333333',elinewidth=0.6,capsize=0.8,capthick=0.6)
            out.append(dict(panel='Fig5d',layer=L_.replace('\n',' '),tier=t+1,positive=k[t],n=nn[t],rate=pr,wilson_lo=lo,wilson_hi=hi,trend_P=p,P_type='two-sided Cochran–Armitage trend, unadjusted'))
        ax.text(i*5+1.5,101,fmt_p(p,'P').replace('= ','=\n'),ha='center',va='bottom',fontsize=FS,linespacing=1.1)
    ax.set_xticks([i*5+1.5 for i in range(5)]); ax.set_xticklabels([l for l,_ in layers]); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0,122); ax.set_yticks(range(0,101,25)); ax.set_ylabel('Supported or replicated (%)'); ax.set_xlim(-0.8,24.0)
    ax.legend(handles=[Patch(color=TIER[t],label=f'Tier {t}') for t in (1,2,3,4)],loc='upper left',bbox_to_anchor=(0,-0.20),ncol=4,fontsize=FS,borderaxespad=0,columnspacing=0.9)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/Fig5')
    return F,rep,pd.DataFrame(out)

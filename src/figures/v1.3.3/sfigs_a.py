from common import *
from scipy.stats import fisher_exact
LAY5={'MSPA':'ms_primary','Seer':'ms_ext_primary','Generation Scotland':'ms_gs_primary','HOLBAEK':'ms_hb_primary'}
def S1(S29):
    out=[]; F=Fig(170,100,'FigS1')
    eps=[('H1_G1_Agtex','Anchor →\nSun 2018\nreplication',C['anchor'],(1.4,6)),('H2a_Agtex','Anchor →\ndeCODE\nreplication',C['anchor'],(1.4,6)),
         ('H2b_PAV_matched','PAV →\ndeCODE repl.\n(matched)',C['pav'],(0.8,5)),('H3i_Aeqtl','Anchor →\nMSPA\nMS support',C['anchor'],(0.25,6)),('H5_Aeqtl','Anchor →\nnon-European\nreplication',C['anchor'],(0.7,2.6))]
    scen=S29[S29.endpoint==eps[0][0]].scenario_label.tolist()
    scen_lab=[S.sup(s.replace('excluding HLA (vacuous)','excluding HLA (no HLA CS)')) for s in scen]
    x0=44; w=(170-x0)/5
    for k,(ep,ttl,col,xl) in enumerate(eps):
        ax=F.panel('abcde'[k],ttl,x0+k*w,0,w,100,ml=1.5,mr=1.5,mt=15,mb=10,row=0)
        F.panels[-1]['T'].set_linespacing(1.1)
        d=S29[S29.endpoint==ep].set_index('scenario_label').loc[scen]; ys=np.arange(len(d))[::-1]
        for y,(lab,r) in zip(ys,d.iterrows()):
            prim=lab=='Primary analysis'
            ax.plot([r.lo,r.hi],[y,y],color=col,lw=1.1); ax.plot(r.OR,y,'o',color=col if prim else 'white',mec=col,ms=4 if prim else 3.4,mew=0.9)
            out.append(dict(panel=f'FigS1{"abcde"[k]}',endpoint=ttl.replace('\n',' ').replace('repl.','replication'),scenario=lab.replace('excluding HLA (vacuous)','excluding HLA (no HLA CS)'),n=int(r.n),OR=r.OR,lo95=r.lo,hi95=r.hi,P=r.P,P_type='unadjusted'))
        ax.axvline(1,color=C['ref'],lw=0.6,ls=(0,(3,2))); ax.axvline(d.loc['Primary analysis','OR'],color=col,lw=0.5,alpha=0.35)
        ax.set_xscale('log'); ax.minorticks_off(); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_xlim(*xl)
        tk=[t for t in ([0.25,1,4] if xl[0]<0.5 else [0.5,1,2,4]) if xl[0]<=t<=xl[1]]; ax.set_xticks(tk); ax.set_xticklabels([f'{t:g}' for t in tk])
        ax.set_yticks(ys); ax.set_yticklabels(scen_lab if k==0 else []); ax.tick_params(axis='y',length=0); ax.set_ylim(-0.6,len(d)-0.4)
        if k>0: ax.spines['left'].set_visible(False)
        if k==2: ax.set_xlabel('Odds ratio (95% CI), log scale')
    F.ax_mm(0,0,1,1).axis('off')
    rep=F.verify(); rep['title_overflow']=[]; F.save(OUT+'/FigS1'); return F,rep,pd.DataFrame(out)
def S2(S30,R):
    out=[]; F=Fig(170,74,'FigS2')
    a=S30[S30.panel_sheet=='FigS2a']
    ax=F.panel('a','564 Seer GWAS files map to UKB-PPP cis genes',0,0,72,74,ml=36,mr=8,mb=20,row=0)
    ys=np.arange(len(a))[::-1]; cols=[MSL['Seer']]*4+['#7FC3B9','#C9E6E1']
    for y,(_,r),c in zip(ys,a.iterrows(),cols):
        ax.barh(y,r.n,0.7,color=c); ax.text(r.n+30,y,f'{int(r.n):,}',va='center',fontsize=FS)
        out.append(dict(panel='FigS2a',step=r.step,n=int(r.n)))
    ax.set_yticks(ys); ax.set_yticklabels(a.step); ax.tick_params(axis='y',length=0); ax.set_xlim(0,2400); ax.set_xticks([0,1000,2000]); ax.set_xlabel('Seer GWAS files'); ax.set_ylim(-0.6,len(a)-0.4)
    ax=F.panel('b','Label outcome differs by MS layer',74,0,52,74,ml=10,mr=1,mb=20,row=0)
    for i,(L_,col) in enumerate(LAY5.items()):
        vc=R[col].value_counts(); sup=int(vc.get('supported',0)); ref=int(vc.get('refuted',0)); ind=int(vc.get('indeterminate',0)+vc.get('indeterminate_underpowered',0)); tot=sup+ref+ind
        b0=0
        for v,c,nm in [(sup,MSL[L_],'supported'),(ref,C['pavneg'],'refuted'),(ind,'#E3E3E3','indeterminate')]:
            ax.bar(i,100*v/tot,0.7,bottom=b0,color=c); b0+=100*v/tot
        ax.text(i,101,f'{tot:,}',ha='center',va='bottom',fontsize=FS)
        out.append(dict(panel='FigS2b',layer=L_,supported=sup,refuted=ref,indeterminate=ind,total=tot))
    ax.set_xticks(range(4)); ax.set_xticklabels(['MSPA','Seer','Gen.\nScotland','HOL-\nBAEK']); ax.tick_params(axis='x',length=0); ax.set_ylim(0,110); ax.set_yticks([0,25,50,75,100]); ax.set_ylabel('CS with a label attempt (%)')
    ax.legend(handles=[Patch(color=MSL['Seer'],label='supported (teal shade = layer)'),Patch(color=C['pavneg'],label='refuted'),Patch(color='#E3E3E3',label='indeterminate')],loc='upper left',bbox_to_anchor=(-0.2,-0.17),ncol=1,fontsize=FS,borderaxespad=0)
    c=S30[S30.panel_sheet=='FigS2c']
    ax=F.panel('c','Seer and MSPA labels agree',128,0,42,74,ml=14,mr=1,mb=20,row=0)
    M=np.array([[c.iloc[0].A1_col0,c.iloc[0].A1_col1],[c.iloc[1].A1_col0,c.iloc[1].A1_col1]])
    cmap=mpl.colors.LinearSegmentedColormap.from_list('teal',['#FFFFFF','#0F4F48'])
    for i in range(2):
        for j in range(2):
            v=M[i,j]; ax.add_patch(Rectangle((j,1-i),1,1,fc=cmap(v/M.max()),ec='white',lw=1)); ax.text(j+0.5,1-i+0.5,f'{int(v)}',ha='center',va='center',fontsize=FS,color='white' if v>30 else C['text'])
            out.append(dict(panel='FigS2c',MSPA=['refuted','supported'][i],Seer=['refuted','supported'][j],n=int(v)))
    ax.set_xlim(0,2); ax.set_ylim(0,2); ax.set_xticks([0.5,1.5]); ax.set_xticklabels(['refuted','supported']); ax.set_yticks([1.5,0.5]); ax.set_yticklabels(['refuted','supported'])
    ax.set_xlabel('Seer label'); ax.set_ylabel('MSPA label'); ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
    ag=(M[0,0]+M[1,1])/M.sum(); ax.text(1,-0.55,f'agreement {100*ag:.1f}%\nκ = 0.89, n = {int(M.sum())}',ha='center',va='top',fontsize=FS,color=C['sub'])
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS2'); return F,rep,pd.DataFrame(out)
def S3(R):
    out=[]; F=Fig(170,150,'FigS3')
    labs=[('MSPA','ms_primary'),('MSPA peptide label','ms_mspa_peptide'),('Seer','ms_ext_primary'),('Generation Scotland','ms_gs_primary'),('HOLBAEK','ms_hb_primary')]
    nm={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}
    h=150/3; w=170/2
    pos=[(0,0),(w,0),(0,h),(w,h),(0,2*h)]
    for k,((L_,col),(px,py)) in enumerate(zip(labs,pos)):
        d=R[R[col].isin(['supported','refuted'])]; y=d[col]=='supported'
        ax=F.panel('abcde'[k],f'{L_} (n = {len(d)})',px,py,w-2,h,ml=11,mr=1,mt=6,mb=8,row=py)
        for j,cc in enumerate(['concordant','discordant','na']):
            for i,(a_,p_) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
                m=(d.platform_concordance==cc)&(d.anchor_eqtl==a_)&(d.pav==p_); kk=int(y[m].sum()); nn=int(m.sum()); x=j*5+i
                if nn==0: ax.text(x,2,'n.d.',ha='center',fontsize=FS,rotation=90); continue
                pr=kk/nn; lo,hi=wilson(kk,nn); ax.bar(x,100*pr,0.82,color=GRP[nm[(a_,p_)]])
                if kk==0: ax.plot([x-0.41,x+0.41],[0.8,0.8],color=GRP[nm[(a_,p_)]],lw=2.2,solid_capstyle='butt')
                ax.errorbar(x,100*pr,yerr=[[100*(pr-lo)],[100*(hi-pr)]],fmt='none',ecolor='#333333',elinewidth=0.6,capsize=0.8,capthick=0.6)
                ax.text(x,103,f'{nn}',ha='center',va='bottom',fontsize=FS,rotation=90)
                out.append(dict(panel=f'FigS3{"abcde"[k]}',label=L_,concordance=cc,eqtl_anchor=a_,PAV=p_,group=nm[(a_,p_)],supported=kk,n=nn,rate=pr,wilson_lo=lo,wilson_hi=hi))
        ax.set_xticks([1.5,6.5,11.5]); ax.set_xticklabels(['concordant','discordant','not assayed']); ax.tick_params(axis='x',length=0)
        ax.set_ylim(0,125); ax.set_yticks([0,25,50,75,100]); ax.set_ylabel('MS-supported (%)'); ax.set_xlim(-0.7,14.7)
    lg=F.ax_mm(w+12,2*h+10,w-14,h-14); lg.axis('off')
    lg.legend(handles=[Patch(color=GRP[g],label=g) for g in GROUPS],loc='upper left',fontsize=FS,title='eQTL anchor × PAV group',title_fontsize=FS,alignment='left')
    lg.text(0,0.30,'Numbers above bars: evaluable CS.\nConcordance refers to Sun 2018 or deCODE\nSomaScan replication; not assayed = gene\nabsent from both SomaScan layers.',fontsize=FS,va='top',color=C['sub'],transform=lg.transAxes)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS3'); return F,rep,pd.DataFrame(out)
def S4(S32):
    out=[]; F=Fig(170,80,'FigS4')
    a=S32[S32.panel_sheet=='FigS4a']; anc=a.ancestry.unique().tolist()
    nm={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}
    ax=F.panel('a','Replication by ancestry: same group ordering in all five',0,0,96,80,ml=10,mr=1,mb=16,row=0)
    for j,an in enumerate(anc):
        for i,(a_,p_) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
            r=a[(a.ancestry==an)&(a.A_eqtl==a_)&(a.PAV==p_)].iloc[0]; x=j*5+i; pr=r.k/r.n; lo,hi=wilson(r.k,r.n)
            ax.bar(x,100*pr,0.82,color=GRP[nm[(a_,p_)]]); ax.errorbar(x,100*pr,yerr=[[100*(pr-lo)],[100*(hi-pr)]],fmt='none',ecolor='#333333',elinewidth=0.6,capsize=0.8,capthick=0.6)
            out.append(dict(panel='FigS4a',ancestry=an,group=nm[(a_,p_)],replicated=int(r.k),n=int(r.n),rate=pr,wilson_lo=lo,wilson_hi=hi))
    ax.set_xticks([j*5+1.5 for j in range(len(anc))]); ax.set_xticklabels([s.replace('Central/South Asian','Central/\nSouth Asian').replace('Middle Eastern','Middle\nEastern').replace('Admixed American','Admixed\nAmerican').replace('East Asian','East\nAsian') for s in anc]); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0,85); ax.set_yticks([0,20,40,60,80]); ax.set_ylabel('Replicated (%)'); ax.set_xlim(-0.7,len(anc)*5-1.3)
    ax.legend(handles=[Patch(color=GRP[g],label=g) for g in GROUPS],loc='upper right',ncol=2,fontsize=FS,borderaxespad=0.1,columnspacing=0.8)
    b=S32[S32.panel_sheet=='FigS4b']
    labmap={'All (pre-specified)':'All CS (pre-specified)','Unweighted OLS':'Unweighted regression','|EUR β| ≤ 1':'|European β| ≤ 1','PAV− only':'PAV− CS only','PAV+ only':'PAV+ CS only','African only':'African only','C/S Asian only':'Central/South Asian only','excl. 204 β-inconsistent CS':'Excluding 204 CS with\ninconsistent β (post hoc)'}
    ax=F.panel('b','Anchored CS have steeper slopes, except PAV+ CS',98,0,72,80,ml=34,mr=2,mb=16,row=0)
    ys=np.arange(len(b))[::-1]
    for y,(_,r) in zip(ys,b.iterrows()):
        ax.plot([r.slope_unanchored,r.slope_anchored],[y,y],color='#CCCCCC',lw=1.0); ax.plot(r.slope_unanchored,y,'o',color=C['neither'],mec='#8C8C8C',ms=4); ax.plot(r.slope_anchored,y,'o',color=C['anchor'],ms=4)
        out.append(dict(panel='FigS4b',analysis=labmap[r.label].replace('\n',' '),slope_unanchored=r.slope_unanchored,slope_anchored=r.slope_anchored,P_interaction=r.P_interaction,n_rows=int(r.n_rows),P_type='unadjusted interaction P'))
    ax.axvline(1,color=C['ref'],lw=0.6,ls=(0,(3,2))); ax.axvline(0,color=C['ref'],lw=0.6)
    ax.set_yticks(ys); ax.set_yticklabels([labmap[l] for l in b.label]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.set_xlim(-0.1,1.9); ax.set_xlabel('Slope, non-European on European β'); ax.set_ylim(-0.6,len(b)-0.4)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=C['anchor'],ms=4,label='anchored'),Line2D([],[],marker='o',ls='',color=C['neither'],mec='#8C8C8C',ms=4,label='unanchored')],loc='upper left',bbox_to_anchor=(-0.05,-0.13),ncol=2,fontsize=FS,borderaxespad=0)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS4'); return F,rep,pd.DataFrame(out)

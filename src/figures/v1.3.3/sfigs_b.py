from common import *
def S5(S33):
    out=[]; F=Fig(170,92,'FigS5')
    a=S33[S33.panel_sheet=='FigS5a']
    ax=F.panel('a','Signal strength gets opposite weights for MS and ancestry targets',0,0,98,92,ml=36,mr=2,mb=20,row=0)
    ys=np.arange(len(a))[::-1]
    for y,(_,r) in zip(ys,a.iterrows()):
        v=[r.coef_std_MSPA_target,r.coef_std_SeerFull_target,r.coef_std_crossancestry_target]
        ax.plot([min(v),max(v)],[y,y],color='#E0E0E0',lw=1.6,zorder=1)
        ax.plot(v[0],y,'o',color=MSL['MSPA'],ms=4,zorder=3); ax.plot(v[1],y,'D',mfc='white',mec=MSL['Seer'],ms=3.6,mew=0.9,zorder=3); ax.plot(v[2],y,'s',color=C['xa'],ms=3.6,zorder=2)
        out.append(dict(panel='FigS5a',feature=r.label,coef_MSPA=v[0],coef_Seer=v[1],coef_nonEuropean=v[2]))
    ax.axvline(0,color=C['ref'],lw=0.6)
    lab=[l.replace('−log10 P (capped)/10','capped −log$_{10}$P / 10') for l in a.label]
    ax.set_yticks(ys); ax.set_yticklabels(lab); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False); ax.set_ylim(-0.6,len(a)-0.4)
    ax.set_xlim(-1.5,1.5); ax.set_xlabel('Standardised coefficient (training chromosomes)')
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=MSL['MSPA'],ms=4,label='MSPA support (n = 99)'),Line2D([],[],marker='D',ls='',mfc='white',mec=MSL['Seer'],ms=3.6,label='Seer support (n = 228)'),
                       Line2D([],[],marker='s',ls='',color=C['xa'],ms=3.6,label='non-European replication (n = 4,572)')],title='training target',title_fontsize=FS,loc='upper left',bbox_to_anchor=(-0.45,-0.14),ncol=3,fontsize=FS,borderaxespad=0,columnspacing=0.8)
    b=S33[S33.panel_sheet=='FigS5b']
    lab={'MSPA target, held-out':'MSPA model → MSPA (held-out)','Seer-full target, held-out':'Seer model → Seer (held-out)','→ held-out cross-ancestry':'MSPA model → non-European\nreplication','→ MSPA peptide label':'MSPA model → MSPA peptide\nlabel','→ Seer discovery overlap':'MSPA model → Seer discovery\noverlap'}
    colr={'MSPA target, held-out':MSL['MSPA'],'Seer-full target, held-out':MSL['Seer'],'→ held-out cross-ancestry':C['xa'],'→ MSPA peptide label':MSL['MSPA'],'→ Seer discovery overlap':MSL['Seer']}
    ax=F.panel('b','The MSPA model fails on other targets',100,0,70,92,ml=41,mr=4,mb=20,row=0)
    ys=np.arange(len(b))[::-1]
    for y,(_,r) in zip(ys,b.iterrows()):
        c=colr[r.evaluation]; ax.plot([0.5,r.AUC],[y,y],color='#CCCCCC',lw=0.8); ax.plot(r.AUC,y,'o',color=c,ms=4.4)
        ax.text(r.AUC+0.04,y+0.28,f'{r.AUC:.2f}',va='center',ha='left',fontsize=FS)
        out.append(dict(panel='FigS5b',evaluation=lab[r.evaluation].replace('\n',' '),AUC=r.AUC))
    ax.axvline(0.5,color=C['ref'],lw=0.6,ls=(0,(3,2))); ax.set_xlim(0,1); ax.set_xticks([0,0.5,1]); ax.set_xticklabels(['0','0.5','1']); ax.set_xlabel('AUC (0.5 = chance)')
    ax.set_yticks(ys); ax.set_yticklabels([lab[e] for e in b.evaluation]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False); ax.set_ylim(-0.6,len(b)-0.4)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS5'); return F,rep,pd.DataFrame(out)
def S6(R):
    out=[]; F=Fig(170,66,'FigS6')
    ax=F.panel('a','Most CS are Tier 3',0,0,46,66,ml=12,mr=1,mb=9,row=0)
    cnt=R.tier.value_counts().sort_index()
    for t in (1,2,3,4):
        ax.bar(t,cnt[t],0.7,color=TIER[t]); ax.text(t,cnt[t]+150,f'{cnt[t]:,}',ha='center',va='bottom',fontsize=FS)
        out.append(dict(panel='FigS6a',tier=t,n_cs=int(cnt[t]),frac=cnt[t]/len(R)))
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.tick_params(axis='x',length=0); ax.set_ylabel('Credible sets'); ax.set_ylim(0,11800); ax.set_yticks([0,2500,5000,7500,10000]); ax.set_yticklabels(['0','2.5k','5k','7.5k','10k'])
    ax=F.panel('b','Tier 4 is commoner among concordant CS',48,0,66,66,ml=19,mr=2,mb=9,row=0)
    rows=[('concordant','Concordant'),('discordant','Discordant'),('na','Not assayed')]
    for y,(k,lab) in zip([2,1,0],rows):
        d=R[R.platform_concordance==k]; left=0
        for t in (1,2,3,4):
            f=(d.tier==t).mean(); ax.barh(y,f,0.62,left=left,color=TIER[t]); 
            if t==4: ax.text(left+f/2,y,f'{100*f:.1f}%',ha='center',va='center',fontsize=FS,color='white' if f>0.1 else C['text'])
            left+=f; out.append(dict(panel='FigS6b',concordance=lab,tier=t,n=int((d.tier==t).sum()),n_total=len(d),frac=f))
    ax.set_yticks([2,1,0]); ax.set_yticklabels([f'{lab}\n(n={int((R.platform_concordance==k).sum()):,})' for k,lab in rows]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.set_xlim(0,1); ax.set_xticks([0,0.5,1]); ax.set_xticklabels(['0','50','100']); ax.set_xlabel('Credible sets (%)'); ax.set_ylim(-0.5,3.1)
    ax.legend(handles=[Patch(color=TIER[t],label=f'T{t}') for t in (1,2,3,4)],loc='upper left',ncol=4,fontsize=FS,borderaxespad=0.1,columnspacing=0.8,handlelength=0.9)
    ax=F.panel('c','Tier 1–2 CS have the strongest signals',116,0,54,66,ml=12,mr=1,mb=9,row=0)
    data=[R.loc[R.tier==t,'lead_mlog10p'].values for t in (1,2,3,4)]
    bp=ax.boxplot(data,positions=[1,2,3,4],widths=0.6,patch_artist=True,showfliers=False,medianprops=dict(color='#1A1A1A',lw=1.0),whiskerprops=dict(lw=0.7),capprops=dict(lw=0.7),boxprops=dict(lw=0.7))
    for patch,t in zip(bp['boxes'],(1,2,3,4)): patch.set_facecolor(TIER[t]); patch.set_edgecolor('#333333')
    for t,dv in zip((1,2,3,4),data): out.append(dict(panel='FigS6c',tier=t,n=len(dv),median=float(np.median(dv)),q1=float(np.percentile(dv,25)),q3=float(np.percentile(dv,75)),min=float(dv.min()),max=float(dv.max())))
    ax.set_yscale('log'); ax.set_yticks([5,10,30,100,300]); ax.set_yticklabels(['5','10','30','100','300']); ax.minorticks_off(); ax.set_ylim(4,420)
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.tick_params(axis='x',length=0); ax.set_ylabel('Lead −log$_{10}$P (capped at 322)')
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS6'); return F,rep,pd.DataFrame(out)
def S7(S35):
    out=[]; F=Fig(170,112,'FigS7')
    c=S35[S35.panel_sheet=='FigS7c']
    lab={'MSPA primary (H3i)':'MSPA primary label','MSPA peptide label':'MSPA peptide label','UKB-PPP rows only (ST10)':'MSPA, UKB-PPP rows only','excl. PAV-peptide rows':'MSPA, excluding PAV-peptide rows','Seer discovery overlap':'Seer discovery-set overlap',
         'wrong-gene eQTL (MSPA)':'Wrong-gene eQTL, MSPA','Seer full (H3i-ext)':'Seer primary label','ext-new (no MSPA overlap)':'Seer, CS without an MSPA label','nominal label (P < 0.05)':'Seer, nominal label (P < 0.05)','CS-wide label':'Seer, CS-wide label','wrong-gene eQTL (Seer full)':'Wrong-gene eQTL, Seer'}
    order=[('mspa','MSPA label variants'),('seer','Seer label variants'),('control','Negative control'),('sens','MSPA, annotation sensitivity')]
    rows=[]
    for g,h in order:
        rows.append(dict(header=True,label=h,color=MSL['MSPA'] if g!='control' else C['ctrl']))
        for _,r in c[c.group==g].iterrows():
            L_=S.sup(lab.get(r.analysis,r.analysis)); col=C['ctrl'] if g=='control' else C['anchor']
            rows.append(dict(label=L_,est=r.OR,lo=r.lo,hi=r.hi,color=col,open=g=='sens',ms=3.4))
            out.append(dict(panel='FigS7a',group=h,analysis=lab.get(r.analysis,r.analysis),n=int(r.n),OR=r.OR,lo95=r.lo,hi95=r.hi,P=r.P,P_type='unadjusted Wald (gene-clustered)'))
    ax=F.panel('a','No anchoring estimate is significantly above 1',0,0,92,112,ml=46,mr=2,mb=10,row=0)
    forest(ax,rows); ax.set_xlim(0.25,5.5); ax.set_xticks([0.25,0.5,1,2,4]); ax.set_xticklabels(['0.25','0.5','1','2','4']); ax.set_xlabel('Odds ratio for MS support (95% CI)')
    d=S35[S35.panel_sheet=='FigS7d']
    labd={'MSPA primary':'MSPA primary','MSPA peptide':'MSPA peptide','MSPA, ST10 only':'MSPA, UKB-PPP rows only','Seer full (A1)':'Seer primary','Seer, ext-new':'Seer, no MSPA label','Seer, nominal':'Seer, nominal','Seer, CS-wide':'Seer, CS-wide'}
    ax=F.panel('b','PAV ORs fall below 1 and concordance ORs above 1',94,0,76,112,ml=30,mr=2,mb=10,row=0)
    labs=list(dict.fromkeys(d.label)); ys=np.arange(len(labs))[::-1]
    for y,L_ in zip(ys,labs):
        for k,(feat,col,mk) in enumerate([('PAV',C['pav'],'s'),('concordant',C['conc'],'o')]):
            r=d[(d.label==L_)&(d.feature==feat)].iloc[0]; yy=y+(0.17 if k==0 else -0.17)
            ax.plot([r.woolf_lo,r.woolf_hi],[yy,yy],color=col,lw=1.1); ax.plot(r.OR,yy,mk,color=col,ms=3.6)
            out.append(dict(panel='FigS7b',label=labd[L_],feature='PAV in 95% CS' if feat=='PAV' else 'platform concordant',OR=r.OR,woolf_lo=r.woolf_lo,woolf_hi=r.woolf_hi,
                            supported_with=int(r.sup_with),not_supported_with=int(r.notsup_with),supported_without=int(r.sup_without),not_supported_without=int(r.notsup_without),P=r.fisher_P,P_type='two-sided Fisher exact, unadjusted'))
    ax.axvline(1,color=C['ref'],lw=0.6,ls=(0,(3,2))); ax.set_xscale('log'); ax.minorticks_off(); ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlim(0.1,6); ax.set_xticks([0.1,0.25,0.5,1,2,4]); ax.set_xticklabels(['0.1','0.25','0.5','1','2','4']); ax.set_xlabel('OR for MS support (Woolf 95% CI)')
    ax.set_yticks(ys); ax.set_yticklabels([labd[l] for l in labs]); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False); ax.set_ylim(-0.6,len(labs)+0.1)
    ax.legend(handles=[Line2D([],[],marker='s',ls='-',color=C['pav'],ms=3.6,label='PAV in 95% CS'),Line2D([],[],marker='o',ls='-',color=C['conc'],ms=3.6,label='concordant')],loc='upper left',ncol=2,fontsize=FS,borderaxespad=0.1,columnspacing=0.8)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS7'); return F,rep,pd.DataFrame(out)
def S8(S36,R):
    out=[]; F=Fig(170,150,'FigS8')
    a=S36[S36.panel_sheet=='FigS8a'].copy()
    a=a.sort_values(['lead_cs_tier','n_cs','symbol'],ascending=[True,False,True])
    ax=F.panel('a','Tier composition of 42 drug targets',0,0,62,150,ml=17,mr=2,mb=10,row=0)
    ys=np.arange(len(a))[::-1]
    for y,(_,r) in zip(ys,a.iterrows()):
        left=0
        for t in (1,2,3,4):
            v=r[f'n_T{t}']
            if v>0: ax.barh(y,v,0.72,left=left,color=TIER[t],ec='white',lw=0.3)
            left+=v
        ax.plot(-0.35,y,'s',color=TIER[int(r.lead_cs_tier)],ms=3.2,clip_on=False)
        out.append(dict(panel='FigS8a',symbol=r.symbol,n_T1=int(r.n_T1),n_T2=int(r.n_T2),n_T3=int(r.n_T3),n_T4=int(r.n_T4),n_cs=int(r.n_cs),lead_cs_tier=int(r.lead_cs_tier),n_cs_tied_at_top_P=int(r.n_cs_tied_at_top_P)))
    ax.set_yticks(ys); ax.set_yticklabels([r.symbol+('*' if r.n_cs_tied_at_top_P>1 else '') for _,r in a.iterrows()],style='italic'); ax.tick_params(axis='y',length=0,pad=8)
    ax.set_ylim(-0.6,len(a)-0.4); ax.set_xlim(0,a.n_cs.max()+0.5); ax.set_xlabel('Credible sets per target'); ax.spines['left'].set_visible(False)
    ax.legend(handles=[Patch(color=TIER[t],label=f'Tier {t}') for t in (1,2,3,4)],loc='lower right',fontsize=FS,borderaxespad=0.1)
    ax.text(1.0,0.115,'square = lead-CS tier\n* tied lead P',transform=ax.transAxes,ha='right',va='bottom',fontsize=FS,color=C['sub'])
    # b matrix
    b=S36[S36.panel_sheet=='FigS8b'].copy()
    Rk=R.set_index(['symbol','lead_variant_id'])
    cols=[('eQTL anchor','anchor_eqtl'),('isoform anchor','anchor_isoform'),('PAV in CS','pav_cs95'),('PAV is lead','pav_lead'),('SomaScan','platform_concordance'),('MSPA','ms_primary'),('Seer','ms_ext_primary'),('Gen. Scotland','gs'),('HOLBAEK','hb'),('non-Eur. repl.','xa_replicated')]
    ax=F.panel('b','Evidence for the eight targets whose lead CS is Tier 4',64,0,106,92,ml=15,mr=1,mt=24,mb=4,row=0)
    def cell(col,v):
        if col in ('anchor_eqtl','anchor_isoform'): return ('yes',C['anchor']) if v==1 else ('no','white')
        if col in ('pav_cs95','pav_lead'): return ('yes',C['pav']) if v==1 else ('no','white')
        if col=='platform_concordance': return {'concordant':('conc.',C['conc']),'discordant':('disc.',C['neither']),'na':('n.a.','#F2F2F2')}[v]
        if col=='xa_replicated': return ('n.a.','#F2F2F2') if pd.isna(v) else (('yes',C['xa']) if v==1 else ('no','white'))
        return {'supported':('sup.',C['ms']),'refuted':('ref.',C['pavneg']),'indeterminate':('ind.','#E3E3E3'),'indeterminate_underpowered':('ind.','#E3E3E3'),'not_assessed':('n.a.','#F2F2F2'),'not_evaluable':('n.a.','#F2F2F2')}[v]
    n=len(b)
    for i,(_,r) in enumerate(b.iterrows()):
        rr=Rk.loc[(r.symbol,r.lead_variant)]; rr=rr.iloc[0] if isinstance(rr,pd.DataFrame) else rr
        vals={'gs':rr.ms_gs_primary,'hb':rr.ms_hb_primary}
        rec=dict(panel='FigS8b',symbol=r.symbol,lead_variant=r.lead_variant,lead_mlog10p=r.lead_mlog10p,tier=int(r.tier),flag=int(r.flag_concordant_unanchored_PAV),best_T1_cs_mlog10p=r.best_T1_cs_mlog10p,best_T1_lead_variant=r.best_T1_lead_variant)
        for j,(cl,key) in enumerate(cols):
            v=vals[key] if key in vals else r[key]; txt,fc=cell(key,v); y=n-1-i
            ax.add_patch(Rectangle((j,y),1,1,fc=fc,ec='#BBBBBB',lw=0.5))
            dark=fc in (C['anchor'],C['pav'],C['conc'],C['ms'],C['pavneg'],C['xa'])
            ax.text(j+0.5,y+0.5,txt,ha='center',va='center',fontsize=FS,color='white' if dark else C['text'])
            rec[cl]=v if not isinstance(v,float) or not np.isnan(v) else None
        ax.text(len(cols)+0.15,n-1-i+0.5,f"{r.best_T1_lead_variant}" if pd.notna(r.best_T1_lead_variant) else '—',va='center',fontsize=FS)
        out.append(rec)
    ax.set_xlim(0,len(cols)+3.9); ax.set_ylim(0,n); ax.set_yticks([n-1-i+0.5 for i in range(n)]); ax.set_yticklabels([s for s in b.symbol],style='italic'); ax.tick_params(length=0)
    ax.set_xticks([j+0.5 for j in range(len(cols))]+[len(cols)+1.9]); ax.set_xticklabels([c for c,_ in cols]+['Tier 1 alternative\n(lead variant)'],rotation=90); ax.xaxis.tick_top()
    [s.set_visible(False) for s in ax.spines.values()]
    # c lead tiers
    ct=a.lead_cs_tier.value_counts().sort_index(); alt=int(S36[S36.panel_sheet=='FigS8c'].T4_lead_with_alternative_T1.iloc[0])
    ax=F.panel('c',f'Lead CS tiers across the 42 targets; {alt} of {int(ct.get(4,0))} Tier 4 leads have a Tier 1 alternative',64,98,106,52,ml=15,mr=1,mb=9,row=1)
    for t in (1,2,3,4):
        ax.bar(t,ct.get(t,0),0.62,color=TIER[t]); ax.text(t,ct.get(t,0)+0.3,str(int(ct.get(t,0))),ha='center',va='bottom',fontsize=FS)
        out.append(dict(panel='FigS8c',lead_tier=t,n_targets=int(ct.get(t,0)),tier4_with_tier1_alternative=alt if t==4 else None))
    ax.set_xticks([1,2,3,4]); ax.set_xticklabels(['T1','T2','T3','T4']); ax.tick_params(axis='x',length=0); ax.set_ylim(0,19); ax.set_yticks([0,5,10,15]); ax.set_ylabel('Targets'); ax.set_xlim(0.4,4.6)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/FigS8'); return F,rep,pd.DataFrame(out)

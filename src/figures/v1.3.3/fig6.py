from common import *
def make(S28,WB):
    out=[]; F=Fig(170,140,'Fig6')
    a=S28[S28.panel=='a']; base=float(WB['A3_P1_proportions'].set_index('endpoint').loc['baseline_resource_tier4','prop'])
    t4=a[a.metric=='tier4_any'].copy(); fl=a[a.metric=='flag_any'].set_index('source')
    t4=pd.concat([t4[t4.source=='pooled'],t4[t4.source!='pooled'].sort_values('n',ascending=False)])
    ax=F.panel('a','15% of matched MR pairs use a Tier 4 instrument',0,0,104,74,ml=40,mr=1,mb=10,row=0)
    ys=np.arange(len(t4))[::-1]
    for y,(_,r) in zip(ys,t4.iterrows()):
        f=fl.loc[r.source]; big=r.source=='pooled'
        ax.plot([r.lo*100,r.hi*100],[y+0.14,y+0.14],color=C['pav'],lw=1.1); ax.plot(r.prop*100,y+0.14,'o',color=C['pav'],ms=4.4 if big else 3.6)
        ax.plot([f.lo*100,f.hi*100],[y-0.16,y-0.16],color=C['pav'],lw=0.8,alpha=0.8); ax.plot(f.prop*100,y-0.16,'D',mfc='white',mec=C['pav'],ms=3.2,mew=0.9)
        out.append(dict(panel='Fig6a',source=r.label,n_pairs=int(r.n),tier4_k=int(r.k),tier4_prop=r.prop,tier4_lo=r.lo,tier4_hi=r.hi,flag_k=int(f.k),flag_prop=f.prop,flag_lo=f.lo,flag_hi=f.hi))
    ax.axvline(base*100,color=C['ref'],lw=0.7,ls=(0,(3,2))); ax.text(base*100+1,len(t4)-0.45,'resource-wide Tier 4 (8.3%)',fontsize=FS,color=C['sub'],va='center')
    def short(l):
        if l.startswith('Pooled'): return 'Pooled, 11 sources'
        if l.startswith('Sun 2023'): return 'Sun 2023 (PCSK9)'
        head=l.split(',')[0]; yr=l.strip()[-4:]; return f'{head} {yr}'
    lab=[f"{short(r.label)} (n={int(r.n):,})" for _,r in t4.iterrows()]
    ax.set_yticks(ys); ax.set_yticklabels(lab); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
    ax.get_yticklabels()[0].set_fontweight('bold')
    ax.set_ylim(-0.6,len(t4)+0.1); ax.set_xlim(-2,102); ax.set_xticks([0,25,50,75,100]); ax.set_xlabel('Pairs using the instrument class (%)')
    ax.legend(handles=[Line2D([],[],marker='o',ls='-',color=C['pav'],ms=3.6,label='any Tier 4 instrument'),Line2D([],[],marker='D',ls='-',mfc='white',mec=C['pav'],color=C['pav'],ms=3.2,label='flagged: concordant, unanchored, PAV')],
              loc='upper right',bbox_to_anchor=(1,0.90),fontsize=FS,borderaxespad=0.1)
    # b
    b=S28[S28.panel=='b']
    ax=F.panel('b','Not enriched among significant results',106,0,64,74,ml=10,mr=1,mb=10,row=0)
    for j,(m,lbl) in enumerate([('tier4_any','any Tier 4'),('flag_any','flagged')]):
        for k,sig in enumerate([0.0,1.0]):
            r=b[(b.metric==m)&(b.significant==sig)].iloc[0]; x=j*2.7+k; p=r.k/r.n; lo,hi=wilson(r.k,r.n)
            fc=(C['pav'] if m=='tier4_any' else 'white') if sig==1 else ('#F3C7A8' if m=='tier4_any' else 'white')
            ax.bar(x,100*p,0.82,color=fc,ec=C['pav'],lw=0.9,hatch=None if sig==1 or m=='tier4_any' else None)
            ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.4,capthick=0.7)
            out.append(dict(panel='Fig6b',metric=lbl,significant=bool(sig),k=int(r.k),n=int(r.n),prop=p,wilson_lo=lo,wilson_hi=hi,OR=r.logit_or,OR_lo=r.logit_lo,OR_hi=r.logit_hi,P_holm=r.holm_p,P_type='Holm-adjusted (A3 family); logistic regression, source fixed effects, gene-clustered SEs'))
        ax.text(j*2.7+0.5,19.5,f'OR {r.logit_or:.2f}\n({r.logit_lo:.2f}–{r.logit_hi:.2f})\nHolm P = {r.holm_p:.2f}',ha='center',va='bottom',fontsize=FS,linespacing=1.15)
    ax.set_xticks([0,1,2.7,3.7]); ax.set_xticklabels(['not sig.\nn=19,365','sig.\nn=662','not sig.','sig.']); ax.tick_params(axis='x',length=0)
    for xx,tl in [(0.5,'any Tier 4'),(3.2,'flagged')]: ax.text(xx,-0.17,tl,transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=FS,color=C['sub'])
    ax.set_ylim(0,28); ax.set_yticks([0,5,10,15,20]); ax.set_ylabel('Pairs using the instrument class (%)'); ax.set_xlim(-0.6,4.3)
    # c
    c=S28[S28.panel=='c'].sort_values(['n_pairs','n_sources','symbols'],ascending=[False,False,True])
    ax=F.panel('c',f"{int(c.alt_T12.sum())} of {len(c)} genes behind significant Tier 4-based results have a Tier 1–2 alternative",0,78,170,62,ml=10,mr=1,mt=7,mb=17,row=1)
    x=np.arange(len(c))
    ax.bar(x,c.n_pairs,0.78,color=[C['anchor'] if v==1 else C['neither'] for v in c.alt_T12])
    for xi,(_,r) in zip(x,c.iterrows()):
        if r.n_pairs>=2: ax.text(xi,r.n_pairs+0.3,f'{int(r.n_sources)}',ha='center',va='bottom',fontsize=FS,color=C['sub'])
        out.append(dict(panel='Fig6c',gene=r.symbols,significant_pairs=int(r.n_pairs),n_sources=int(r.n_sources),tier12_alternative=bool(r.alt_T12),flagged=bool(r.flag)))
    ax.set_xticks(x); ax.set_xticklabels(c.symbols,rotation=90,style='italic'); ax.tick_params(axis='x',length=0)
    ax.set_xlim(-0.7,len(c)-0.3); ax.set_ylim(0,21); ax.set_yticks([0,5,10,15,20]); ax.set_ylabel('Significant pairs')
    ax.legend(handles=[Patch(color=C['anchor'],label='Tier 1–2 alternative CS exists'),Patch(color=C['neither'],label='no Tier 1–2 alternative')],loc='upper right',fontsize=FS,borderaxespad=0.1)
    ax.text(0.995,0.70,'numbers above bars: source studies',transform=ax.transAxes,ha='right',va='top',fontsize=FS,color=C['sub'])
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/Fig6')
    return F,rep,pd.DataFrame(out)

from common import *
from scipy.stats import fisher_exact
from sklearn.metrics import cohen_kappa_score as kap
LAY={'Generation Scotland':'ms_gs_primary','HOLBAEK':'ms_hb_primary','MSPA':'ms_primary','Seer':'ms_ext_primary'}
SHORT={'Generation Scotland':'Gen.\nScotland','HOLBAEK':'HOLBAEK','MSPA':'MSPA','Seer':'Seer'}
def make(R):
    out=[]; F=Fig(170,144,'Fig4')
    lab={k:R[v].map({'supported':1,'refuted':0}) for k,v in LAY.items()}
    names=list(LAY)
    # a kappa
    ax=F.panel('a','Seer-study labels barely match cohorts',0,0,54,64,ml=15,mr=1,mt=6,mb=13,row=0)
    cmap=mpl.colors.LinearSegmentedColormap.from_list('teal',['#FFFFFF','#0F4F48'])
    for i in range(4):
        for j in range(4):
            if j>=i: continue
            m=lab[names[i]].notna()&lab[names[j]].notna(); kv=kap(lab[names[i]][m],lab[names[j]][m]); n=int(m.sum()); ag=(lab[names[i]][m]==lab[names[j]][m]).mean()
            col=cmap(max(kv,0)); ax.add_patch(Rectangle((j,3-i),1,1,fc=col,ec='white',lw=1))
            ax.text(j+0.5,3-i+0.5,f'{kv:.2f}\nn={n}',ha='center',va='center',fontsize=FS,color='white' if kv>0.5 else C['text'])
            out.append(dict(panel='Fig4a',layer_1=names[i],layer_2=names[j],n=n,kappa=kv,agreement=ag))
    ax.set_xlim(0,3); ax.set_ylim(0,3); ax.set_aspect('auto')
    ax.set_xticks([0.5,1.5,2.5]); ax.set_xticklabels([SHORT[n] for n in names[:3]]); ax.set_yticks([0.5,1.5,2.5]); ax.set_yticklabels([SHORT[n] for n in names[1:]][::-1])
    ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
    ax.text(0.5,-0.25,"Cohen's κ (n = CS labelled in both)",transform=ax.transAxes,ha='center',va='top',fontsize=FS,color=C['sub'])
    def feat_panel(ax,feat,colon,coloff,labon,laboff,pname,ly=-0.30):
        for i,L_ in enumerate(names):
            d=R[R[LAY[L_]].isin(['supported','refuted'])]; y=(d[LAY[L_]]=='supported').astype(int)
            if feat=='conc': ok=d.platform_concordance!='na'; d=d[ok]; y=y[ok]; m=d.platform_concordance=='concordant'
            elif feat=='anchor': m=d.anchor_eqtl==1
            else: m=d.pav==1
            k1,n1,k0,n0=int(y[m].sum()),int(m.sum()),int(y[~m].sum()),int((~m).sum()); P=fisher_exact([[k1,n1-k1],[k0,n0-k0]])[1]
            for k,(kk,nn,col) in enumerate([(k1,n1,colon),(k0,n0,coloff)]):
                x=i*2.6+k; p=kk/nn; lo,hi=wilson(kk,nn)
                ax.bar(x,100*p,0.82,color=col); ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.2,capthick=0.7)
            ax.text(i*2.6+0.5,101,fmt_p(P,'P').replace('= ','=\n'),ha='center',va='bottom',fontsize=FS,linespacing=1.1)
            out.append(dict(panel=pname,layer=L_,feature=feat,k_with=k1,n_with=n1,rate_with=k1/n1,k_without=k0,n_without=n0,rate_without=k0/n0,P=P,P_type='two-sided Fisher exact (descriptive)'))
        ax.set_xticks([i*2.6+0.5 for i in range(4)]); ax.set_xticklabels([SHORT[n] for n in names]); ax.tick_params(axis='x',length=0)
        ax.set_ylim(0,122); ax.set_yticks(range(0,101,25)); ax.set_xlim(-0.6,4*2.6-0.9); ax.axvline(2*2.6-0.8,color='#BBBBBB',lw=0.6,ls=(0,(2,2)))
        ax.legend(handles=[Patch(color=colon,label=labon),Patch(color=coloff,label=laboff)],loc='lower left',bbox_to_anchor=(0,ly),ncol=2,fontsize=FS,borderaxespad=0,columnspacing=1.0)
    ax=F.panel('b','Concordance tracks independent labels',54,0,59,64,ml=11,mr=1,mt=6,mb=13,row=0)
    feat_panel(ax,'conc',C['conc'],C['neither'],'platform-concordant','discordant','Fig4b'); ax.set_ylabel('MS-supported (%)')
    ax=F.panel('c','PAV CS supported less often everywhere',113,0,57,64,ml=9,mr=1,mt=6,mb=13,row=0)
    feat_panel(ax,'pav',C['pav'],C['pavneg'],'PAV+','PAV−','Fig4c')
    ax=F.panel('d','Anchoring does not predict MS support',0,70,60,74,ml=11,mr=1,mt=6,mb=17,row=1)
    feat_panel(ax,'anchor',C['anchor'],C['neither'],'eQTL-anchored','unanchored','Fig4d',ly=-0.27); ax.set_ylabel('MS-supported (%)')
    # e: concordant CS by group
    ax=F.panel('e','Concordant CS: no consistent deficit for unanchored PAV-containing CS',62,70,108,74,ml=9,mr=1,mt=6,mb=17,row=1)
    G=[('Anchored',C['anchor']),('No anchor, no PAV',C['neither']),('No anchor, PAV',C['pav'])]
    for i,L_ in enumerate(names):
        d=R[R[LAY[L_]].isin(['supported','refuted'])&(R.platform_concordance=='concordant')]; y=(d[LAY[L_]]=='supported').astype(int)
        ms={'Anchored':d.anchor_eqtl==1,'No anchor, no PAV':(d.anchor_eqtl==0)&(d.pav==0),'No anchor, PAV':(d.anchor_eqtl==0)&(d.pav==1)}
        kn={g:(int(y[m].sum()),int(m.sum())) for g,m in ms.items()}
        for k,(g,col) in enumerate(G):
            kk,nn=kn[g]; x=i*3.8+k; p=kk/nn; lo,hi=wilson(kk,nn)
            ax.bar(x,100*p,0.84,color=col); ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.2,capthick=0.7)
            if p>0.3: ax.text(x,4,f'{kk}/{nn}',ha='center',va='bottom',fontsize=FS,rotation=90,color='white' if col!=C['neither'] else C['text'])
            else: ax.text(x,100*hi+2,f'{kk}/{nn}',ha='center',va='bottom',fontsize=FS,rotation=90,color=C['text'])
            out.append(dict(panel='Fig4e',layer=L_,group=g,supported=kk,n=nn,rate=p,wilson_lo=lo,wilson_hi=hi))
        (a1,n1),(a0,n0)=kn['No anchor, PAV'],kn['Anchored']
        orr,P=fisher_exact([[a1,n1-a1],[a0,n0-a0]])
        ortxt='OR not finite' if not np.isfinite(orr) else f'OR {orr:.2f}'
        ax.text(i*3.8+1,103,f'{ortxt}\n{fmt_p(P,"P")}',ha='center',va='bottom',fontsize=FS)
        out[-1].update(CNP_vs_CA_OR=orr,P=P,P_type='two-sided Fisher exact (pre-specified test for MSPA, Seer, GS; secondary for HOLBAEK)')
    ax.set_xticks([i*3.8+1 for i in range(4)]); ax.set_xticklabels([SHORT[n].replace('\n',' ') for n in names]); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0,122); ax.set_yticks(range(0,101,25)); ax.set_xlim(-0.7,4*3.8-1.1); ax.axvline(2*3.8-0.9,color='#BBBBBB',lw=0.6,ls=(0,(2,2)))
    ax.legend(handles=[Patch(color=c,label=g) for g,c in G],loc='lower left',bbox_to_anchor=(0,-0.27),ncol=3,fontsize=FS,borderaxespad=0)
    ax.text(0.5*(2*3.8-0.9-0.7)/(4*3.8-1.1+0.7),1.0,'',transform=ax.transAxes)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/Fig4')
    return F,rep,pd.DataFrame(out)

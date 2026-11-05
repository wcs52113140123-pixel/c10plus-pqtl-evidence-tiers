from common import *
import json
def make(WB,S25):
    g=lambda sh,k: WB[sh].set_index('key').loc[k,'value']
    out=[]; F=Fig(170,134,'Fig3')
    L={'Sun 2018':'H1_detail','deCODE':'H2_detail'}
    # a: r by PAV
    ax=F.panel('a','PAV+ pairs agree less in effect size',0,0,56,62,ml=10,mr=2,row=0)
    for j,lay in enumerate(['deCODE','Sun 2018']):
        for k,(pv,col) in enumerate([('PAV−',C['pavneg']),('PAV+',C['pav'])]):
            r=S25[(S25.panel_sheet=='Fig3a')&(S25.layer==lay)&(S25.PAV==pv)].iloc[0]; x=j*2.4+k*0.9
            ax.plot([x,x],[r.fisherCI_lo,r.fisherCI_hi],color=col,lw=1.1); ax.plot(x,r.pearson_r,'o',color=col,ms=4.2)
            ax.plot(x+0.28,r.r_same_direction_only,'o',mfc='white',mec=col,ms=3.6,mew=0.9)
            out.append(dict(panel='Fig3a',layer=lay,PAV=pv,n_pairs=int(r.n_pairs),pearson_r=r.pearson_r,fisherCI_lo=r.fisherCI_lo,fisherCI_hi=r.fisherCI_hi,r_same_direction_only=r.r_same_direction_only,n_same_direction=int(r.n_same_direction)))
        z,pz=json.loads(g(L[lay],'effect_corr.fisher_z_PAV1_lt_PAV0')); ax.text(j*2.4+0.6,0.98,'one-sided\n'+fmt_p(pz,'P'),ha='center',va='bottom',fontsize=FS)
        out[-1].update(P=pz,P_type='one-sided Fisher z (PAV+ < PAV−), unadjusted')
    ax.set_xticks([0.45,2.85]); ax.set_xticklabels(['deCODE','Sun 2018']); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0.2,1.2); ax.set_yticks([0.2,0.4,0.6,0.8,1.0]); ax.set_ylabel('Cross-platform β correlation (r)'); ax.set_xlim(-0.5,3.9)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=C['pavneg'],ms=4,label='PAV−'),Line2D([],[],marker='o',ls='',color=C['pav'],ms=4,label='PAV+'),
                       Line2D([],[],marker='o',ls='',mfc='white',mec='#555555',ms=3.6,label='same-direction pairs only')],loc='lower left',fontsize=FS,borderaxespad=0.1)
    # b: opposite sign
    ax=F.panel('b','PAV+ pairs more often flip sign',58,0,50,62,ml=10,mr=2,row=0)
    for j,lay in enumerate(['deCODE','Sun 2018']):
        for k,(pv,col) in enumerate([('PAV−',C['pavneg']),('PAV+',C['pav'])]):
            r=S25[(S25.panel_sheet=='Fig3b')&(S25.layer==lay)&(S25.PAV==pv)].iloc[0]; x=j*2.4+k*0.9; p=r.n_opposite/r.n_pairs; lo,hi=wilson(r.n_opposite,r.n_pairs)
            ax.bar(x,100*p,0.8,color=col); ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.5,capthick=0.7)
            ax.text(x,100*hi+0.6,f'{int(r.n_opposite)}/{int(r.n_pairs)}',ha='center',va='bottom',fontsize=FS,rotation=90)
            out.append(dict(panel='Fig3b',layer=lay,PAV=pv,n_opposite=int(r.n_opposite),n_pairs=int(r.n_pairs),frac=p,wilson_lo=lo,wilson_hi=hi))
        pf=float(g(L[lay],'effect_corr.opp_fisher_PAV1_gt.P_one_sided')); ax.text(j*2.4+0.45,35.5,'one-sided\n'+fmt_p(pf,'P'),ha='center',va='bottom',fontsize=FS)
        out[-1].update(P=pf,P_type='one-sided Fisher exact (PAV+ > PAV−), unadjusted')
    ax.set_xticks([0.45,2.85]); ax.set_xticklabels(['deCODE','Sun 2018']); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0,44); ax.set_yticks(range(0,31,10)); ax.set_ylabel('Opposite-sign pairs (%)'); ax.set_xlim(-0.6,3.9)
    # c: by anchor x PAV
    ax=F.panel('c','deCODE: reversals mostly unanchored PAV',110,0,60,62,ml=10,mr=1,row=0)
    nm={(0,0):'Neither',(0,1):'PAV only',(1,0):'Anchor only',(1,1):'Anchor + PAV'}
    for j,lay in enumerate(['deCODE','Sun 2018']):
        sub=S25[(S25.panel_sheet=='Fig3c')&(S25.layer==lay)]
        for i,(a_,p_) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
            r=sub[(sub.A_eqtl==a_)&(sub.PAV.astype(int)==p_)].iloc[0]; x=j*4.8+i; p=r.n_opposite/r.n_pairs; lo,hi=wilson(r.n_opposite,r.n_pairs)
            ax.bar(x,100*p,0.8,color=GRP[nm[(a_,p_)]]); ax.errorbar(x,100*p,yerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.2,capthick=0.7)
            ax.text(x,100*hi+0.6,f'{int(r.n_opposite)}/{int(r.n_pairs)}',ha='center',va='bottom',fontsize=FS,rotation=90)
            out.append(dict(panel='Fig3c',layer=lay,group=nm[(a_,p_)],n_opposite=int(r.n_opposite),n_pairs=int(r.n_pairs),frac=p,wilson_lo=lo,wilson_hi=hi))
    ax.set_xticks([1.5,6.3]); ax.set_xticklabels(['deCODE','Sun 2018']); ax.tick_params(axis='x',length=0)
    ax.set_ylim(0,40); ax.set_yticks(range(0,31,10)); ax.set_ylabel('Opposite-sign pairs (%)'); ax.set_xlim(-0.7,8.5)
    ax.legend(handles=[Patch(color=GRP[k],label=k) for k in GROUPS],loc='upper left',fontsize=FS,ncol=1,borderaxespad=0.1)
    # d: cohort forest
    d=S25[S25.panel_sheet=='Fig3d']; rows=[]
    coh={'deCODE (H2)':'deCODE','ARIC EA':'ARIC, European ancestry','Fenland':'Fenland','ARIC AA (secondary)':'ARIC, African ancestry','Random-effects meta (deCODE+ARIC EA+Fenland)':'Meta-analysis (3 European cohorts)'}
    for ep,hdr,col in [('H2a: eQTL anchor','eQTL anchor',C['anchor']),('H2b: PAV (power-matched)','PAV, after power matching',C['pav'])]:
        rows.append(dict(header=True,label=hdr,color=col))
        for _,r in d[d.endpoint==ep].iterrows():
            mk='D' if 'meta' in r.cohort else ('s' if col==C['pav'] else 'o')
            rows.append(dict(label=coh[r.cohort],est=r.OR,lo=r.lo,hi=r.hi,color=col,marker=mk,open='AA' in r.cohort,ms=3.8 if mk=='D' else 3.4))
            out.append(dict(panel='Fig3d',endpoint=hdr,cohort=coh[r.cohort],OR=r.OR,lo95=r.lo,hi95=r.hi))
    ax=F.panel('d','Both associations reproduce in every SomaScan v4 cohort',0,66,88,68,ml=44,mr=1,mb=10,row=1)
    ys=forest(ax,rows); ax.set_xlim(1,5.5); ax.set_xticks([1,2,3,4,5]); ax.set_xticklabels(['1','2','3','4','5']); ax.set_xlabel('Odds ratio for replication (95% CI)')
    # e, f per cohort: horizontal, shared cohort axis
    ef=S25[S25.panel_sheet=='Fig3ef']; cohorts=['deCODE','ARIC EA','Fenland','ARIC AA']; yl=['deCODE','ARIC, European','Fenland','ARIC, African']
    axe=F.panel('e','Lower r for PAV+ pairs…',90,66,44,68,ml=21,mr=1,mb=10,row=1)
    axf=F.panel('f','…and flip sign more',136,66,34,68,ml=2,mr=1,mb=10,row=1)
    for i,cn in enumerate(cohorts):
        y0=3-i
        for k,(pv,col) in enumerate([(0,C['pavneg']),(1,C['pav'])]):
            r=ef[(ef.cohort==cn)&(ef.PAV.astype(int)==pv)].iloc[0]; y=y0+(0.5-k)*0.36; op=cn=='ARIC AA'
            axe.plot([r.fisherCI_lo,r.fisherCI_hi],[y,y],color=col,lw=1.1); axe.plot(r.pearson_r,y,'o',color=col,ms=3.8,mfc='white' if op else col,mew=0.9)
            out.append(dict(panel='Fig3e',cohort=cn,PAV=['PAV−','PAV+'][pv],n_pairs=int(r.n_pairs),pearson_r=r.pearson_r,fisherCI_lo=r.fisherCI_lo,fisherCI_hi=r.fisherCI_hi))
            p=r.n_opposite/r.n_pairs; lo,hi=wilson(r.n_opposite,r.n_pairs)
            axf.barh(y,100*p,0.34,color='white' if op else col,ec=col,lw=0.9); axf.errorbar(100*p,y,xerr=[[100*(p-lo)],[100*(hi-p)]],fmt='none',ecolor='#333333',elinewidth=0.7,capsize=1.0,capthick=0.7)
            out.append(dict(panel='Fig3f',cohort=cn,PAV=['PAV−','PAV+'][pv],n_opposite=int(r.n_opposite),n_pairs=int(r.n_pairs),frac=p,wilson_lo=lo,wilson_hi=hi))
    for a_ in (axe,axf): a_.set_ylim(-0.6,5.5); a_.set_yticks(range(4)); a_.tick_params(axis='y',length=0)
    axe.set_yticklabels(yl[::-1]); axf.set_yticklabels([]); axe.set_xlim(0.25,0.9); axe.set_xticks([0.3,0.5,0.7,0.9]); axe.set_xlabel('β correlation (r)')
    axf.set_xlim(0,27); axf.set_xticks([0,10,20]); axf.set_xlabel('Opposite sign (%)')
    axf.legend(handles=[Patch(color=C['pavneg'],label='PAV−'),Patch(color=C['pav'],label='PAV+'),Patch(fc='white',ec='#555555',label='open: African\nancestry')],loc='upper left',fontsize=FS,ncol=1,borderaxespad=0.1)
    rep=F.verify(); rep['title_overflow']=title_fit(F); F.save(OUT+'/Fig3')
    return F,rep,pd.DataFrame(out)

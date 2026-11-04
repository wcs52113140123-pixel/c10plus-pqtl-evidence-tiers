"""H3: mass-spectrometry truth layer (Suhre 2025 MSPA ST10/ST9 primary; Seer MS discovery secondary)."""
import sys; sys.path.insert(0,'.')
from c10common import *
def h3_block(M,label,adj=('mlp',),allow_stop=True):
    o={}
    T=M[M[label].notna()].copy(); T['y']=T[label].astype(int)
    o['n_evaluable']=int(len(T)); o['n_supported']=int(T.y.sum()); o['n_genes']=int(T.ga.nunique())
    o['H3i_fisher']=fisher(T,'A_eqtl','y')
    try: o['H3i_logit'],_=glm(T,'y',['A_eqtl',*adj])
    except Exception as e: o['H3i_logit']={'error':str(e)}
    try: o['H3i_negctrl_wrong'],_=glm(T,'y',['W_gtex',*adj])
    except Exception as e: o['H3i_negctrl_wrong']={'error':str(e)}
    o['H3i_negctrl_wrong_fisher']=fisher(T,'W_gtex','y')
    try: o['A_gtex_only_logit'],_=glm(T,'y',['A_gtex',*adj])
    except Exception as e: o['A_gtex_only_logit']={'error':str(e)}
    o['A_iso_fisher']=fisher(T,'A_iso','y'); o['PAV_fisher']=fisher(T,'PAV','y')
    o['concordant_fisher']=fisher(T[T.conc_cat!='na'],'concordant','y')
    C=T[T.concordant==1].copy()
    C['grp']=np.select([C.A_eqtl==1,(C.A_eqtl==0)&(C.PAV==1),(C.A_eqtl==0)&(C.PAV==0)],['CA','CNP','CN0'],'x')
    o['concordant_groups']=C.groupby('grp').y.agg(['sum','size','mean']).reset_index().to_dict('records')
    def cmp(df,g1,g0):
        d=df[df.grp.isin([g1,g0])].copy(); d['g']=(d.grp==g1).astype(int)
        r=fisher(d,'g','y'); r['comparison']=f'{g1} vs {g0}'
        return r
    o['H3ii_CNP_vs_CA']=cmp(C,'CNP','CA'); o['CN0_vs_CA']=cmp(C,'CN0','CA')
    n_cnp=int((C.grp=='CNP').sum()); n_ca=int((C.grp=='CA').sum())
    o['H3ii_stop_rule_triggered']=bool(allow_stop and (len(T)<60 or n_cnp<10 or n_ca<10))
    d=C[C.grp.isin(['CNP','CA'])].copy(); d['CNP']=(d.grp=='CNP').astype(int)
    try: o['H3ii_logit_adj'],_=glm(d,'y',['CNP',*adj])
    except Exception as e: o['H3ii_logit_adj']={'error':str(e)}
    # discordant & anchored vs CNP
    T['grp2']=np.select([(T.conc_cat=='discordant')&(T.A_eqtl==1),(T.concordant==1)&(T.A_eqtl==0)&(T.PAV==1)],['DA','CNP'],'x')
    d2=T[T.grp2.isin(['DA','CNP'])].copy(); d2['g']=(d2.grp2=='DA').astype(int); o['DA_vs_CNP']=fisher(d2,'g','y') if len(d2) else {}
    # direction-aware concordance sensitivity
    Cd=T[T.concordant_dir==1].copy(); Cd['grp']=np.select([Cd.A_eqtl==1,(Cd.A_eqtl==0)&(Cd.PAV==1),(Cd.A_eqtl==0)&(Cd.PAV==0)],['CA','CNP','CN0'],'x')
    o['H3ii_CNP_vs_CA_samedir']=cmp(Cd,'CNP','CA')
    # full cross-classification (descriptive)
    o['grid']=T.groupby(['conc_cat','A_eqtl','PAV']).y.agg(['sum','size','mean']).reset_index().to_dict('records')
    return o,T
if __name__=="__main__":
    M=derive(load()); res={}
    res['primary'],T=h3_block(M,'ms_primary')
    res['mspa'],_=h3_block(M,'ms_mspa',allow_stop=False)
    res['seer'],_=h3_block(M,'ms_seer',adj=('mlp','absb','maf'),allow_stop=False)
    # Seer restricted to the strongest UKB cs per gene (power-comparable)
    top=M.sort_values('mlog10p',ascending=False).drop_duplicates('ga')
    res['seer_leadcs'],_=h3_block(top,'ms_seer',adj=('mlp','absb','maf'),allow_stop=False)
    Mx=M[~(M.ms_pep_pav.fillna(False).astype(bool))]
    res['primary_excl_pep_pav'],_=h3_block(Mx,'ms_primary',allow_stop=False)
    res['primary_ST10_only'],_=h3_block(M[M.ms_sheet=='ST10'],'ms_primary',allow_stop=False)
    T.to_csv(R+'h3_ms_primary_table.csv',index=False)
    dump(res,'h3.json')
    p=res['primary']; print(json.dumps({k:p[k] for k in ['n_evaluable','n_supported','H3i_fisher','H3i_logit','concordant_groups','H3ii_CNP_vs_CA','H3ii_stop_rule_triggered']},indent=1,default=str))

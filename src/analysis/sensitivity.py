"""Step 5 sensitivity analyses (pre-specified list, prereg section 5) on H1-G1, H2a, H2b, H3i/H3ii, H5 primary estimates."""
import sys; sys.path.insert(0,'.')
from c10common import *
def est(M):
    o={}
    H1=M[M.v3]; H2=M[M.v4&~M.v3]
    def g(D,y,xs,**k):
        try: r,_=glm(D,y,xs,**k); v=list(r.values())[0]; return [v['OR'],v['lo'],v['hi'],v['P'],v['n']]
        except Exception as e: return [np.nan]*4+[len(D)]
    o['H1_G1_Agtex']=g(H1,'rep_sun',['A_gtex','mlp','logn95','PAV','sun_has_cis'])
    try:
        _,mc=glm(H1,'rep_sun',['A_gtex','W_gtex','mlp','logn95','PAV','sun_has_cis'],terms=['A_gtex']); o['H1_joint_diff_P']=diff_test(mc,'A_gtex','W_gtex')['P_one_sided']
    except Exception: o['H1_joint_diff_P']=np.nan
    o['H2a_Agtex']=g(H2,'rep_decode',['A_gtex','mlp','logn95','PAV','decode_has_cis'])
    try:
        _,mc=glm(H2,'rep_decode',['A_gtex','W_gtex','mlp','logn95','PAV','decode_has_cis'],terms=['A_gtex']); o['H2a_joint_diff_P']=diff_test(mc,'A_gtex','W_gtex')['P_one_sided']
    except Exception: o['H2a_joint_diff_P']=np.nan
    Mt=pav_match(H2,'rep_decode','decode_has_cis')
    o['H2b_PAV_matched']=g(Mt,'rep_decode',['PAV'],cluster='set_id') if Mt is not None else [np.nan]*5
    T=M[M.ms_primary.notna()].copy(); T['y']=T.ms_primary.astype(int)
    o['H3i_Aeqtl']=g(T,'y',['A_eqtl','mlp'])
    C=T[T.concordant==1]; cnp=C[(C.A_eqtl==0)&(C.PAV==1)]; ca=C[C.A_eqtl==1]
    if len(cnp) and len(ca):
        from scipy.stats import fisher_exact
        OR,p=fisher_exact([[cnp.y.sum(),len(cnp)-cnp.y.sum()],[ca.y.sum(),len(ca)-ca.y.sum()]])
        o['H3ii_CNP_vs_CA']=[OR,np.nan,np.nan,p,len(cnp)+len(ca)]; o['H3ii_rates']=[float(cnp.y.mean()),len(cnp),float(ca.y.mean()),len(ca)]
    D=M[M.xa_rep.notna()].copy(); D['logse']=np.log(D.xa_se)
    o['H5_Aeqtl']=g(D,'xa_rep',['A_eqtl','mlp','absb','maf','logse','logn95'])
    return o
if __name__=="__main__":
    base=load(); res={}
    sc={'primary':lambda: derive(base),
        'h4_0.9':lambda: derive(base,thr=0.9),'h4_0.95':lambda: derive(base,thr=0.95),
        'prior_p12_5e-6':lambda: derive(base,var='_p12_5e6'),'prior_p12_1e-6':lambda: derive(base,var='_p12_1e6'),
        'coloc_intersection':lambda: derive(base,var='_inter'),
        'window_lead_TSS_le_1Mb':lambda: derive(base)[lambda d: d.dist_tss<=1e6],
        'window_lead_TSS_le_500kb':lambda: derive(base)[lambda d: d.dist_tss<=5e5],
        'lead_cs_per_gene':lambda: derive(base).sort_values('mlog10p',ascending=False).drop_duplicates('ga'),
        'exclude_HLA':lambda: derive(base)[lambda d: d.hla==0],
        'exclude_trans_hotspot':lambda: derive(base)[lambda d: d.hotspot==0],
        'PAV_canonical':lambda: derive(base).assign(PAV=lambda d: d.PAV_canon),
        'PAV_ppsum_ge_0.2':lambda: derive(base).assign(PAV=lambda d: d.PAV_pp02),
        'PAV_lead_variant':lambda: derive(base).assign(PAV=lambda d: d.PAV_lead)}
    for k,f in sc.items():
        M=f(); res[k]=est(M); res[k]['n_cs']=int(len(M)); print(k,res[k]['H2a_Agtex'][:2],res[k]['H5_Aeqtl'][:2],flush=True)
    dump(res,'sensitivity.json')

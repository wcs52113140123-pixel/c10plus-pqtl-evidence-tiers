"""H1 (Sun 2018 internal reproduction) and H2 (deCODE v4 confirmatory)."""
import sys; sys.path.insert(0,'.')
from c10common import *
from scipy.stats import pearsonr
import statsmodels.formula.api as smf
def anchor_block(D,rep,has,label):
    cov=['mlp','logn95','PAV',has]; o={}
    o['n_cs']=int(len(D)); o['n_genes']=int(D.ga.nunique()); o['rep_rate']=float(D[rep].mean()); o['A_gtex_rate']=float(D.A_gtex.mean()); o['A_nongtex_rate']=float(D.A_nongtex.mean()); o['PAV_rate']=float(D.PAV.mean())
    o['main'],_=glm(D,rep,['A_gtex',*cov])
    o['negctrl_wrong'],_=glm(D,rep,['W_gtex',*cov])
    j,mc=glm(D,rep,['A_gtex','W_gtex',*cov],terms=['A_gtex','W_gtex']); o['joint']=j; o['joint_diff']=diff_test(mc,'A_gtex','W_gtex')
    # restrict wrong-gene to CS that have a non-cognate overlapping eQTL CS (more stringent control, descriptive)
    o['nongtex'],_=glm(D,rep,['A_nongtex',*cov])
    S=D[(D.mlog10p>30)&(D[has]==1)]
    o['power_matched_subset'],_=glm(S,rep,['A_gtex','mlp','logn95','PAV'])
    o['PAV_adj'],_=glm(D,rep,['PAV','mlp','logn95',has,'A_gtex'])
    o['four_group']=D.groupby(['A_gtex','PAV'])[rep].agg(['mean','size']).reset_index().to_dict('records')
    o['four_group_Aeqtl']=D.groupby(['A_eqtl','PAV'])[rep].agg(['mean','size']).reset_index().to_dict('records')
    # PAV power matching (H2b method)
    Mt=pav_match(D,rep,has)
    if Mt is not None:
        o['pav_matched_n_sets']=int(Mt.set_id.nunique()); o['pav_matched_n']=int(len(Mt))
        o['pav_matched'],_=glm(Mt,rep,['PAV'],cluster='set_id')
        o['pav_matched_rates']={'PAV1':float(Mt[Mt.PAV==1][rep].mean()),'PAV0':float(Mt[Mt.PAV==0][rep].mean())}
        o['pav_smd_before']=smd(D.dropna(subset=['mlp','absb','maf','logn95']),'PAV'); o['pav_smd_after']=smd(Mt,'PAV')
        o['pav_n_treated_total']=int(D.PAV.sum())
    # has_cis perfectly separates (no cis signal on the gene -> no replication possible), so the pre-specified model with
    # the has_cis term does not converge; the equivalent MLE is obtained on the has_cis==1 subset (documented deviation).
    d=D[D[has]==1].dropna(subset=['mlp','absb','maf','logn95']).copy()
    f=f"{rep} ~ PAV + cr(mlp, df=4) + cr(absb, df=4) + cr(maf, df=4) + cr(logn95, df=4)"
    m=smf.glm(f,d,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(d.ga)[0]})
    o['pav_spline_adj']=dict(OR=float(np.exp(m.params['PAV'])),lo=float(np.exp(m.conf_int().loc['PAV',0])),hi=float(np.exp(m.conf_int().loc['PAV',1])),P=float(m.pvalues['PAV']),n=int(len(d)))
    return o
def eff_corr(D,bu,bo,label):
    E=D.dropna(subset=[bu,bo]).copy(); E=E[(E[bu]!=0)&(E[bo]!=0)]
    E['opp']=(np.sign(E[bu])!=np.sign(E[bo])).astype(int)
    o={'n_pairs':int(len(E))}
    for k,g in E.groupby('PAV'):
        r=pearsonr(g[bu],g[bo])[0] if len(g)>3 else np.nan
        gs=g[g.opp==0]; rs=pearsonr(gs[bu],gs[bo])[0] if len(gs)>3 else np.nan
        o[f'PAV{k}']=dict(n=int(len(g)),r=float(r),r_samedir=float(rs),n_samedir=int(len(gs)),opp_frac=float(g.opp.mean()),n_opp=int(g.opp.sum()))
    if 'PAV0' in o and 'PAV1' in o:
        o['fisher_z_PAV1_lt_PAV0']=fisher_z_diff(o['PAV1']['r'],o['PAV1']['n'],o['PAV0']['r'],o['PAV0']['n'])
        o['fisher_z_samedir_PAV1_lt_PAV0']=fisher_z_diff(o['PAV1']['r_samedir'],o['PAV1']['n_samedir'],o['PAV0']['r_samedir'],o['PAV0']['n_samedir'])
        t=[[o['PAV1']['n_opp'],o['PAV1']['n']-o['PAV1']['n_opp']],[o['PAV0']['n_opp'],o['PAV0']['n']-o['PAV0']['n_opp']]]
        OR,p=fisher_exact(t,alternative='greater'); o['opp_fisher_PAV1_gt']=dict(table=t,OR=float(OR),P_one_sided=float(p))
    # by anchor too (descriptive)
    o['by_anchor_PAV']=E.groupby(['A_eqtl','PAV']).agg(n=('opp','size'),opp_frac=('opp','mean')).reset_index().to_dict('records')
    return o,E
if __name__=="__main__":
    M=derive(load())
    res={}
    H1=M[M.v3].copy(); H2=M[M.v4&~M.v3].copy()
    res['H1']=anchor_block(H1,'rep_sun','sun_has_cis','H1')
    res['H1']['effect_corr'],E1=eff_corr(H1[H1.rep_sun==1],'sun_b_ukb','sun_b_sun','sun')
    res['H2']=anchor_block(H2,'rep_decode','decode_has_cis','H2')
    res['H2']['effect_corr'],E2=eff_corr(H2[(H2.rep_decode==1)|(H2.rep_decode_99==1)],'dec_b_ukb','dec_b_dec','decode')
    # H2 replication-rule sensitivities
    for alt in ['rep_decode_leadonly','rep_decode_99']:
        res['H2']['main_'+alt],_=glm(H2,alt,['A_gtex','mlp','logn95','PAV','decode_has_cis'])
        Mt=pav_match(H2,alt,'decode_has_cis'); res['H2']['pav_matched_'+alt],_=glm(Mt,alt,['PAV'],cluster='set_id')
    # pilot-overlap disclosure: also deCODE on the Sun genes (v3 & v4) as an extra internal check (descriptive)
    H2s=M[M.v4&M.v3].copy(); res['H2_on_Sun_genes_descriptive'],_=glm(H2s,'rep_decode',['A_gtex','mlp','logn95','PAV','decode_has_cis'])
    E1.to_csv(R+'h1_effect_pairs.csv',index=False); E2.to_csv(R+'h2_effect_pairs.csv',index=False)
    dump(res,'h1_h2.json'); print(json.dumps({k:{kk:res[k][kk] for kk in ['n_cs','rep_rate','main','pav_matched'] if kk in res[k]} for k in ['H1','H2']},indent=1,default=str))

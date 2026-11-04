"""H5: cross-ancestry replication of UKB-PPP EUR cis credible sets in pooled non-European UKB-PPP ancestries."""
import sys; sys.path.insert(0,'.')
from c10common import *
def h5_block(M,rep='xa_rep',se='xa_se',expo='A_eqtl'):
    D=M[M[rep].notna()&M[se].notna()].copy(); D['logse']=np.log(D[se])
    cov=['mlp','absb','maf','logse','logn95']
    o={'n_evaluable':int(len(D)),'n_genes':int(D.ga.nunique()),'rep_rate':float(D[rep].mean())}
    if len(D)<500: o['stop_rule_triggered']=True
    else: o['stop_rule_triggered']=False
    o['main'],_=glm(D,rep,[expo,*cov])
    o['rates']=D.groupby([expo,'PAV'])[rep].agg(['mean','size']).reset_index().to_dict('records')
    return o,D
if __name__=="__main__":
    M=derive(load()); res={}
    res['primary'],D=h5_block(M)
    res['negctrl_wrong'],_=h5_block(M,expo='W_gtex')
    j,mc=glm(D.assign(logse=np.log(D.xa_se)),'xa_rep',['A_eqtl','W_gtex','mlp','absb','maf','logse','logn95'],terms=['A_eqtl','W_gtex']); res['joint']=j; res['joint_diff']=diff_test(mc,'A_eqtl','W_gtex')
    res['A_gtex_only'],_=h5_block(M,expo='A_gtex')
    res['samedir'],_=h5_block(M,rep='xa_samedir')
    res['PAV_term'],_=glm(D.assign(logse=np.log(D.xa_se)),'xa_rep',['PAV','A_eqtl','mlp','absb','maf','logse','logn95'])
    # per ancestry (AFR, CSA) with ancestry-specific SE
    IV=pd.read_parquet(W+'nonEUR_index_variants.parquet')
    for anc,tag in [('African','afr'),('Central_South Asian','csa'),('East Asian','eas'),('Middle East','mid'),('American','amr')]:
        e=IV[IV.ancestry==anc].set_index('a')
        Mx=M.set_index('a'); Mx[f'xa_rep_{tag}2']=((np.sign(e.b_h)==np.sign(e.eur_b))&(2*norm.sf(np.abs(e.b_h/e.nse))<0.05)).astype(float).reindex(Mx.index)
        Mx[f'se_{tag}']=e.nse.reindex(Mx.index); Mx=Mx.reset_index()
        res[f'ancestry_{tag}'],_=h5_block(Mx,rep=f'xa_rep_{tag}2',se=f'se_{tag}')
    # effect-consistency slope (per ancestry rows; WLS through origin with A_eqtl interaction, gene-clustered)
    E=IV.merge(M[['a','ga','A_eqtl','PAV']],on='a'); E['eur_b_x_A']=E.eur_b*E.A_eqtl
    w=1/E.nse**2
    m=sm.WLS(E.b_h,E[['eur_b','eur_b_x_A']],weights=w).fit(cov_type='cluster',cov_kwds={'groups':pd.factorize(E.ga)[0]})
    res['slope']=dict(slope_unanchored=float(m.params['eur_b']),slope_anchored=float(m.params['eur_b']+m.params['eur_b_x_A']),interaction=float(m.params['eur_b_x_A']),P_interaction=float(m.pvalues['eur_b_x_A']),n_rows=int(len(E)))
    for k,g in E.groupby('A_eqtl'):
        mm=sm.WLS(g.b_h,g[['eur_b']],weights=1/g.nse**2).fit(); res['slope'][f'slope_A{k}_separate']=float(mm.params['eur_b'])
    # PAV sensitivities
    for k in [0,1]: res[f'PAV{k}_stratum'],_=h5_block(M[M.PAV==k])
    res['excl_mafdiff_gt0.1'],_=h5_block(M[M.xa_maxmafdiff.fillna(0)<=0.1])
    res['lead_only_variant'],_=h5_block(M,rep='xa_rep_leadonly',se='xa_se_leadonly')
    dump(res,'h5.json'); print(json.dumps({k:res[k] for k in ['primary','negctrl_wrong','joint_diff','slope']},indent=1,default=str)[:4000])

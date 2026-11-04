"""Post hoc robustness analyses (Additional file 6).
Input: resource v2/v2.1 parquet (Additional file 3). Requires c10common.py from the analysis code (pav_match, glm, fisher).
Usage: python revision_analyses.py <resource.parquet> <path/to/analysis_code_dir> <out.json>
"""
import sys, json, numpy as np, pandas as pd
from scipy.stats import norm, mannwhitneyu
import statsmodels.api as sm
PQ, CODE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, CODE); import c10common as cc
M = pd.read_parquet(PQ)
lp = M.lead_mlog10p.abs(); M['mlog10p'] = lp; M['mlp'] = np.minimum(lp, 300) / 10
M['absb'] = M.lead_beta.abs(); M['maf'] = np.minimum(M.lead_maf_gnomad_nfe, 1 - M.lead_maf_gnomad_nfe)
M['logn95'] = np.log(np.maximum(M.n_variants_cs95.fillna(1), 1)); M['ga'] = M.gene_id
M['PAV'] = M.pav_cognate_cs95_vep116_alltx.fillna(0).astype(int)
M['A_gtex'] = (M.h4max_eqtl_gtex_v10 >= 0.8).astype(int); M['W_gtex'] = (M.h4max_wronggene_eqtl_gtex >= 0.8).astype(int)
M['A_eqtl'] = M.anchor_eqtl.fillna(0).astype(int); M['A2'] = (M.n_anchor_eqtl_datasets.fillna(0) >= 2).astype(int)
M['v3'] = M.somascan_v3_sun2018_assayed == 1; M['v4'] = M.somascan_v4_decode_assayed == 1
M['n_cs_gene'] = M.groupby('ga').studyLocusId.transform('size'); M['log_ncs'] = np.log(M.n_cs_gene)
lab = {'supported': 1.0, 'refuted': 0.0}
for c in ['ms_primary', 'ms_ext_primary']: M[c + '_y'] = M[c].map(lab)
out = {}
# --- H2 reconstruction: gene-level 'any deCODE cis signal' approximated as >=1 replicated CS in gene
H2 = M[M.v4 & ~M.v3].copy(); H2['rep'] = H2.rep_decode_v4; H2['has'] = H2.groupby('ga').rep.transform('max').astype(int)
cal = {}
for c in [0.2, 0.5, 1.0, 99]:
    Mt = cc.pav_match(H2, 'rep', 'has', caliper=c); o, _ = cc.glm(Mt, 'rep', ['PAV'], cluster='set_id'); cal[str(c)] = dict(n_sets=int(Mt.set_id.nunique()), **o['PAV'])
out['H2b_caliper'] = cal
Mt = cc.pav_match(H2, 'rep', 'has', caliper=0.5); Mt['set_gene'] = Mt.groupby('set_id').ga.transform('first')
a, b = Mt[Mt.PAV == 1].n_cs_gene, Mt[Mt.PAV == 0].n_cs_gene
out['ncs_balance'] = dict(med_PAV1=float(a.median()), med_PAV0=float(b.median()), SMD_log=float((np.log(a).mean() - np.log(b).mean()) / np.sqrt((np.log(a).var() + np.log(b).var()) / 2)), P=float(mannwhitneyu(a, b).pvalue))
out['H2b_cluster_gene'] = cc.glm(Mt, 'rep', ['PAV'], cluster='ga')[0]['PAV']
out['H2b_adj_logncs'] = cc.glm(Mt, 'rep', ['PAV', 'log_ncs'], cluster='set_id')[0]['PAV']
Mt5 = cc.pav_match(H2, 'rep', 'has', covs=('mlp', 'absb', 'maf', 'logn95', 'log_ncs'))
out['H2b_match_ncs'] = dict(n_sets=int(Mt5.set_id.nunique()), **cc.glm(Mt5, 'rep', ['PAV'], cluster='set_id')[0]['PAV'])
rng = np.random.default_rng(20261103); g = {k: v for k, v in Mt.groupby('set_gene')}; keys = list(g); bl = []
for _ in range(1000):
    d = pd.concat([g[keys[j]] for j in rng.choice(len(keys), len(keys), replace=True)])
    bl.append(sm.GLM(d.rep.astype(float), sm.add_constant(d.PAV.astype(float)), family=sm.families.Binomial()).fit().params['PAV'])
out['H2b_gene_bootstrap'] = dict(OR=float(np.exp(np.median(bl))), lo=float(np.exp(np.quantile(bl, .025))), hi=float(np.exp(np.quantile(bl, .975))))
d = H2[H2.has == 1]
out['H2a_recon'] = cc.glm(d, 'rep', ['A_gtex', 'mlp', 'logn95', 'PAV'], cluster='ga')[0]['A_gtex']
out['H2a_adj_ncs'] = cc.glm(d, 'rep', ['A_gtex', 'mlp', 'logn95', 'PAV', 'log_ncs'], cluster='ga')[0]['A_gtex']
out['H2a_anchor_ge2'] = cc.glm(d, 'rep', ['A2', 'mlp', 'logn95', 'PAV'], cluster='ga')[0]['A2']
D5 = M[M.xa_replicated.notna() & M.xa_meta_se.notna()].copy(); D5['logse'] = np.log(D5.xa_meta_se); cov = ['mlp', 'absb', 'maf', 'logse', 'logn95']
out['H5'] = cc.glm(D5, 'xa_replicated', ['A_eqtl', *cov])[0]['A_eqtl']
out['H5_adj_ncs'] = cc.glm(D5, 'xa_replicated', ['A_eqtl', *cov, 'log_ncs'])[0]['A_eqtl']
out['H5_anchor_ge2'] = cc.glm(D5, 'xa_replicated', ['A2', *cov])[0]['A2']
# --- MS layer
def fz(df, gcol, y):
    df = df[df[y].notna()]; r = cc.fisher(df.assign(yy=df[y].astype(int)), gcol, 'yy'); return dict(table=r['table'], OR=r['OR'], P=r['P'])
pep = M.ms_pav_peptide_flag.map({True: 1, False: 0})
out['PAV_MSPA'] = fz(M, 'PAV', 'ms_primary_y'); out['PAV_MSPA_excl_peptide'] = fz(M[pep != 1], 'PAV', 'ms_primary_y'); out['PAV_A1'] = fz(M, 'PAV', 'ms_ext_primary_y')
T = M[M.ms_ext_primary_y.notna()]
out['A1_wronggene_adj'] = cc.glm(T.assign(y=T.ms_ext_primary_y), 'y', ['W_gtex', 'mlp'])[0]['W_gtex']
out['A1_wronggene_leadpergene'] = fz(T.sort_values('lead_mlog10p', ascending=False).drop_duplicates('ga'), 'W_gtex', 'ms_ext_primary_y')
alpha = 0.05 / 1779; zc = norm.isf(alpha / 2); E = M[M.seer_ext_zexp.notna()].copy()
cap = E[E.seer_ext_zexp.abs() >= 6.5]; rr = (cap.seer_ext_beta_aligned / cap.seer_ext_se).abs() / cap.seer_ext_zexp
out['equal_effect_ratio'] = dict(n=int(len(cap)), median=float(rr.median()), q25=float(rr.quantile(.25)), q75=float(rr.quantile(.75)), frac_ge1=float((rr >= 1).mean()))
sens = {}
for s in [0.8, 0.9, 0.95, 1.0, 1.1, 1.25, 1.5]:
    z = E.seer_ext_zexp.abs() * s; pw = norm.cdf(z - zc) + norm.cdf(-z - zc)
    Ms = M.copy(); Ms['lab_s'] = np.nan; Ms.loc[E.index, 'lab_s'] = np.where(E.ms_ext_primary_y == 1, 1.0, np.where(pw >= 0.99, 0.0, np.nan))
    sens[str(s)] = dict(n_refuted=int((Ms.lab_s == 0).sum()), anchor=fz(Ms, 'A_eqtl', 'lab_s') if (Ms.lab_s == 0).sum() >= 10 else None)
out['equal_effect_sensitivity'] = sens
json.dump(out, open(OUT, 'w'), indent=1, default=float)
print(json.dumps({k: (v.get('OR') if isinstance(v, dict) and 'OR' in v else None) for k, v in out.items()}, default=float))

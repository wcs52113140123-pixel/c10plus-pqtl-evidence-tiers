import json;d=json.load(open('../results/a4.json')); s=d['secondary']
def h3s(k):
    h=s[k]; a=h['H3i_logit'].get('A_eqtl',{})
    return dict(n=h['n_evaluable'],sup=h['n_supported'],H3i_OR=a.get('OR'),H3i_P=a.get('P'),H3ii=(h['H3ii_CNP_vs_CA'].get('OR'),h['H3ii_CNP_vs_CA'].get('P'),h['H3ii_CNP_vs_CA'].get('n_g1'),h['H3ii_CNP_vs_CA'].get('n_g0')),PAV=(h['PAV_fisher'].get('OR'),h['PAV_fisher'].get('P'),h['PAV_fisher'].get('rate_g1'),h['PAV_fisher'].get('rate_g0')),conc=(h['concordant_fisher'].get('OR'),h['concordant_fisher'].get('P')))
for k in ['gs_nominal_H3','gs_cswide_H3','gs_excl_beta_flag_H3','gs_multi_gene_H3','gs_new_genes_H3','hb_nominal_H3','hb_cswide_H3','hb_excl_beta_flag_H3','hb_new_genes_H3','combined_H3','S_PAVdir_H3','HB_tier1_H3','HB_tierflag_H3']:
    try: print(k,h3s(k))
    except Exception as e: print(k,'ERR',e)
for k in ['gs_PAV_fisher','hb_PAV_fisher','S_PAVdir_PAV_fisher','HB_tier1_PAV_fisher','HB_tierflag_PAV_fisher','S_PAVdir','HB_tier']: print(k,s[k])
for k in ['gs_tiers','hb_tiers']: print(k,s[k])
for k,v in s['kappa'].items(): print(k,v)
print(d['info_gs_multi'])

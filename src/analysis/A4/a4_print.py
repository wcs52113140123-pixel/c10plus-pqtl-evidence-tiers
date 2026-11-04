import json;d=json.load(open('../results/a4.json'))
for l in ['gs','hb']:
    x=d[l]; h=x['H3']
    print(l,'n_eval',h['n_evaluable'],'sup',h['n_supported'],'genes',h['n_genes'],'stop100',x['layer_stop_lt100'])
    print(' H3i_logit',h['H3i_logit'].get('A_eqtl')); print(' H3i_fisher',h['H3i_fisher']); print(' negctrl',h['H3i_negctrl_wrong'].get('W_gtex'))
    print(' conc groups',h['concordant_groups']); print(' H3ii',h['H3ii_CNP_vs_CA'],'stop',h['H3ii_stop_rule_triggered']); print(' H3ii adj',h['H3ii_logit_adj'].get('CNP'))
    print(' PAV_fisher',h['PAV_fisher']); print(' conc_fisher',h['concordant_fisher'])
    H4=x['H4']; print(' H4',{k:H4[k] for k in H4 if k not in('delong_full_vs_B1',)})
    print(' H4 vs B1',H4.get('delong_full_vs_B1'))
    print(' H4 allchr',x['H4_allchr'].get('delong_full_vs_B0'),x['H4_allchr'].get('n_test'))
    print(' holm',x['holm'])
print(d['verdicts'],d['verdicts_hb'])

"""Build the per-credible-set robustness grading resource and the cardiometabolic drug-target application table."""
import sys; sys.path.insert(0,'.')
from c10common import *
M=derive(load()); S=pd.read_parquet(R+'h4_scores.parquet'); M=M.merge(S,on='a',how='left')
def lab_ms(x): return 'not_assessed' if pd.isna(x) else ('supported' if x==1 else 'refuted')
def dirn(bu,bo,rep):
    if rep!=1: return 'NA'
    if pd.isna(bu) or pd.isna(bo) or bu==0 or bo==0: return 'undetermined'
    return 'same' if np.sign(bu)==np.sign(bo) else 'opposite'
T=pd.DataFrame({
 'studyLocusId':M.a,'studyId':M.studyId,'gene_id':M.ga,'symbol':M.symbol,'olink_oid':M.oid,'chrom':M.chrom,'lead_variant_id':M.lead_variant_id,'pos_b38':M.pos_b38,
 'lead_beta':M.beta,'lead_mlog10p':M.mlog10p,'lead_maf_gnomad_nfe':M.maf,'n_variants_cs95':M.n95,'purity_min_r2':M.purityMinR2,'dist_lead_to_tss':M.dist_tss,
 'h4max_eqtl_gtex_v10':M.gtex_h4,'h4max_eqtl_eqtlcatalogue':M.nongtex_h4,'h4max_sqtl_tuqtl':M.iso_h4,'h4max_wronggene_eqtl_gtex':M.wrong_h4,
 'n_anchor_eqtl_datasets':M.n_anchor_datasets.astype(int),'anchor_eqtl':M.A_eqtl,'anchor_eqtl_gtex':M.A_gtex,'anchor_eqtl_nongtex':M.A_nongtex,'anchor_isoform':M.A_iso,
 'pav_cognate_cs95_vep116_alltx':M.PAV,'pav_pip_sum':M.pav_pp,'pav_cognate_cs95_canonical':M.PAV_canon,'pav_lead_variant':M.PAV_lead,
 'somascan_v3_sun2018_assayed':M.v3.astype(int),'h4max_sun2018':np.where(M.v3,M.sun_h4,np.nan),'rep_sun2018':np.where(M.v3,M.rep_sun,np.nan),
 'dir_sun2018':[dirn(a,b,r) for a,b,r in zip(M.sun_b_ukb,M.sun_b_sun,np.where(M.v3,M.rep_sun,0))],
 'somascan_v4_decode_assayed':M.v4.astype(int),'rep_decode_v4':np.where(M.v4,M.rep_decode,np.nan),
 'dir_decode_v4':[dirn(a,b,r) for a,b,r in zip(M.dec_b_ukb,M.dec_b_dec,np.where(M.v4,M.rep_decode,0))],
 'platform_concordance':M.conc_cat,
 'ms_mspa_source':M.ms_sheet,'ms_primary':[lab_ms(x) for x in M.ms_primary],'ms_mspa_peptide':[lab_ms(x) for x in M.ms_mspa],'ms_pav_peptide_flag':M.ms_pep_pav,
 'ms_seer_discovery':np.where(M.ms_seer.isna(),'not_assessed',np.where(M.ms_seer==1,'overlap','no_overlap')),
 'xa_n_ancestries':M.xa_nanc,'xa_meta_beta_aligned':M.xa_b,'xa_meta_se':M.xa_se,'xa_meta_p':M.xa_p,'xa_replicated':M.xa_rep,'xa_same_direction':M.xa_samedir,
 'in_hla':M.hla,'n_other_genes_trans_cs_shared':M.n_trans_genes_shared.astype(int),
 'tier':M.tier,'p_ms_support_model':M.p_ms_support_model,'p_xa_rep_model':M.p_xa_rep_model})
T['tier_label']=T.tier.map({1:'T1_transcript_anchored_noPAV',2:'T2_anchored_withPAV_or_isoform_only',3:'T3_unresolved_noAnchor_noPAV',4:'T4_caution_noAnchor_withPAV'})
q=T.p_ms_support_model.quantile([1/3,2/3]).values
T['model_confidence']=np.where(T.p_ms_support_model>=q[1],'upper_tertile',np.where(T.p_ms_support_model>=q[0],'middle_tertile','lower_tertile'))
T['flag_concordant_unanchored_PAV']=((T.platform_concordance=='concordant')&(T.anchor_eqtl==0)&(T.pav_cognate_cs95_vep116_alltx==1)).astype(int)
T.to_parquet(R+'C10plus_cis_pqtl_robustness_grades_v1.parquet',index=False)
T.to_csv(R+'C10plus_cis_pqtl_robustness_grades_v1.tsv.gz',sep='\t',index=False,compression='gzip')
GENES=['PCSK9','HMGCR','NPC1L1','APOB','LPA','ANGPTL3','ANGPTL4','CETP','APOC3','LPL','LDLR','IL6R','IL6','IL1B','F2','F10','F11','ACE','AGT','REN','DPP4','GCG','LEP','FGF21','GDF15','PLA2G7','MPO','CRP','ASGR1','APOA1','APOE','LPA','PCSK9','ANGPTL8','IL1RN','TNF','CD40','SELE','SELP','MMP12','LGALS3','NPPB','GP6','KLKB1','F12','PLG','SERPINE1','CCL2','IL18','ADIPOQ','LEPR','SORT1','LIPC','LCAT','PLTP','APOA5','GHR','IGF1','TFRC','HP']
A=T[T.symbol.isin(set(GENES))].sort_values(['symbol','lead_mlog10p'],ascending=[True,False])
A.to_csv(R+'application_cardiometabolic_targets.csv',index=False)
g=A.groupby('symbol').agg(n_cs=('tier','size'),n_T1=('tier',lambda s:(s==1).sum()),n_T2=('tier',lambda s:(s==2).sum()),n_T3=('tier',lambda s:(s==3).sum()),n_T4=('tier',lambda s:(s==4).sum()),
    top_cs_tier=('tier','first'),top_cs_mlog10p=('lead_mlog10p','first'),top_cs_PAV=('pav_cognate_cs95_vep116_alltx','first'),top_cs_concordance=('platform_concordance','first'),top_cs_ms=('ms_primary','first'),top_cs_xa=('xa_replicated','first'))
g.to_csv(R+'application_cardiometabolic_summary.csv')
print(len(T), T.tier.value_counts().sort_index().to_dict(), 'flag',T.flag_concordant_unanchored_PAV.sum(), 'app genes',A.symbol.nunique(), len(A))
print(g.to_string())

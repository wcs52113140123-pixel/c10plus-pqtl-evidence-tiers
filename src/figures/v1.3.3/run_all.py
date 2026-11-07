"""Regenerate C10+ v1.3.3 main Figs 1-6 and Supplementary Figs S1-S8, their source data and the alignment report.
Inputs: results workbook (Additional file 5, v6), resource v2.3 TSV, Additional file 2 (figure source sheets used as inputs).
Usage: python run_all.py --config figure_inputs.json   (see figure_inputs.example.json), or --workbook/--resource/--af2."""
import sys, os, re, json, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import common, fig1v3 as fig1, fig2, fig3, fig4, fig5, fig6, sfigs_a, sfigs_b
def run(WB_PATH,RES_PATH,AF2_PATH):
    WB=pd.read_excel(WB_PATH,sheet_name=None); X=pd.read_excel(AF2_PATH,sheet_name=None)
    R=pd.read_csv(RES_PATH,sep='\t',low_memory=False); R['pav']=R.pav_cognate_cs95_vep116_alltx.astype(float)
    figs={}; src={}; rep={}
    F,r,s=fig1.make(WB); figs['Fig1']=F; rep['Fig1']=r
    old=X['Table S23']; s1a=old[old.panel_sheet=='Fig1a'].dropna(axis=1,how='all').copy()
    src['Table S23']=pd.concat([s1a,s.rename(columns={'panel':'panel_sheet'})],ignore_index=True)
    for nm,call,key in [('Fig2',lambda: fig2.make(R,WB,X['Table S24']),'Table S24'),('Fig3',lambda: fig3.make(WB,X['Table S25']),'Table S25'),('Fig4',lambda: fig4.make(R),'Table S26'),
                        ('Fig5',lambda: fig5.make(R,WB,X['Table S27']),'Table S27'),('Fig6',lambda: fig6.make(X['Table S28'],WB),'Table S28'),
                        ('FigS1',lambda: sfigs_a.S1(X['Table S29']),'Table S29'),('FigS2',lambda: sfigs_a.S2(X['Table S30'],R),'Table S30'),('FigS3',lambda: sfigs_a.S3(R),'Table S31'),
                        ('FigS4',lambda: sfigs_a.S4(X['Table S32']),'Table S32'),('FigS5',lambda: sfigs_b.S5(X['Table S33']),'Table S33'),('FigS6',lambda: sfigs_b.S6(R),'Table S34'),
                        ('FigS7',lambda: sfigs_b.S7(X['Table S35']),'Table S35'),('FigS8',lambda: sfigs_b.S8(X['Table S36'],R),'Table S36')]:
        F,r,s=call(); figs[nm]=F; rep[nm]=r; src[key]=s.rename(columns={'panel':'panel_sheet'})
    for k in src: 
        cols=['panel_sheet']+[c for c in src[k].columns if c!='panel_sheet']; src[k]=src[k][cols]
    return figs,src,rep,WB,X,R
def fig_texts(F):
    r=F.fig.canvas.get_renderer(); out=[]
    hidden=set()
    for ax in F.fig.axes:
        for axis,lim in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            lo,hi=min(lim),max(lim)
            for tk in axis.get_major_ticks()+axis.get_minor_ticks():
                if not (lo-1e-9<=tk.get_loc()<=hi+1e-9): hidden.add(tk.label1)
    ticks=set()
    for ax in F.fig.axes: ticks|=set(ax.get_xticklabels()+ax.get_yticklabels())
    for t in F.fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text().strip() and t not in hidden: out.append((t.get_text(), t in ticks))
    return out


def main(argv=None):
    """Command-line entry point. Input paths come from --config (JSON) and/or explicit arguments; nothing is hard-coded."""
    import argparse
    ap = argparse.ArgumentParser(description="Regenerate the v1.3.3 main figures (Fig. 1-6) and Supplementary Figs S1-S8 into ./out/")
    ap.add_argument("--config", help="JSON file with keys 'workbook', 'resource', 'af2'; relative paths are resolved against the config file's folder")
    ap.add_argument("--workbook", help="results workbook (Additional file 5, v6)")
    ap.add_argument("--resource", help="resource TSV(.gz) v2.3 (resource/C10plus_cis_pqtl_robustness_grades_v2.3.tsv.gz)")
    ap.add_argument("--af2", help="Additional file 2 workbook (sheets Table S23-S36 are read as figure inputs)")
    a = ap.parse_args(argv)
    cfg = {}
    if a.config:
        base = os.path.dirname(os.path.abspath(a.config))
        cfg = {k: (v if os.path.isabs(v) else os.path.join(base, v)) for k, v in json.load(open(a.config)).items()}
    for k in ("workbook", "resource", "af2"):
        if getattr(a, k): cfg[k] = getattr(a, k)
    miss = [k for k in ("workbook", "resource", "af2") if not cfg.get(k)]
    if miss: ap.error("missing input path(s): " + ", ".join(miss))
    figs, src, rep, WB, X, R = run(cfg["workbook"], cfg["resource"], cfg["af2"])
    json.dump(rep, open(os.path.join(common.OUT, "figure_layout_report.json"), "w"), indent=1, default=str)
    print("wrote", len(figs), "figures to", os.path.relpath(common.OUT))

if __name__ == "__main__":
    main()

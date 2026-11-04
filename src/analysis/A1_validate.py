import gzip, glob, os, sys, pandas as pd, json
SE="./seer_ext/"
ST=pd.read_csv(SE+'seer_study_windows.csv',dtype={'chrom':str}).set_index('gcst')
bad={}; ok=0
for f in sorted(glob.glob(SE+'all/*.cis.tsv.gz')):
    g=os.path.basename(f).split('.')[0]
    try:
        with gzip.open(f,'rt') as h: lines=h.read().split('\n')
        if lines[-1]=='': lines=lines[:-1]
        hdr=lines[0].split('\t'); assert len(hdr)==12 and hdr[0]=='chromosome', 'header'
        rows=[l.split('\t') for l in lines[1:]]
        assert all(len(r)==12 for r in rows), 'ncol'
        ch=ST.loc[g,'chrom']; w0=ST.loc[g,'w0_37']; w1=ST.loc[g,'w1_37']
        pos=[int(r[1]) for r in rows]
        assert all(r[0]==ch for r in rows), 'chrom'
        assert all(w0<=p<=w1 for p in pos), 'window'
        assert pos==sorted(pos), 'unsorted'
        assert len(rows)>0, 'empty'
        ok+=1
    except Exception as e: bad[g]=str(e)[:80]
json.dump(bad,open(SE+'validation_bad.json','w')); print(ok, len(bad), list(bad.items())[:10])

#!/bin/bash
# POST HOC (L1): copy frozen analysis scripts into posthoc_L1L2/, redirect outputs, add AUTOSOME_ONLY switch to load()
set -e
W=.; PH=$W/posthoc_L1L2
mkdir -p $PH/results_full $PH/results_auto
for f in c10common.py h1_h2.py h5.py h3.py h4.py a2_analysis.py; do cp $W/$f $PH/$f; done
cd $PH
sed -i "s#sys.path.insert(0,'.')#sys.path.insert(0,'$PH')#" h1_h2.py h5.py h3.py h4.py a2_analysis.py
python3 - <<'PY'
p='c10common.py'; s=open(p).read()
s=s.replace('R="./results/"','import os\nR="./posthoc_L1L2/results_"+os.environ.get("PH_TAG","full")+"/"')
s=s.replace("    return M\ndef derive","    if os.environ.get('AUTOSOME_ONLY')=='1': M=M[M.chrom!='X'].copy()   # POST HOC L1\n    return M\ndef derive",1)
open(p,'w').write(s)
p='a2_analysis.py'; s=open(p).read(); s=s.replace("W+'a2_cs_fields.parquet'","R+'a2_cs_fields.parquet'"); open(p,'w').write(s)
PY
grep -n "AUTOSOME_ONLY\|PH_TAG" c10common.py; grep -n "a2_cs_fields" a2_analysis.py; grep -n "sys.path" h1_h2.py

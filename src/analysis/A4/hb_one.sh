#!/bin/bash
gcst=$1; chrom=$4; w0=$5; w1=$6; url=$7
for t in 1 2 3; do
  (curl -s -r 0-100000 "$url" | gzip -dc 2>/dev/null | head -1; tabix "$url" ${chrom}:${w0}-${w1}) > hb_cis/$gcst.cis.tsv.tmp 2> hb_cis/$gcst.err && [ -s hb_cis/$gcst.cis.tsv.tmp ] && [ $(wc -l < hb_cis/$gcst.cis.tsv.tmp) -gt 1 ] && { gzip -c hb_cis/$gcst.cis.tsv.tmp > hb_cis/$gcst.cis.tsv.gz; echo -e "$gcst\tok\t$(($(wc -l < hb_cis/$gcst.cis.tsv.tmp)-1))\t$(sha256sum hb_cis/$gcst.cis.tsv.gz|cut -c1-64)"; rm hb_cis/$gcst.cis.tsv.tmp; exit 0; }
  sleep 5
done
echo -e "$gcst\tfailed"

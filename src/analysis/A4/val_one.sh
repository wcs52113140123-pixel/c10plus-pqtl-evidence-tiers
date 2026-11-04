#!/bin/bash
gcst=$1; chrom=$4; w0=$5; w1=$6; url=$7
curl -s -o hb_val/$gcst.h.tsv.gz "$url"; md5=$(md5sum hb_val/$gcst.h.tsv.gz | cut -c1-32); exp=$(curl -s "${url%/*}/md5sum.txt" | grep " $gcst.h.tsv.gz$" | cut -c1-32)
a=$(zcat hb_val/$gcst.h.tsv.gz | awk -F'\t' -v c=$chrom -v a=$w0 -v b=$w1 'NR==1 || ($1==c && $2>=a && $2<=b)' | tr -d '\r' | md5sum | cut -c1-32)
b=$(zcat hb_cis/$gcst.cis.tsv.gz | tr -d '\r' | md5sum | cut -c1-32)
echo -e "$gcst\tfile_md5=$md5\texpected=$exp\tmd5_ok=$([ "$md5" == "$exp" ] && echo yes || echo no)\twindow_identical=$([ "$a" == "$b" ] && echo yes || echo no)"
rm -f hb_val/$gcst.h.tsv.gz

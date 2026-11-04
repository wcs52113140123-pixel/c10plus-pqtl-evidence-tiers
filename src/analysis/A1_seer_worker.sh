#!/bin/bash
# Seer full-sumstats streaming cis-window extractor (Addendum A1). Usage: seer_worker.sh <tasks.tsv> <outdir> <parallel>
# task line: gcst url chrom w0_37 w1_37 ; keeps header + rows with chromosome==chrom (or 23 for X) and w0<=pos<=w1.
TASKS=$1; OUT=$2; PAR=$3; mkdir -p $OUT/done $OUT/fail
one() {
  set -o pipefail
  g=$1; u=$2; ch=$3; s=$4; e=$5; O=$OUT
  [ -f $O/done/$g.ok ] && return 0
  d=$(dirname $u)
  exp=$(curl -s -m 120 --retry 5 $d/md5sum.txt | awk -v f="$g.tsv.gz" '$2==f||$2=="./"f{print $1}')
  for att in 1 2 3 4 5; do
    curl -sS --fail -m 14400 --retry 3 --speed-time 300 --speed-limit 1000 "$u" | tee >(md5sum | cut -d' ' -f1 > $O/$g.md5) | gunzip -c | \
      mawk -F'\t' -v c="$ch" -v s="$s" -v e="$e" 'NR==1{print;next} ($1==c || (c=="X" && ($1=="23"||$1=="X"))) && $2>=s && $2<=e' | gzip -6 > $O/$g.cis.tsv.gz.tmp
    rc=$?; sleep 1; got=$(cat $O/$g.md5 2>/dev/null)
    if [ $rc -eq 0 ] && [ -n "$got" ] && { [ -z "$exp" ] || [ "$got" == "$exp" ]; }; then
      mv $O/$g.cis.tsv.gz.tmp $O/$g.cis.tsv.gz
      n=$(gunzip -c $O/$g.cis.tsv.gz | wc -l); sh=$(sha256sum $O/$g.cis.tsv.gz | cut -d' ' -f1)
      echo -e "$g\t$got\t$exp\t$n\t$sh\t$att\t$(date -u +%FT%TZ)" > $O/done/$g.ok; rm -f $O/$g.md5; return 0
    fi
    echo -e "$g\tattempt$att\trc=$rc\tgot=$got\texp=$exp" >> $O/fail/$g.log; sleep $((att*30))
  done
  return 1
}
export -f one; export OUT
awk -F'\t' '{print $1, $2, $3, $4, $5}' $TASKS | xargs -P $PAR -n 5 bash -c 'one "$@"' _
echo ALLDONE $(ls $OUT/done | wc -l) > $OUT/ALLDONE

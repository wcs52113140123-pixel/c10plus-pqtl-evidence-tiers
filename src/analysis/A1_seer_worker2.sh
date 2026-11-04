#!/bin/bash
# v2: early-stop streaming extractor (files verified sorted by chromosome then position). Usage: seer_worker2.sh <tasks.tsv> <outdir> <parallel>
TASKS=$1; OUT=$2; PAR=$3; mkdir -p $OUT/done $OUT/fail
one() {
  g=$1; u=$2; ch=$3; s=$4; e=$5; O=$OUT
  [ -f $O/done/$g.ok ] && return 0
  d=$(dirname $u)
  exp=$(curl -s -m 120 --retry 5 $d/md5sum.txt | awk -v f="$g.tsv.gz" '$2==f||$2=="./"f{print $1}')
  for att in $(seq 1 30); do
    rm -f $O/$g.mark $O/$g.md5
    curl -sS --fail -m 14400 --speed-time 300 --speed-limit 1000 "$u" 2>/dev/null | tee >(md5sum | cut -d' ' -f1 > $O/$g.md5) | gunzip -c 2>/dev/null | \
      mawk -F'\t' -v c="$ch" -v s="$s" -v e="$e" -v mk="$O/$g.mark" '
        NR==1{print;next}
        { isc = ($1==c || (c=="X" && ($1=="23"||$1=="X"))) }
        isc { seen=1; if ($2+0>e) {stop=1; exit} if ($2+0>=s) print; next }
        seen { stop=1; exit }
        END { print (stop?"early":"eof") > mk }' | gzip -6 > $O/$g.cis.tsv.gz.tmp
    rcs=("${PIPESTATUS[@]}"); sleep 1
    mark=$(cat $O/$g.mark 2>/dev/null); got=$(cat $O/$g.md5 2>/dev/null)
    ok=0; meth=""
    if [ "$mark" == "early" ]; then ok=1; meth="early_stop"; got="NA"; fi
    if [ "$mark" == "eof" ] && [ "${rcs[0]}" == "0" ] && [ "${rcs[2]}" == "0" ] && [ -n "$got" ] && { [ -z "$exp" ] || [ "$got" == "$exp" ]; }; then ok=1; meth="full_md5"; fi
    if [ $ok -eq 1 ]; then
      mv $O/$g.cis.tsv.gz.tmp $O/$g.cis.tsv.gz
      n=$(gunzip -c $O/$g.cis.tsv.gz | wc -l); sh=$(sha256sum $O/$g.cis.tsv.gz | cut -d' ' -f1)
      echo -e "$g\t$got\t$exp\t$n\t$sh\t$att\t$(date -u +%FT%TZ)\t$meth" > $O/done/$g.ok; rm -f $O/$g.md5 $O/$g.mark; return 0
    fi
    echo -e "$g\tattempt$att\trcs=${rcs[*]}\tmark=$mark\tgot=$got\texp=$exp" >> $O/fail/$g.log; sleep $(( att<6 ? att*30 : 180 ))
  done
  return 1
}
export -f one; export OUT
awk -F'\t' '{print $1, $2, $3, $4, $5}' $TASKS | xargs -P $PAR -n 5 bash -c 'one "$@"' _
echo ALLDONE $(ls $OUT/done | wc -l) > $OUT/ALLDONE2

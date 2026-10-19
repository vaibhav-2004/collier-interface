#!/bin/bash
EXE=${1:-./loopdrv}
echo "checking $EXE against reference.txt"
while read -r p2 m02 m12 mu2 re im rest; do
  case "$p2" in \#*|"") continue;; esac
  out=$($EXE B0 $p2 $m02 $m12 $mu2 2>/dev/null | grep -E '^ *[-+0-9]' | tail -1)
  [ -z "$out" ] && { printf "**  no output for p2=%s\n" "$p2"; continue; }
  set -- $out
  awk -v a="$re" -v b="$1" -v c="$im" -v d="$2" -v e1="$3" -v f1="$5" \
      -v p="$p2" -v m="$m02" -v n="$m12" 'BEGIN{
    finok = ((a-b)^2 + (c-d)^2 < 1e-12);
    polok = (e1-1)^2 < 1e-16 && f1*f1 < 1e-16;
    printf "%s p2=%-10s m02=%-12s fin %+12.8f %+12.8f   1/e %s\n",
           (finok && polok ? "OK " : "** "), p, m, b, d,
           (polok ? "= 1 (correct)" : "= " e1 " ** expected 1") }'
done < reference.txt

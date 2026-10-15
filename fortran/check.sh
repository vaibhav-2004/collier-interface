#!/bin/bash
EXE=${1:-./b0ref}
echo "checking $EXE against reference.txt"
while read p2 m02 m12 mu2 re im; do
  out=$($EXE $p2 $m02 $m12 $mu2 2>/dev/null | grep -E '^ *[-+0-9]' | tail -1)
  gr=$(echo $out | awk '{print $1}'); gi=$(echo $out | awk '{print $2}')
  awk -v a=$re -v b=$gr -v c=$im -v d=$gi -v p=$p2 -v m=$m02 -v n=$m12 \
    'BEGIN{printf "%s p2=%-9s m02=%-11s m12=%-11s ref %12.8f %12.8f   got %12.8f %12.8f\n",
     (((a-b)^2+(c-d)^2 < 1e-12) ? "OK " : "** "), p, m, n, a, c, b, d}'
done < reference.txt

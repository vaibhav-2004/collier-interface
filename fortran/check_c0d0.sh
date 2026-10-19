#!/bin/bash
EXE=${1:-./loopdrv}
TOL=${2:-1e-6}
echo "checking $EXE against reference_c0d0.txt  (rel tol $TOL on finite part)"
grep -v '^#' reference_c0d0.txt | grep -v '^$' | while read -r typ rest; do
  n=$(echo $rest | wc -w)
  args=$(echo $rest | cut -d' ' -f1-$((n-2)))
  re=$(echo $rest | awk '{print $(NF-1)}')
  im=$(echo $rest | awk '{print $NF}')
  out=$($EXE $typ $args 2>/dev/null | grep -E '^ *[-+0-9]' | tail -1)
  [ -z "$out" ] && { printf "**  %-3s  no output\n" "$typ"; continue; }
  set -- $out
  awk -v a="$re" -v b="$1" -v c="$im" -v d="$2" \
      -v e1="$3" -v e2="$4" -v f1="$5" -v f2="$6" \
      -v t="$TOL" -v ty="$typ" -v ar="$(echo $args | cut -c1-40)" 'BEGIN{
    num=sqrt((a-b)^2+(c-d)^2); den=sqrt(a*a+c*c);
    rel=(den>0)?num/den:num;
    pol=(e1*e1+e2*e2+f1*f1+f2*f2 > 1e-30) ? "  [has poles]" : "  [finite]";
    printf "%s %-3s %-42s ref %+.8e  got %+.8e%s\n",
           (rel<t?"OK ":"** "), ty, ar, a, b, pol }'
done

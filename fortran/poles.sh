#!/bin/bash
EXE=${1:-./loopdrv}
MU=29756.25
S=1000000;  T=-391997.844;  U=-548131.156;  M=29929
show() { printf "%-32s " "$1"; shift
         $EXE "$@" $MU 2>/dev/null | grep -E '^ *[-+0-9]' | tail -1 | \
         awk '{printf "fin %+11.4e %+11.4e   1/e %+11.4e %+11.4e   1/e2 %+11.4e %+11.4e\n",$1,$2,$3,$4,$5,$6}'; }

echo "=== A0, B0 (UV divergent) ==="
show "A0(mt^2)"                A0 $M
show "B0(s; 0,0)"              B0 $S 0 0
show "B0(s; mt,mt)"            B0 $S $M $M
show "B0(t; 0,mt)"             B0 $T 0 $M
show "B0(u; 0,mt)"             B0 $U 0 $M
show "B0(mt^2; 0,mt)"          B0 $M 0 $M
show "B0(0; mt,mt)"            B0 0 $M $M
echo
echo "=== C0 (vertex corrections) ==="
show "C0(0,0,s; 0,0,0)"        C0 0 0 $S 0 0 0
show "C0(0,0,s; mt,mt,mt)"     C0 0 0 $S $M $M $M
show "C0(mt2,mt2,s; mt,0,mt)"  C0 $M $M $S $M 0 $M
show "C0(0,mt2,t; 0,0,mt)"     C0 0 $M $T 0 0 $M
show "C0(0,mt2,u; 0,0,mt)"     C0 0 $M $U 0 0 $M
echo
echo "=== D0 (boxes) ==="
show "D0 t-ch, 2 massless"     D0 0 0 $M $M $S $T 0 0 $M $M
show "D0 u-ch, 2 massless"     D0 0 0 $M $M $S $U 0 0 $M $M
show "D0 t-ch, 1 massless"     D0 0 0 $M $M $S $T 0 $M $M $M

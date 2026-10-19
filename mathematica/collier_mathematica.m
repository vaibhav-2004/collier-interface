(* ::Package:: *)

(* ==================================================================== *)
(*  COLLIER <-> MATHEMATICA INTERFACE                                   *)
(*                                                                      *)
(*  Wraps the loopdrv executable built in ~/collier_work.               *)
(*  Provides A0, B0, C0, D0 with memoisation.                          *)
(*                                                                      *)
(*  Conventions (set inside loopdrv.f90):                               *)
(*    SetDeltaUV_cll(0)      -> UV pole dropped, finite part returned   *)
(*    SetDeltaIR_cll(0,0)    -> IR poles dropped likewise               *)
(*    SetMuUV2 = SetMuIR2 = mu2                                         *)
(*                                                                      *)
(*  NOTE ON mu:  B0 is UV divergent and D0 is IR divergent, so BOTH     *)
(*  finite parts depend on mu.  C0 (with a massive internal line) is    *)
(*  finite in both senses and genuinely mu-independent.  Always state   *)
(*  mu alongside any number you quote.                                  *)
(* ==================================================================== *)

Clear[colExe, colWSL, colRun, colParse, finOf, c1Of, c2Of,
      A0Full, B0Full, C0Full, D0Full, A0cF, B0cF, C0cF, D0cF,
      A0col, B0col, C0col, D0col,
      A0c, B0c, C0c, D0c, colClearCache,
      colValidate, colRefs, colSelfTest];


(* ==================================================================== *)
(* 1.  SETUP  --  EDIT THESE TWO LINES                                  *)
(* ==================================================================== *)

colExe = "/home/vaibhav/collier_work/loopdrv";

(* True  if Mathematica runs on WINDOWS and the driver lives in WSL
   False if Mathematica runs inside WSL (or native Linux)              *)
colWSL = True;

Print["driver : ", colExe];
Print["via WSL: ", colWSL];
Print["OS     : ", $OperatingSystem];


(* ==================================================================== *)
(* 2.  LOW-LEVEL CALL                                                   *)
(*                                                                      *)
(*  Collier prints a startup banner to stdout, so we keep only lines    *)
(*  that start with a number/sign and take the last one.                *)
(* ==================================================================== *)

(* The driver prints SIX numbers:                                       *)
(*    Re(fin) Im(fin)  Re(c1) Im(c1)  Re(c2) Im(c2)                     *)
(* where  I = c2/eps^2 + c1/eps + fin.                                  *)
(* colParse returns an association so nothing is silently discarded.    *)
(* Collier also prints a startup banner, hence the numeric-line filter. *)

colParse[out_String] := Module[{lines, good, nums},
  lines = StringSplit[out, "\n"];
  good = Select[lines, StringMatchQ[#, RegularExpression[" *[-+0-9].*"]] &];
  If[good === {},
     Print["** no numeric line in driver output:"]; Print[out];
     Return[$Failed]];
  nums = ToExpression /@ StringSplit[
           StringReplace[Last[good], {"E" -> "*^", "D" -> "*^"}]];
  If[Length[nums] =!= 6,
     Print["** expected 6 numbers, got ", Length[nums], ": ", Last[good]];
     Return[$Failed]];
  <|"fin" -> nums[[1]] + I nums[[2]],
    "c1"  -> nums[[3]] + I nums[[4]],
    "c2"  -> nums[[5]] + I nums[[6]]|>
];

(* Convenience: pull one piece out of a result, or $Failed through.     *)
finOf[r_] := If[r === $Failed, $Failed, r["fin"]];
c1Of[r_]  := If[r === $Failed, $Failed, r["c1"]];
c2Of[r_]  := If[r === $Failed, $Failed, r["c2"]];

colRun[what_String, args_List] := Module[{cmd, res},
  cmd = Join[
     If[colWSL, {"wsl", "-d", "Ubuntu-22.04", "-u", "vaibhav", "-e", colExe},
                {colExe}],
     {what},
     ToString[N[#], CForm] & /@ args];
  res = RunProcess[cmd];
  If[res["ExitCode"] =!= 0,
     Print["** driver exit ", res["ExitCode"], ": ", res["StandardError"]];
     Return[$Failed]];
  colParse[res["StandardOutput"]]
];


(* ==================================================================== *)
(* 3.  THE FOUR INTEGRALS                                               *)
(*                                                                      *)
(*  Argument conventions, read off Collier's own demo.f90:              *)
(*    pij  means  (p_i - p_j)^2                                         *)
(*                                                                      *)
(*    A0(m02)                                                           *)
(*    B0(p10; m02, m12)                                                 *)
(*    C0(p10, p21, p20; m02, m12, m22)                                  *)
(*    D0(p10, p21, p32, p30, p20, p31; m02, m12, m22, m32)              *)
(*                                                                      *)
(*  For D0 the first FOUR momenta are the external legs and the last    *)
(*  TWO are the diagonals (the Mandelstams).  Easy to get wrong.        *)
(* ==================================================================== *)

(* *Full  -> association with "fin", "c1", "c2"                        *)
(* plain  -> the finite part only, which is what most code wants        *)

A0Full[m02_, mu2_] := colRun["A0", {m02, mu2}];
B0Full[p10_, m02_, m12_, mu2_] := colRun["B0", {p10, m02, m12, mu2}];
C0Full[p10_, p21_, p20_, m02_, m12_, m22_, mu2_] :=
   colRun["C0", {p10, p21, p20, m02, m12, m22, mu2}];
D0Full[p10_, p21_, p32_, p30_, p20_, p31_,
       m02_, m12_, m22_, m32_, mu2_] :=
   colRun["D0", {p10, p21, p32, p30, p20, p31, m02, m12, m22, m32, mu2}];

A0col[a__] := finOf[A0Full[a]];
B0col[a__] := finOf[B0Full[a]];
C0col[a__] := finOf[C0Full[a]];
D0col[a__] := finOf[D0Full[a]];


(* ==================================================================== *)
(* 4.  MEMOISED VERSIONS  --  USE THESE                                 *)
(*                                                                      *)
(*  WHY THIS MATTERS MORE THAN ANYTHING ELSE HERE.                      *)
(*                                                                      *)
(*  Each scalar integral depends only on the PHASE-SPACE POINT, not on  *)
(*  the helicity configuration.  B0(t,0,mt^2) is the same number for    *)
(*  all sixteen helicity amplitudes.  Without caching you evaluate it   *)
(*  sixteen times, each spawning a process.                             *)
(*                                                                      *)
(*  This is exactly the mistake that took the LO notebook from a cloud  *)
(*  timeout to 1.6 seconds.  It is available again here, at larger      *)
(*  scale.                                                              *)
(* ==================================================================== *)

(* memoised: the Full versions cache, the plain ones read from them,    *)
(* so a point is evaluated once however you ask for it                  *)
A0cF[a__] := A0cF[a] = A0Full[a];
B0cF[a__] := B0cF[a] = B0Full[a];
C0cF[a__] := C0cF[a] = C0Full[a];
D0cF[a__] := D0cF[a] = D0Full[a];

A0c[a__] := finOf[A0cF[a]];
B0c[a__] := finOf[B0cF[a]];
C0c[a__] := finOf[C0cF[a]];
D0c[a__] := finOf[D0cF[a]];

colClearCache[] := (
   Clear[A0cF, B0cF, C0cF, D0cF];
   A0cF[a__] := A0cF[a] = A0Full[a];
   B0cF[a__] := B0cF[a] = B0Full[a];
   C0cF[a__] := C0cF[a] = C0Full[a];
   D0cF[a__] := D0cF[a] = D0Full[a];
   Print["cache cleared"]);


(* ==================================================================== *)
(* 5.  VALIDATION  --  the 15 B0 points                                 *)
(* ==================================================================== *)

Module[{mH = 125.25, mt = 172.5},
 colRefs = {
  {"p2=62500,  mH2, mt2   <-- TASK",  62500., mH^2, mt^2,  0.992443465,  0.},
  {"p2=250,    mH2, mt2   alt",         250., mH^2, mt^2,  0.288037486,  0.},
  {"p2=62500,  mt2, mH2   symmetry",  62500., mt^2, mH^2,  0.992443465,  0.},
  {"p2=1e5,    mH2, mt2   uneq above",100000., mH^2, mt^2, 2.038711720,  1.046280287},
  {"p2=5e5,    mH2, mt2   uneq above",500000., mH^2, mt^2,-0.443817841,  2.843126532},
  {"p2=0,      mH2, mt2",                  0., mH^2, mt^2, 0.286164934,  0.},
  {"p2=mt2,    0,   mt2   on shell",    mt^2,   0.,  mt^2, 2.000000000,  0.},
  {"p2=0,      mt2, mt2",                  0., mt^2, mt^2, 0.,           0.},
  {"p2=5e4,    mt2, mt2   eq below",   50000., mt^2, mt^2, 0.343010838,  0.},
  {"p2=2e5,    mt2, mt2   eq above",  200000., mt^2, mt^2, 1.043099340,  1.998988765},
  {"p2=1e6,    mt2, mt2   eq above",   10.^6,  mt^2, mt^2,-1.240391193,  2.948707299},
  {"p2=-5e4,   0,   mt2   spacelike", -50000.,   0.,  mt^2, 0.427308974, 0.},
  {"p2=-1e6,   mt2, mt2   spacelike",  -10.^6, mt^2, mt^2,-1.778323671,  0.},
  {"p2=1e6,    0,   0     massless",   10.^6,   0.,   0.,  -1.514716085, 3.141592654},
  {"p2=29929,  0,   mt2   just above", 29929.,  0.,  mt^2,  2.029719701, 0.018133253}
 }];

colValidate[tol_: 10^-7] := Module[{mu2 = 172.5^2, got, ref, nbad = 0},
  Print["=================================================================="];
  Print[" B0 via Collier,  mu^2 = mt^2 = ", mu2];
  Print["=================================================================="];
  Do[
    ref = t[[5]] + I t[[6]];
    got = B0c[t[[2]], t[[3]], t[[4]], mu2];
    If[got === $Failed,
       nbad++; Print["** ", t[[1]], "   driver failed"],
       If[Abs[got - ref] > tol, nbad++];
       Print[StringPadRight[If[Abs[got - ref] > tol, "** ", "OK "] <> t[[1]], 36],
             StringPadLeft[ToString[NumberForm[N[ref], {12, 9}]], 24],
             StringPadLeft[ToString[NumberForm[N[got], {12, 9}]], 24]]],
    {t, colRefs}];
  Print[If[nbad == 0, "  ALL 15 AGREE", "  " <> ToString[nbad] <> " MISMATCHES"]];
  nbad == 0];


(* ==================================================================== *)
(* 6.  SELF-TEST FOR C0 AND D0                                          *)
(*                                                                      *)
(*  No analytic reference for these, so we test what we can:            *)
(*    - C0 with a massive internal line is UV and IR finite, so it is   *)
(*      genuinely mu-INDEPENDENT.                                       *)
(*    - D0 with two massless gluon lines is IR divergent, so its        *)
(*      finite part DOES depend on mu.  That is correct, not a bug.     *)
(*    - above threshold both must have Im != 0; spacelike, Im = 0.      *)
(* ==================================================================== *)

colSelfTest[] := Module[{mu2 = 172.5^2, mt2 = 29929., s, t, c1, c2, d1, d2},
  s = 1000000.;  t = -391840.886785;

  Print["\n--- C0 (g t tbar vertex): must be mu-INDEPENDENT ---"];
  c1 = C0c[mt2, s, mt2, mt2, 0., mt2, mu2];
  c2 = C0c[mt2, s, mt2, mt2, 0., mt2, 2 mu2];
  Print["   mu^2   : ", N[c1, 10]];
  Print["   2 mu^2 : ", N[c2, 10]];
  Print["   diff   : ", Chop[c2 - c1], If[Chop[c2 - c1] === 0, "   OK", "   ** unexpected"]];

  Print["\n--- D0 (box): IR divergent, so mu-dependence is EXPECTED ---"];
  d1 = D0c[0., 0., mt2, mt2, s, t, 0., 0., mt2, mt2, mu2];
  d2 = D0c[0., 0., mt2, mt2, s, t, 0., 0., mt2, mt2, 2 mu2];
  Print["   mu^2   : ", N[d1, 10]];
  Print["   2 mu^2 : ", N[d2, 10]];
  Print["   d/dLog[mu^2] = ", N[(d2 - d1)/Log[2], 8]];
  Print["   (equals minus the 1/eps coefficient - that is the check)"];

  Print["\n--- imaginary parts ---"];
  Print["   C0 above threshold, Im != 0 : ", Abs[Im[c1]] > 10^-16];
  Print["   D0 above threshold, Im != 0 : ", Abs[Im[d1]] > 10^-16];
  Print["   D0 spacelike,       Im == 0 : ",
        Abs[Im[D0c[0., 0., mt2, mt2, -200000., t, 0., 0., mt2, mt2, mu2]]] < 10^-16];
];


(* ==================================================================== *)
(* 7.  RUN                                                              *)
(* ==================================================================== *)


colValidate[];
colSelfTest[];



(* ==================================================================== *)
(* 8.  HOW TO USE IN THE AMPLITUDE CODE                                 *)
(*                                                                      *)
(*  Build the arguments from the ACTUAL phase-space point, never by     *)
(*  hand - the six D0 invariants are not independent, and inconsistent  *)
(*  values put Collier into a degenerate configuration (COLI and DD     *)
(*  then disagree, which is how we found this).                         *)
(*                                                                      *)
(*    s  = (k1+k2)^2                                                    *)
(*    t  = (k2-pt)^2                                                    *)
(*    mt2 = pt.pt                                                       *)
(*                                                                      *)
(*  then e.g. the t-channel box:                                        *)
(*    D0c[0, 0, mt2, mt2, s, t, 0, 0, mt2, mt2, mu2]                    *)
(*                                                                      *)
(*  and always call the MEMOISED versions (A0c, B0c, C0c, D0c), not     *)
(*  the raw ones.                                                       *)
(* ==================================================================== *)

colLoopArgs[k1_, k2_, pt_] := Module[{g, dot, s, t, mt2},
  g = DiagonalMatrix[{1, -1, -1, -1}];
  dot[a_, b_] := a . g . b;
  s   = dot[k1 + k2, k1 + k2];
  t   = dot[k2 - pt, k2 - pt];
  mt2 = dot[pt, pt];
  <|"s" -> s, "t" -> t, "u" -> 2 mt2 - s - t, "mt2" -> mt2|>];


(* ==================================================================== *)
(* 9.  TIMING  --  check whether RunProcess is fast enough yet          *)
(* ==================================================================== *)

colTiming[] := Module[{mu2 = 172.5^2, t1, t2},
  colClearCache[];
  t1 = First@AbsoluteTiming[
     Do[B0col[62500. + n, 15687.5625, 29756.25, mu2], {n, 20}]];
  t2 = First@AbsoluteTiming[
     Do[B0c[62500., 15687.5625, 29756.25, mu2], {n, 20}]];
  Print["\n20 uncached calls : ", t1, " s   (", t1/20, " s each)"];
  Print["20 cached calls   : ", t2, " s"];
  Print["speedup from cache: ", Round[t1/Max[t2, 10.^-6]], "x"];
  Print["\nIf the per-call time is above ~10 ms, LibraryLink will be"];
  Print["needed for a phase-space scan.  Validate this interface first."];
];


(* ==================================================================== *)
(* 10.  POLE STRUCTURE                                                  *)
(*                                                                      *)
(*  The driver returns  I = c2/eps^2 + c1/eps + fin.  These are the IR  *)
(*  poles that must cancel against real emission, so they are needed    *)
(*  for the virtual correction, not just nice to have.                  *)
(*                                                                      *)
(*  Useful sanity rules:                                                *)
(*    B0  -> c1 = 1 exactly, c2 = 0        (B0 = Delta_eps + fin)       *)
(*    A0(m^2) -> c1 = m^2, c2 = 0                                       *)
(*    C0 with all internal lines massive -> finite, c1 = c2 = 0         *)
(*    D0 with two massless internal lines -> both poles present         *)
(*    D0 with one massless internal line  -> finite                     *)
(*  and every c2 should come out purely REAL: the leading soft-collinear*)
(*  pole carries no phase.                                              *)
(* ==================================================================== *)

colPoles[r_] := If[r === $Failed, Print["** failed"],
  Print["   fin   = ", ScientificForm[N[r["fin"]], 8]];
  Print["   1/eps = ", ScientificForm[N[r["c1"]], 8]];
  Print["   1/e^2 = ", ScientificForm[N[r["c2"]], 8]];
  If[Abs[Im[r["c2"]]] > 10^-16,
     Print["   ** c2 has an imaginary part - unexpected"]];
  r];

colPoleTable[mu2_: 172.5^2] := Module[
  {s = 1000000., t = -391997.844, u = -548131.156, m = 29929., row},
  row[lbl_, res_] := If[res === $Failed, {lbl, "failed", "", ""},
    {lbl, ScientificForm[N[res["fin"]], 5],
          ScientificForm[N[res["c1"]], 5],
          ScientificForm[N[res["c2"]], 5]}];
  Print[Grid[Join[
    {{"integral", "finite", "1/eps", "1/eps^2"}},
    {row["A0(mt^2)",              A0cF[m, mu2]]},
    {row["B0(s; 0,0)",            B0cF[s, 0., 0., mu2]]},
    {row["B0(s; mt,mt)",          B0cF[s, m, m, mu2]]},
    {row["B0(t; 0,mt)",           B0cF[t, 0., m, mu2]]},
    {row["B0(u; 0,mt)",           B0cF[u, 0., m, mu2]]},
    {row["B0(mt^2; 0,mt)",        B0cF[m, 0., m, mu2]]},
    {row["C0(0,0,s; 0,0,0)",      C0cF[0., 0., s, 0., 0., 0., mu2]]},
    {row["C0(0,0,s; mt,mt,mt)",   C0cF[0., 0., s, m, m, m, mu2]]},
    {row["C0(mt2,mt2,s; mt,0,mt)",C0cF[m, m, s, m, 0., m, mu2]]},
    {row["C0(0,mt2,t; 0,0,mt)",   C0cF[0., m, t, 0., 0., m, mu2]]},
    {row["D0 t-ch, 2 massless",   D0cF[0., 0., m, m, s, t, 0., 0., m, m, mu2]]},
    {row["D0 u-ch, 2 massless",   D0cF[0., 0., m, m, s, u, 0., 0., m, m, mu2]]},
    {row["D0 t-ch, 1 massless",   D0cF[0., 0., m, m, s, t, 0., m, m, m, mu2]]}],
   Frame -> All, Alignment -> Left]];
];

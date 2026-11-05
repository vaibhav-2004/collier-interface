(* ::Package:: *)

(* ==================================================================== *)
(*                                                                      *)
(*   gg -> t tbar  SPIN-DENSITY MATRIX                                  *)
(*   LO (complete) and NLO (infrastructure + assembly)                  *)
(*                                                                      *)
(*   ------------------------------------------------------------------ *)
(*   CONTENTS                                                            *)
(*                                                                       *)
(*     PART 1   setup: assumptions, kinematics, spin basis,             *)
(*              polarisations, Dirac algebra, spinors                   *)
(*     PART 2   amplitudes, colour, the Ward identity                   *)
(*     PART 3   the trace test, the 4x4 matrix, Bloch-Fano, targets     *)
(*     PART 4   numerics and evaluation at any phase-space point        *)
(*     PART 5   NLO conventions and the analytic B0                     *)
(*     PART 6   the Collier interface                                   *)
(*     PART 7   the one-loop integral list        (Task 1a)             *)
(*     PART 8   the colour matrix on the extended basis                 *)
(*     PART 9   assembling R at one loop                                *)
(*     PART 10  the physics check: T-odd coefficients switch on         *)
(*                                                                       *)
(*   ------------------------------------------------------------------ *)
(*   HOW TO RUN                                                          *)
(*                                                                       *)
(*   Run Quit[] first, then evaluate PART 1 to PART 4 as a single block. *)
(*   Do not leave long gaps between cells: a recycled kernel loses every *)
(*   definition and the failure is SILENT - you get plausible-looking    *)
(*   output full of unevaluated symbols.  The sanity gate at the end of  *)
(*   PART 2 catches this in two seconds.                                 *)
(*                                                                       *)
(*   PARTS 5-10 need the NLO inputs:                                     *)
(*     - the Collier interface (PART 6) calls an external program        *)
(*     - the integral list (PART 7) is read from a .mx file produced in  *)
(*       WSL by ttbar_stages.wls, because FormCalc CANNOT run in         *)
(*       Windows Mathematica (MathLink will not launch a process across  *)
(*       the Windows/WSL boundary).  Reading a .mx across that boundary  *)
(*       is a plain file read, so it works fine.                         *)
(*                                                                       *)
(*   ------------------------------------------------------------------ *)
(*   WHAT THE CALCULATION DOES, IN ONE PARAGRAPH                         *)
(*                                                                       *)
(*   We want the spin state of the t tbar pair, not just the rate.  So   *)
(*   instead of summing over final spins (which gives one number) we     *)
(*   keep all 16 helicity amplitudes explicitly, sum and average ONLY    *)
(*   over the initial gluon helicities and colours, and leave the two    *)
(*   top spin labels open.  Those become the row and column of a 4x4     *)
(*   matrix R.  Its trace is the ordinary |M|^2; the other fifteen       *)
(*   numbers are new information that a cross section cannot see.        *)
(*                                                                       *)
(*   ------------------------------------------------------------------ *)
(*   DESIGN CHOICES THAT MATTER                                          *)
(*                                                                       *)
(*     - $Assumptions declares th, mt, gs REAL.  Without this every      *)
(*       uBar carries Conjugate[th], Conjugate[Sqrt[...]] and            *)
(*       FullSimplify wastes its time unwinding them.  gs especially:    *)
(*       colorProduct takes a Conjugate, so without gs real the trace    *)
(*       comes out as gs^2 Conjugate[gs]^2 and never matches gs^4.       *)
(*     - Everything expensive is MEMOIZED:  f[x_] := f[x] = ...          *)
(*       Building R calls the amplitude routine 128 times; without       *)
(*       caching each call re-runs FullSimplify.  This alone took the    *)
(*       notebook from "cloud timeout" to 1.6 seconds.                   *)
(*     - Energy is scaled out (Et = 1).  A 2->2 amplitude has mass       *)
(*       dimension 4 - 4 = 0, so this is EXACT, and it removes Sqrt[s]   *)
(*       from every spinor and propagator.                               *)
(*     - The TRACE is computed before the 4x4 matrix: it needs only      *)
(*       four diagonal sums and is the decisive normalisation test.      *)
(*     - R is built element by element (upper triangle only), not as     *)
(*       one giant Table plus a global FullSimplify.                     *)
(*     - The spin-projection operators are NOT sigma_1,2,3.  See the     *)
(*       derivation in PART 3.                                           *)
(* ==================================================================== *)

Quit[]   (* run this first, then evaluate PART 1 onwards in one go *)


(* ==================================================================== *)
(*  PART 1  --  SETUP                                                   *)
(* ==================================================================== *)

(* -------------------------------------------------------------------- *)
(* 1.0  ASSUMPTIONS                                                      *)
(*                                                                       *)
(* Mathematica treats every symbol as COMPLEX unless told otherwise.     *)
(* uBar[h] = Conjugate[uSpin[h]].gam0 therefore generates Conjugate[th], *)
(* Conjugate[Sqrt[1+mt]] and so on, and those survive into every         *)
(* amplitude.  Declaring the symbols real here is the single biggest     *)
(* speed improvement in the whole notebook.  Note gs in particular:      *)
(* colorProduct takes Conjugate[b], so without gs real the answer comes  *)
(* out as gs^2 Conjugate[gs]^2 and will NOT match gs^4 |M|^2.            *)
(* -------------------------------------------------------------------- *)

Clear[assumptions, simp, simpN];

assumptions = {
   Element[{th, gs}, Reals],
   Element[bet, Reals],         (* bet is the FREE symbol, mt is derived *)
   0 < bet < 1,                 (* top velocity: 0 at threshold, 1 at    *)
   0 < th < Pi,                 (* high energy                           *)
   gs > 0
};

(* WHY bet is the variable of choice:                                    *)
(*   bet = Sqrt[1 - 4 mt^2/s] is the top VELOCITY in the ttbar rest      *)
(*   frame.  It runs 0 -> 1 from threshold to high energy, it is what    *)
(*   the Appendix-B reference formulas are written in, and it is what    *)
(*   bernreuther_rmatrix.py and the rest of the literature use.  Keeping *)
(*   results in bet means no translation step when comparing.            *)
(*                                                                       *)
(* TRADE-OFF (known):  with bet free, mt = Sqrt[1-bet^2], so the spinor  *)
(*   norms become Sqrt[1 +- Sqrt[1-bet^2]] - a NESTED radical that       *)
(*   FullSimplify cannot see through.  The spinor residuals in 1.5 then  *)
(*   print as long expressions that are secretly zero.  Check them       *)
(*   numerically instead (see 1.5c).                                     *)

$Assumptions = assumptions;

simp[x_]  := FullSimplify[x, assumptions];
simpN[x_] := Simplify[x, assumptions];   (* cheaper; use inside hot loops *)


(* -------------------------------------------------------------------- *)
(* 1.1  KINEMATICS  (partonic CM frame = t tbar rest frame at LO)        *)
(*                                                                       *)
(*  In 4 dimensions an amplitude with n external legs has mass dimension *)
(*  4 - n.  For a 2 -> 2 process n = 4, so the amplitude is EXACTLY      *)
(*  dimensionless - it cannot depend on the overall energy scale, only   *)
(*  on ratios.  (Check it directly: the two spinors give Sqrt[E].Sqrt[E],*)
(*  the propagator numerator gives E, the denominator E^2, net E^0.)     *)
(*                                                                       *)
(*  So we are free to fix the CM energy per particle to 1:               *)
(*      Et = Sqrt[s]/2 = 1   =>   s = (2 Et)^2 = 4                       *)
(*  This is EXACT, not an approximation - verified entry by entry        *)
(*  against the version carrying s and mt explicitly.  The payoff is     *)
(*  that Sqrt[s] disappears from every spinor normalisation and every    *)
(*  propagator, so everything depends on bet and th alone.               *)
(*                                                                       *)
(*  In these units mt is the RATIO 2 m_t^phys / Sqrt[s^phys], in (0,1):  *)
(*      mt = Sqrt[1 - bet^2],   mt -> 1 at threshold, 0 at high energy.  *)
(*  PART 4 converts back to physical s and m_t.                          *)
(* -------------------------------------------------------------------- *)

Clear[Et, sMand, mt, bet, eta, mdot, k1v, k2v, ptv, ptbv];

Et    = 1;                       (* CM energy per particle, scaled to 1   *)
sMand = 4;                       (* therefore s = (2 Et)^2 = 4            *)
mt    = Sqrt[1 - bet^2];         (* rescaled mass, DERIVED from bet       *)

eta = DiagonalMatrix[{1, -1, -1, -1}];          (* metric (+,-,-,-)     *)
mdot[a_List, b_List] := a . eta . b;            (* Minkowski product    *)

(* incoming gluons along +-z *)
k1v = Et {1, 0, 0,  1};
k2v = Et {1, 0, 0, -1};

(* outgoing top / antitop, scattering angle th in the x-z plane *)
ptv  = Et {1,  bet Sin[th], 0,  bet Cos[th]};
ptbv = Et {1, -bet Sin[th], 0, -bet Cos[th]};

kinCheck = simp[{
   mdot[k1v, k1v], mdot[k2v, k2v],
   mdot[ptv, ptv] - mt^2, mdot[ptbv, ptbv] - mt^2,
   k1v + k2v - ptv - ptbv
}];
Print["1.1 kinematics {k1^2, k2^2, pt^2-mt^2, ptb^2-mt^2, cons} = ", kinCheck];


(* -------------------------------------------------------------------- *)
(* 1.2  SPIN BASIS  {r, n, k}                                            *)
(*                                                                       *)
(* We need three axes to write the spin correlations against.  Any       *)
(* orthonormal triad works mathematically; this one is chosen so the     *)
(* zeros are obvious.  k = top direction, r = the other in-plane         *)
(* direction, n = normal to the plane containing the beam and the top.   *)
(* At LO every coefficient with an ODD number of n indices vanishes      *)
(* (they are T-odd and need an absorptive part), so this split makes     *)
(* four free checks visible.  In a random basis they would be smeared    *)
(* across many coefficients and the checks would be lost.                *)
(*                                                                       *)
(* CONVENTION WARNING.  n = r x k gives n = (0,1,0); n = k x r gives     *)
(* the opposite.  At LO in the SM this choice is INVISIBLE, because      *)
(* C_rn = C_nr = C_nk = C_kn = 0 and C_nn is even in n.  It becomes      *)
(* visible only for CP-odd observables, i.e. exactly the BSM physics     *)
(* this project is aimed at.  Settle it against bernreuther_rmatrix.py   *)
(* before the BSM stage, and record the choice in the repository.        *)
(* -------------------------------------------------------------------- *)

Clear[pHat, kHat, rHat, nHat];

pHat = {0, 0, 1};                             (* beam direction          *)
kHat = {Sin[th], 0, Cos[th]};                 (* top direction           *)
rHat = simp[(pHat - Cos[th] kHat)/Sin[th]];   (* -> {-Cos[th],0,Sin[th]} *)
nHat = simp[Cross[rHat, kHat]];               (* -> {0,1,0}              *)

basisCheck = simp[{
   kHat . kHat, rHat . rHat, nHat . nHat,     (* want 1,1,1   *)
   kHat . rHat, kHat . nHat, rHat . nHat,     (* want 0,0,0   *)
   Cross[rHat, kHat] - nHat                   (* want {0,0,0} *)
}];
Print["1.2 basis {norms, orthogonality, handedness} = ", basisCheck];


(* -------------------------------------------------------------------- *)
(* 1.3  INCOMING GLUON POLARIZATION VECTORS                              *)
(*                                                                       *)
(* Gluon 1 moves along +z, gluon 2 along -z.  eps2 is obtained from      *)
(* eps1 by a rotation of pi about the y axis (Jacob-Wick), which flips   *)
(* the x component.  Incoming particles use eps, not eps*.               *)
(* -------------------------------------------------------------------- *)

Clear[eps1, eps2];

eps1[+1] = -(1/Sqrt[2]) {0,  1,  I, 0};
eps1[-1] =  (1/Sqrt[2]) {0,  1, -I, 0};
eps2[+1] = -(1/Sqrt[2]) {0, -1,  I, 0};
eps2[-1] =  (1/Sqrt[2]) {0, -1, -I, 0};

(* transversality, normalisation, and the fact that summing over the     *)
(* two physical helicities reproduces the ghost-free projector used in   *)
(* the old spin-summed calculation.                                      *)
pol1FromHelicities = Sum[Outer[Times, eps1[h], Conjugate[eps1[h]]], {h, {1, -1}}];
pol2FromHelicities = Sum[Outer[Times, eps2[h], Conjugate[eps2[h]]], {h, {1, -1}}];
pol1Physical = -eta + (Outer[Times, k1v, k2v] + Outer[Times, k2v, k1v])/mdot[k1v, k2v];
pol2Physical = -eta + (Outer[Times, k2v, k1v] + Outer[Times, k1v, k2v])/mdot[k2v, k1v];

polCheck = simp[{
   Table[mdot[k1v, eps1[h]], {h, {1, -1}}],                  (* want 0,0   *)
   Table[mdot[k2v, eps2[h]], {h, {1, -1}}],                  (* want 0,0   *)
   Table[mdot[eps1[h], Conjugate[eps1[h]]], {h, {1, -1}}],   (* want -1,-1 *)
   Flatten[pol1FromHelicities - pol1Physical],               (* want all 0 *)
   Flatten[pol2FromHelicities - pol2Physical]                (* want all 0 *)
}];
Print["1.3 polarizations {k.e, k.e, e.e*, complete1, complete2} = ", polCheck];


(* -------------------------------------------------------------------- *)
(* 1.4  DIRAC ALGEBRA  (Dirac representation, metric (+,-,-,-))          *)
(*                                                                       *)
(* The only property of the gamma matrices that matters is the Clifford  *)
(* algebra {gam^m, gam^n} = 2 g^mn; everything else follows from it, so  *)
(* that is what we check.  slash4[p] = gam^mu p_mu takes CONTRAVARIANT   *)
(* components and the metric supplies the minus signs.                   *)
(* -------------------------------------------------------------------- *)

Clear[id2, zero2, sig1, sig2, sig3, gam0, gam1, gam2, gam3,
      gammaList, id4, slash4, sigDot, Sigma4, lower4];

id2   = IdentityMatrix[2];
zero2 = ConstantArray[0, {2, 2}];

sig1 = {{0, 1}, {1, 0}};
sig2 = {{0, -I}, {I, 0}};
sig3 = {{1, 0}, {0, -1}};

gam0 = ArrayFlatten[{{id2, zero2}, {zero2, -id2}}];
gam1 = ArrayFlatten[{{zero2, sig1}, {-sig1, zero2}}];
gam2 = ArrayFlatten[{{zero2, sig2}, {-sig2, zero2}}];
gam3 = ArrayFlatten[{{zero2, sig3}, {-sig3, zero2}}];

gammaList = {gam0, gam1, gam2, gam3};
id4       = IdentityMatrix[4];

(* one-line definitions on purpose: a line break after a syntactically   *)
(* complete expression silently terminates the input in Mathematica      *)
slash4[p_List] := gam0 p[[1]] - gam1 p[[2]] - gam2 p[[3]] - gam3 p[[4]];
sigDot[v_List] := v[[1]] sig1 + v[[2]] sig2 + v[[3]] sig3;
Sigma4[v_List] := ArrayFlatten[{{sigDot[v], zero2}, {zero2, sigDot[v]}}];
lower4[v_List] := eta . v;

cliffordCheck = simp[Union[Flatten[
   Table[gammaList[[m]] . gammaList[[n]] + gammaList[[n]] . gammaList[[m]]
         - 2 eta[[m, n]] id4, {m, 4}, {n, 4}]]]];
Print["1.4 Clifford {gam^m,gam^n}=2g^mn residual (want {0}) = ", cliffordCheck];


(* -------------------------------------------------------------------- *)
(* 1.5  TOP AND ANTITOP HELICITY SPINORS                                 *)
(*                                                                       *)
(* chiT[h]  : eigenspinor of sigma.kHat  with eigenvalue h               *)
(* chiTb[h] : eigenspinor of sigma.(-kHat) with eigenvalue h             *)
(*            (the antitop moves along -kHat)                            *)
(*                                                                       *)
(* u_h = ( a chi , h b chi ) , v_h = ( h b chibar , a chibar )           *)
(* with a = Sqrt[Et+mt], b = Sqrt[Et-mt].  Note b = |p|/a, which is why  *)
(* these solve the Dirac equation.                                       *)
(*                                                                       *)
(* IMPORTANT: h on vSpin is the label carried by chiTb, NOT the physical *)
(* antitop spin.  See PART 3 for the correct projection operator.        *)
(* -------------------------------------------------------------------- *)

Clear[chiT, chiTb, aSpin, bSpin, uSpin, vSpin, uBar, vBar];

chiT[+1]  = { Cos[th/2],  Sin[th/2]};
chiT[-1]  = {-Sin[th/2],  Cos[th/2]};
chiTb[+1] = { Sin[th/2], -Cos[th/2]};
chiTb[-1] = { Cos[th/2],  Sin[th/2]};

aSpin = Sqrt[Et + mt];
bSpin = Sqrt[Et - mt];

uSpin[h_] := uSpin[h] = Join[aSpin chiT[h],    h bSpin chiT[h]];
vSpin[h_] := vSpin[h] = Join[h bSpin chiTb[h], aSpin chiTb[h]];
uBar[h_]  := uBar[h]  = simp[Conjugate[uSpin[h]]] . gam0;
vBar[h_]  := vBar[h]  = simp[Conjugate[vSpin[h]]] . gam0;

spinorCheck = simp[Union[Flatten[Table[{
   (slash4[ptv]  - mt id4) . uSpin[h],        (* Dirac eqn, want 0     *)
   (slash4[ptbv] + mt id4) . vSpin[h],        (* Dirac eqn, want 0     *)
   Sigma4[kHat]  . uSpin[h] - h uSpin[h],     (* helicity eigenvalue   *)
   Sigma4[-kHat] . vSpin[h] - h vSpin[h]      (* helicity eigenvalue   *)
}, {h, {1, -1}}]]]];
normCheck = simp[Table[
   {uBar[h] . uSpin[h], vBar[h] . vSpin[h], Conjugate[uSpin[h]] . uSpin[h]},
   {h, {1, -1}}]];
Print["1.5a spinor residuals (want {0}) = ", spinorCheck];
Print["1.5b {ubar u, vbar v, u^dag u} (want {2mt, -2mt, 2Et}) = ", normCheck];

(* 1.5c  NUMERICAL spinor check                                          *)
(*                                                                       *)
(* With bet as the free symbol the residuals in 1.5a are proportional to *)
(*   aSpin (1 - mt) - bet bSpin ,                                        *)
(* which vanishes because 1 - mt = bSpin^2 and bet = aSpin bSpin.  But   *)
(* FullSimplify cannot prove it through the nested Sqrt[1 +- Sqrt[1-b^2]]*)
(* so 1.5a prints long expressions that are secretly zero.  Check them   *)
(* numerically instead - this is conclusive.                             *)

spinorCheckNum = Chop[Max[Abs[Flatten[
   N[spinorCheck /. {bet -> 0.8693, th -> 1.1}]]]]];
Print["1.5c spinor residuals, NUMERICAL (want 0) = ", spinorCheckNum];


(* ==================================================================== *)
(*  PART 2  --  AMPLITUDES, COLOUR, THE WARD IDENTITY                   *)
(* ==================================================================== *)

(* -------------------------------------------------------------------- *)
(* 2.1  FIXED-HELICITY AMPLITUDES                                        *)
(*                                                                       *)
(* Colour is stripped off and handled separately in a 3-dim basis        *)
(*     c1 = T^a T^b  (t-channel), c2 = T^b T^a  (u-channel),             *)
(*     c3 = f^{abc} T^c  (s-channel)                                     *)
(* so At, Au, As below are the COLOUR-STRIPPED coefficients.             *)
(*                                                                       *)
(* NOTE for the NLO part: this tree basis is NOT the one-loop basis.     *)
(* See PART 8 - f^abc T^c = -i(c1 - c2), so as matrices the tree basis   *)
(* spans only 2 dimensions, and at one loop delta^ab*1 appears as a      *)
(* genuinely new structure.  PART 8 converts between the two.            *)
(*                                                                       *)
(* Internal momenta follow the fermion arrow:  qT = pt - k1 satisfies    *)
(* qT + k1 = pt at the top-side vertex.  Same for qU = pt - k2.          *)
(*                                                                       *)
(* Every definition is memoized (f[x_] := f[x] = ...).  Without this,    *)
(* building R calls Avec 128 times and re-runs FullSimplify each time.   *)
(* -------------------------------------------------------------------- *)

Clear[qTv, qUv, denT, denU, propT, propU,
      epsSlash1, epsSlash2, At, Au, As, Avec,
      p3v, V3Tensor, v3Contract, qCurrent];

qTv = ptv - k1v;
qUv = ptv - k2v;

denT = simp[mdot[qTv, qTv] - mt^2];    (* = t - mt^2 = -(s/2)(1 - b Cos) *)
denU = simp[mdot[qUv, qUv] - mt^2];    (* = u - mt^2 = -(s/2)(1 + b Cos) *)
Print["2.1a {denT, denU} = ", {denT, denU}];

propT = simp[(slash4[qTv] + mt id4)/denT];
propU = simp[(slash4[qUv] + mt id4)/denU];

epsSlash1[h_] := epsSlash1[h] = slash4[eps1[h]];
epsSlash2[h_] := epsSlash2[h] = slash4[eps2[h]];

(* --- t channel: gluon 1 attaches at the ubar end --- *)
At[h1_, h2_, ht_, hb_] := At[h1, h2, ht, hb] = simp[
   uBar[ht] . (I gs epsSlash1[h1]) . (I propT) . (I gs epsSlash2[h2]) . vSpin[hb]];

(* --- u channel: the two gluon insertions are exchanged --- *)
Au[h1_, h2_, ht_, hb_] := Au[h1, h2, ht, hb] = simp[
   uBar[ht] . (I gs epsSlash2[h2]) . (I propU) . (I gs epsSlash1[h1]) . vSpin[hb]];

(* --- s channel: three-gluon vertex, all momenta incoming --- *)
p3v = -(k1v + k2v);

V3Tensor[m_, n_, sg_] := gs (
     eta[[m, n]]  (k1v - k2v)[[sg]]
   + eta[[n, sg]] (k2v - p3v)[[m]]
   + eta[[sg, m]] (p3v - k1v)[[n]]);

(* contract the two external gluon indices; lower4 supplies eps_mu so    *)
(* the result X^sigma is contravariant                                   *)
v3Contract[h1_, h2_] := v3Contract[h1, h2] = simp[Table[
   Sum[lower4[eps1[h1]][[m]] lower4[eps2[h2]][[n]] V3Tensor[m, n, sg],
       {m, 4}, {n, 4}], {sg, 4}]];

(* the fixed-helicity top current J^rho = ubar gam^rho v.                *)
(* simp[] here matters: without it Conjugate[] wrappers propagate into   *)
(* all 16 s-channel amplitudes.                                          *)
qCurrent[ht_, hb_] := qCurrent[ht, hb] = simp[
   Table[uBar[ht] . gammaList[[r]] . vSpin[hb], {r, 4}]];

(* quark vertex (I gs) times gluon propagator (-I g_rs / s)              *)
As[h1_, h2_, ht_, hb_] := As[h1, h2, ht, hb] = simp[
   (I gs) (-I)/sMand mdot[qCurrent[ht, hb], v3Contract[h1, h2]]];

Avec[h1_, h2_, ht_, hb_] := Avec[h1, h2, ht, hb] =
   {At[h1, h2, ht, hb], Au[h1, h2, ht, hb], As[h1, h2, ht, hb]};

helicities  = {+1, -1};
finalStates = {{+1, +1}, {+1, -1}, {-1, +1}, {-1, -1}};

(* check the three-gluon contraction against the closed form that        *)
(* follows from eps1.k2 = eps2.k1 = 0 in this frame                      *)
v3Check = simp[Union[Flatten[Table[
   v3Contract[h1, h2] - gs mdot[eps1[h1], eps2[h2]] (k1v - k2v),
   {h1, helicities}, {h2, helicities}]]]];
Print["2.1b three-gluon contraction residual (want {0}) = ", v3Check];

(* build and cache all 16 amplitude vectors ONCE, up front *)
Print["2.1c caching 16 amplitudes, timing: ",
   AbsoluteTiming[Do[Avec[h1, h2, ht, hb],
      {h1, helicities}, {h2, helicities},
      {ht, helicities}, {hb, helicities}];][[1]], " s"];


(* -------------------------------------------------------------------- *)
(* 2.2  COLOUR AT TREE LEVEL                                             *)
(*                                                                       *)
(* Colour is kept completely separate from the Dirac algebra.  Putting   *)
(* SU(3) generators inside a Dirac chain means two different             *)
(* non-commutative structures in one product, which is what broke the    *)
(* original FeynCalc attempt ("failed to isolate the colour algebraic    *)
(* expression").  Here c1, c2, c3 are ordinary COMMUTING symbols, so     *)
(* each amplitude is just a 3-vector and the Dirac algebra never sees    *)
(* them.  The three diagrams are never added together: the bilinear      *)
(* form A_i^T . C . A_j^* automatically produces all nine diagram-pair   *)
(* interference terms with the right colour weights.                     *)
(*                                                                       *)
(* C_mn = Sum_colours c_m c_n^dagger with the quark indices traced:      *)
(*   C11 = C22 = 16/3 , C12 = C21 = -2/3 , C33 = 12 ,                    *)
(*   C13 = -C31 = 6i , C23 = -C32 = -6i                                  *)
(* using Tr(T^a T^b) = d^ab/2, T^a T^b T^a = -T^b/6, f^abc f^abc = 24,   *)
(* Tr(T^a T^b T^c) = (d^abc + i f^abc)/4.                                *)
(* -------------------------------------------------------------------- *)

Clear[Ccolor, colorProduct];

Ccolor = {
   {16/3, -2/3,  6 I},
   {-2/3, 16/3, -6 I},
   {-6 I,  6 I,   12}
};

(* NOTE the _List patterns.  If Avec is ever undefined, its head is      *)
(* Avec rather than List, this pattern silently fails to match, and the  *)
(* whole calculation returns unevaluated symbols with NO error message.  *)
(* The guard below turns that silent failure into a loud one.            *)
colorProduct[a_List, b_List] := a . Ccolor . Conjugate[b];
colorProduct::badarg = "colorProduct received non-list arguments: `1`, `2`. \
Definitions are missing - the kernel has probably been restarted.";
colorProduct[a_, b_] := (Message[colorProduct::badarg, a, b]; Abort[]);

Print["2.2a Ccolor Hermitian: ", ConjugateTranspose[Ccolor] === Ccolor];
Print["2.2b colorProduct[++;++] (expect 0, c1-c2 = i c3 cancellation) = ",
      simp[colorProduct[Avec[1, 1, 1, 1], Avec[1, 1, 1, 1]]]];


(* -------------------------------------------------------------------- *)
(* 2.3  WARD IDENTITY  --  the strongest single check in the notebook    *)
(*                                                                       *)
(* Gauge invariance says: replace an external gluon's polarisation       *)
(* vector by its own momentum, and the amplitude must vanish.            *)
(*                                                                       *)
(*      eps1 -> k1     =>     k1_mu M^{mu nu} = 0                        *)
(*                                                                       *)
(* WHY THIS IS WORTH DOING: it only holds for the SUM of the three       *)
(* diagrams, so one number simultaneously tests                          *)
(*    - the relative signs of the t, u and s channels                    *)
(*    - the three-gluon vertex convention (signs and momentum routing)   *)
(*    - the internal propagator momenta q_T = pt - k1 etc.               *)
(*    - the external polarisation treatment                              *)
(* Nothing else in the notebook probes those four things together.       *)
(*                                                                       *)
(* NOTE ON THE FORM: because c1 - c2 = i c3 (the SU(3) commutator), the  *)
(* three colour-stripped pieces do NOT each vanish separately.  The      *)
(* gauge-invariant combination is  i (At - Au) + As = 0.                 *)
(* -------------------------------------------------------------------- *)

Clear[epsWard, AtWard, AuWard, v3Ward, AsWard, wardCheck];

epsWard = k1v;                    (* replace eps1 by k1 *)

AtWard[h2_, ht_, hb_] := AtWard[h2, ht, hb] = simp[
   uBar[ht] . (I gs slash4[epsWard]) . (I propT)
            . (I gs epsSlash2[h2]) . vSpin[hb]];

AuWard[h2_, ht_, hb_] := AuWard[h2, ht, hb] = simp[
   uBar[ht] . (I gs epsSlash2[h2]) . (I propU)
            . (I gs slash4[epsWard]) . vSpin[hb]];

v3Ward[h2_] := v3Ward[h2] = simp[Table[
   Sum[lower4[epsWard][[m]] lower4[eps2[h2]][[n]] V3Tensor[m, n, sg],
       {m, 4}, {n, 4}], {sg, 4}]];

AsWard[h2_, ht_, hb_] := AsWard[h2, ht, hb] = simp[
   (I gs) (-I)/sMand mdot[qCurrent[ht, hb], v3Ward[h2]]];

wardCheck = simp[Union[Flatten[Table[
   I (AtWard[h2, ht, hb] - AuWard[h2, ht, hb]) + AsWard[h2, ht, hb],
   {h2, helicities}, {ht, helicities}, {hb, helicities}]]]];

Print["2.3 Ward identity  i(At-Au)+As  (want {0}) = ", wardCheck];


(* -------------------------------------------------------------------- *)
(* --- SANITY GATE: run this before anything expensive ---               *)
(* If any line below is wrong, STOP.  The kernel has lost its state and  *)
(* everything downstream will be meaningless unevaluated symbols.        *)
(* -------------------------------------------------------------------- *)

Print["GATE simp defined       (want True)  : ", Head[simp[1]] =!= simp];
Print["GATE Avec head          (want List)  : ", Head[Avec[1, 1, 1, 1]]];
Print["GATE helicities         (want {1,-1}): ", helicities];
Print["GATE colorProduct works (want 0)     : ",
      simp[colorProduct[Avec[1, 1, 1, 1], Avec[1, 1, 1, 1]]]];


(* ==================================================================== *)
(*  PART 3  --  THE TRACE, THE 4x4 MATRIX, BLOCH-FANO, THE TARGETS      *)
(* ==================================================================== *)

(* -------------------------------------------------------------------- *)
(* 3.1  THE TRACE  --  the decisive test, and the cheap one              *)
(*                                                                       *)
(*   Tr R = (1/256) Sum_{h1,h2,ht,htb} A^T C A^*                         *)
(*                                                                       *)
(* By construction this is exactly the old spin-summed, initial-state    *)
(* averaged |M|^2.  Prove this BEFORE building any off-diagonal entry:   *)
(* it needs only the 4 diagonal sums, not the full 4x4.                  *)
(*                                                                       *)
(* 1/4  = two physical helicities for each of two incoming gluons        *)
(* 1/64 = 8 x 8 incoming adjoint colours                                 *)
(* -------------------------------------------------------------------- *)

Clear[traceR, t1, u1, M2ref, Fg, Atilde];

traceR = simp[(1/(4*64)) Sum[
   colorProduct[Avec[h1, h2, ht, hb], Avec[h1, h2, ht, hb]],
   {h1, helicities}, {h2, helicities},
   {ht, helicities}, {hb, helicities}]];

t1 = -(sMand/2)(1 - bet Cos[th]);      (* t - mt^2 *)
u1 = -(sMand/2)(1 + bet Cos[th]);      (* u - mt^2 *)

M2ref = gs^4 (1/(6 t1 u1) - 3/(8 sMand^2)) *
        (t1^2 + u1^2 + 4 mt^2 sMand (1 - mt^2 sMand/(t1 u1)));

Fg     = (7 + 9 bet^2 Cos[th]^2)/(192 (1 - bet^2 Cos[th]^2)^2);
Atilde = Fg (1 + 2 bet^2 Sin[th]^2 - bet^4 (1 + Sin[th]^4));

Print["3.1a Tr R                  = ", traceR];
Print["3.1b Tr R - |M|^2_standard (want 0) = ", simp[traceR - M2ref]];
Print["3.1c Tr R - 4 gs^4 Atilde  (want 0) = ", simp[traceR - 4 gs^4 Atilde]];


(* -------------------------------------------------------------------- *)
(* 3.2  THE 4x4 DENSITY MATRIX, ELEMENT BY ELEMENT                       *)
(*                                                                       *)
(*   R_ij = (1/256) Sum_{h1,h2} A_i^T C A_j^*                            *)
(*                                                                       *)
(* Basis order  |++>, |+->, |-+>, |-->  with the first label the top.    *)
(* Initial helicities and colours are summed and averaged; the final     *)
(* top spins are kept OPEN.  That is what turns a scalar rate into a     *)
(* 4x4 spin matrix.                                                      *)
(*                                                                       *)
(* R is Hermitian, so only the 10 upper-triangle entries are             *)
(* independent; at LO it is also real, and Bose symmetry of the two      *)
(* gluons plus parity cuts it further to just FIVE independent numbers.  *)
(* Computing one element at a time keeps each FullSimplify small - far   *)
(* better than one big Table followed by a global simplification.        *)
(* -------------------------------------------------------------------- *)

Clear[Rel, Rhel, RTarget];

Rel[i_, j_] := Rel[i, j] = simp[(1/(4*64)) Sum[
   colorProduct[Avec[h1, h2, Sequence @@ finalStates[[i]]],
                Avec[h1, h2, Sequence @@ finalStates[[j]]]],
   {h1, helicities}, {h2, helicities}]];

Do[Print["3.2a R[", i, ",", j, "] = ", Rel[i, j]], {i, 4}, {j, i, 4}];

Rhel = Table[If[j >= i, Rel[i, j], Conjugate[Rel[j, i]]], {i, 4}, {j, 4}];

Print["3.2b Hermiticity residual (want {0}) = ",
      simp[Union[Flatten[Rhel - ConjugateTranspose[Rhel]]]]];
Print["3.2c trace consistency    (want 0)   = ", simp[Tr[Rhel] - traceR]];

(* closed form to compare against, derived from the Appendix-B           *)
(* coefficients and verified numerically at four (bet,th) points         *)
RTarget = gs^4 2 Fg {
  { bet^2 Sin[th]^2 (1 + Cos[th]^2),
   -bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
   -bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
    bet^2 Sin[th]^4},
  {-bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
   (1-bet^2)(1 + bet^2 + bet^2 Sin[th]^4),
   (1-bet^2)(bet^2 (1 + Sin[th]^4) - 1),
    bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th]},
  {-bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
   (1-bet^2)(bet^2 (1 + Sin[th]^4) - 1),
   (1-bet^2)(1 + bet^2 + bet^2 Sin[th]^4),
    bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th]},
  { bet^2 Sin[th]^4,
    bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
    bet^2 Sqrt[1-bet^2] Sin[th]^3 Cos[th],
    bet^2 Sin[th]^2 (1 + Cos[th]^2)}
};
Print["3.2d Rhel - RTarget (want {0}) = ",
      simp[Union[Flatten[Rhel - RTarget]]]];


(* -------------------------------------------------------------------- *)
(* 3.3  BLOCH-FANO PROJECTION                                            *)
(*                                                                       *)
(*   R = Atil (1 x 1) + B+_i (sig_i x 1) + B-_i (1 x sig_i)              *)
(*                     + C_ij (sig_i x sig_j)                            *)
(*                                                                       *)
(* Using Tr sig_i = 0 and Tr(sig_i sig_j) = 2 delta_ij:                  *)
(*   Atil = Tr R /4 ,  C_ij = Tr[R (Sig_i x Sig_j)]/4 , etc.             *)
(*                                                                       *)
(* THE OPERATORS ARE NOT sig1, sig2, sig3.  Two separate effects:        *)
(*                                                                       *)
(*  (a) The spinor basis {chiT[+1], chiT[-1]} is the helicity basis,     *)
(*      so k is the quantisation axis.  In that basis                    *)
(*         sigma.kHat -> sig3 ,  sigma.nHat -> sig2 ,                    *)
(*         sigma.rHat -> -sig1                                           *)
(*      The minus sign on r is because rHat maps to -x when kHat is      *)
(*      rotated onto z.  Assuming sig_r = sig1 flips C_rk and C_kr.      *)
(*                                                                       *)
(*  (b) For the antitop the label carried by vSpin is related to the     *)
(*      physical spin by charge conjugation, so the correct operator is  *)
(*      -(sigma.a)^* , not sigma.a and not -(sigma.a).                   *)
(*      Evaluating -(sigma.a)^* in the chiTb basis happens to give the   *)
(*      SAME three matrices as for the top, which is why both slots      *)
(*      below use SigR, SigN, SigK.                                      *)
(*                                                                       *)
(* ** THE FINDING WORTH REMEMBERING **                                   *)
(*                                                                       *)
(* Tr R = |M|^2 is NECESSARY BUT NOT SUFFICIENT.  All three antitop      *)
(* prescriptions reproduce the trace EXACTLY while giving different      *)
(* correlation coefficients:                                             *)
(*                                                                       *)
(*     +sigma.a        wrong Crr, Ckk, Crk    (Cnn ok)                   *)
(*     -sigma.a        wrong Cnn              (rest ok)                  *)
(*     -(sigma.a)*     ALL CORRECT                                       *)
(*                                                                       *)
(* The right one follows from charge conjugation acting on the           *)
(* antiparticle spinor label, using sigma2 sigma* sigma2 = -sigma.       *)
(*                                                                       *)
(* LESSON: for a tomography calculation, validation means matching the   *)
(* COEFFICIENTS, not the trace.                                          *)
(* -------------------------------------------------------------------- *)

Clear[SigR, SigN, SigK, SigList, axisName, Atil, Bplus, Bminus, Cij];

SigR = -sig1;
SigN =  sig2;
SigK =  sig3;

SigList  = {SigR, SigN, SigK};
axisName = {"r", "n", "k"};

Atil   = simp[Tr[Rhel]/4];
Bplus  = Table[simp[Tr[Rhel . KroneckerProduct[SigList[[i]], id2]]/4], {i, 3}];
Bminus = Table[simp[Tr[Rhel . KroneckerProduct[id2, SigList[[i]]]]/4], {i, 3}];
Cij    = Table[simp[Tr[Rhel . KroneckerProduct[SigList[[i]], SigList[[j]]]]/4],
               {i, 3}, {j, 3}];

Print["3.3a Atil = ", Atil];
Print["3.3b B+ (want {0,0,0}) = ", Bplus];
Print["3.3c B- (want {0,0,0}) = ", Bminus];
Do[Print["3.3d C_", axisName[[i]], axisName[[j]], " = ", Cij[[i, j]]],
   {i, 3}, {j, 3}];


(* -------------------------------------------------------------------- *)
(* 3.4  MATCH THE APPENDIX-B TARGETS                                     *)
(*                                                                       *)
(* Five independent structures at LO in the SM:                          *)
(*   Atil, C_rr, C_nn, C_kk, C_rk = C_kr                                 *)
(* plus the zeros B+_i = B-_i = C_rn = C_nr = C_nk = C_kn = 0.           *)
(*                                                                       *)
(* THE EIGHT ZEROS AND WHY THEY ARE ZERO                                 *)
(*   B+_i = B-_i = 0                parity: QCD conserves it             *)
(*   C_rn = C_nr = C_nk = C_kn = 0  T-odd: need an absorptive part,      *)
(*                                  which tree level does not have       *)
(* These come back at NLO - see PART 10.                                 *)
(* -------------------------------------------------------------------- *)

Clear[CrrT, CnnT, CkkT, CrkT];

CrrT = -Fg (1 - bet^2 (2 - bet^2)(1 + Sin[th]^4));
CnnT = -Fg (1 - 2 bet^2 + bet^4 (1 + Sin[th]^4));
CkkT = -Fg (1 - (bet^2/2) Sin[2 th]^2 - bet^4 (1 + Sin[th]^4));
CrkT =  Fg Sqrt[1 - bet^2] bet^2 Sin[2 th] Sin[th]^2;

Print["3.4a Atil - gs^4 Atilde (want 0) = ", simp[Atil - gs^4 Atilde]];
Print["3.4b Crr  - gs^4 CrrT   (want 0) = ", simp[Cij[[1, 1]] - gs^4 CrrT]];
Print["3.4c Cnn  - gs^4 CnnT   (want 0) = ", simp[Cij[[2, 2]] - gs^4 CnnT]];
Print["3.4d Ckk  - gs^4 CkkT   (want 0) = ", simp[Cij[[3, 3]] - gs^4 CkkT]];
Print["3.4e Crk  - gs^4 CrkT   (want 0) = ", simp[Cij[[1, 3]] - gs^4 CrkT]];
Print["3.4f Ckr  - gs^4 CrkT   (want 0) = ", simp[Cij[[3, 1]] - gs^4 CrkT]];
Print["3.4g off-plane zeros (want {0}) = ",
      simp[Union[{Cij[[1,2]], Cij[[2,1]], Cij[[2,3]], Cij[[3,2]]}]]];


(* --- save the state so a kernel restart costs a second, not an hour -- *)

DumpSave["stage5.mx", {simp, simpN, assumptions, mt, bet, eta, mdot,
   k1v, k2v, ptv, ptbv, kHat, rHat, nHat, eps1, eps2,
   gam0, gam1, gam2, gam3, gammaList, id4, sig1, sig2, sig3,
   slash4, sigDot, Sigma4, lower4, chiT, chiTb, uSpin, vSpin, uBar, vBar,
   At, Au, As, Avec, qCurrent, v3Contract, Ccolor, colorProduct,
   helicities, finalStates, Rhel, Atil, Bplus, Bminus, Cij,
   SigR, SigN, SigK, SigList, traceR, M2ref, Fg, Atilde}];
Print["state saved to stage5.mx  (restore with  << \"stage5.mx\" )"];


(* ==================================================================== *)
(*  PART 4  --  NUMERICS AND EVALUATION AT ANY PHASE-SPACE POINT        *)
(* ==================================================================== *)

(* -------------------------------------------------------------------- *)
(* 4.1  NUMERICAL SPOT CHECKS                                            *)
(*                                                                       *)
(* Purpose: separate a PHYSICS error from a FullSimplify timeout.  If    *)
(* the numbers here are right but the symbolic check above is not zero,  *)
(* the algebra is fine and the simplifier just gave up.                  *)
(*                                                                       *)
(* Run these FIRST whenever something is slow or suspicious.             *)
(* Reference values (gs = 1), computed independently:                    *)
(*   bet = 0.86930, th = 1.1  ->  Tr R = 0.3113740409                    *)
(*   bet = 0.63938, th = 2.3  ->  Tr R = 0.3317768301                    *)
(* -------------------------------------------------------------------- *)

Clear[numCheck];

numCheck[b0_, th0_] := Module[{rule, Rn},
   rule = {bet -> b0, th -> th0, gs -> 1};
   Rn = N[Table[(1/(4*64)) Sum[
        colorProduct[Avec[h1, h2, Sequence @@ finalStates[[i]]],
                     Avec[h1, h2, Sequence @@ finalStates[[j]]]],
        {h1, helicities}, {h2, helicities}] /. rule, {i, 4}, {j, 4}]];
   {
    "TrR"        -> Chop[Tr[Rn]],
    "TrR-|M|^2"  -> Chop[Tr[Rn] - N[M2ref /. rule]],
    "Hermitian"  -> Chop[Max[Abs[Flatten[Rn - ConjugateTranspose[Rn]]]]],
    "diag>=0"    -> (Min[Re[Diagonal[Rn]]] >= -10^-10),
    "R-RTarget"  -> Chop[Max[Abs[Flatten[Rn - N[RTarget /. rule]]]]]
   }];

Print["4.1a ", numCheck[0.86930, 1.1]];
Print["4.1b ", numCheck[0.63938, 2.3]];
Print["4.1c ", numCheck[0.40000, 0.7]];


(* -------------------------------------------------------------------- *)
(* 4.2  RESTORING PHYSICAL UNITS                                         *)
(*                                                                       *)
(* Everything above used Et = 1.  R and all Bloch-Fano coefficients are  *)
(* dimensionless, so nothing needs rescaling.  To express results in     *)
(* s and mt instead of bet, substitute                                   *)
(*     bet -> Sqrt[1 - 4 mtPhys^2/sPhys]                                 *)
(* and, if you want Mandelstam variables,                                *)
(*     t -> mtPhys^2 - sPhys/2 + sPhys bet Cos[th]/2                     *)
(*     u -> mtPhys^2 - sPhys/2 - sPhys bet Cos[th]/2                     *)
(* -------------------------------------------------------------------- *)

toPhysical = {bet -> Sqrt[1 - 4 mtPhys^2/sPhys]};
Print["4.2 Atil in (s, mt): ", Simplify[Atil /. toPhysical,
      {sPhys > 4 mtPhys^2, mtPhys > 0, 0 < th < Pi, gs > 0}]];


(* -------------------------------------------------------------------- *)
(* 4.3  BUILDING A PHASE-SPACE POINT                                     *)
(*                                                                       *)
(* The four momenta are NOT independent: they must satisfy momentum      *)
(* conservation and put both tops on shell.  Typing them by hand risks   *)
(* violating one silently.  makePoint builds them from (sqrt(s), mt,     *)
(* theta, phi) so both hold by construction.                             *)
(*                                                                       *)
(*   makePoint[1000., 173., 1.7383393226, 3.66 Degree]                   *)
(*                                                                       *)
(* gives Marcel's reference point.                                       *)
(* -------------------------------------------------------------------- *)

Clear[makePoint];

makePoint[rs_, mtPh_, th0_, phi0_: 0.] := Module[{E0, be, p3},
  E0 = rs/2;
  be = Sqrt[1 - 4 mtPh^2/rs^2];
  p3 = E0 be {Sin[th0] Cos[phi0], Sin[th0] Sin[phi0], Cos[th0]};
  N[{{E0, 0, 0, E0}, {E0, 0, 0, -E0},
     Prepend[p3, E0], Prepend[-p3, E0]}]];


(* -------------------------------------------------------------------- *)
(* 4.4  EXTRACT INVARIANTS AND SPIN BASIS FROM FOUR MOMENTA              *)
(*       Input: {k1, k2, pt, ptbar}, each {E, px, py, pz}                *)
(* -------------------------------------------------------------------- *)

Clear[psExtract, psEvaluate, psChecks, psEntanglement, psReport,
      psSelfTest, psWhichFree, numAssert];

psExtract[mom_List] := Module[
  {k1, k2, pt, ptb, g, dot, sVal, mtVal, betVal, p3, kH, pH, cth, thVal,
   rH, nH, phiVal, tot, ok},

  {k1, k2, pt, ptb} = N[mom];
  g = DiagonalMatrix[{1, -1, -1, -1}];
  dot[a_, b_] := a . g . b;

  ok = True;
  tot = k1 + k2 - pt - ptb;
  If[Max[Abs[tot]] > 10^-8,
     Print["** momentum NOT conserved: ", tot]; ok = False];
  If[Abs[dot[k1, k1]] > 10^-8 || Abs[dot[k2, k2]] > 10^-8,
     Print["** incoming gluons not massless"]; ok = False];
  If[Abs[dot[pt, pt] - dot[ptb, ptb]] > 10^-6,
     Print["** pt^2 != ptbar^2"]; ok = False];
  If[Max[Abs[(k1 + k2)[[2 ;; 4]]]] > 10^-8,
     Print["** NOT in the CM frame. Boost first - the (r,n,k) basis is ",
           "defined in the ttbar rest frame."]; ok = False];

  sVal   = dot[k1 + k2, k1 + k2];
  mtVal  = Sqrt[dot[pt, pt]];
  betVal = Sqrt[1 - 4 mtVal^2/sVal];

  p3    = pt[[2 ;; 4]];
  kH    = p3/Norm[p3];
  pH    = {0., 0., 1.};
  cth   = kH . pH;
  thVal = ArcCos[cth];
  rH    = (pH - cth kH)/Sin[thVal];
  nH    = Cross[rH, kH];
  phiVal = ArcTan[p3[[1]], p3[[2]]];

  <|"ok" -> ok, "sqrtS" -> Sqrt[sVal], "s" -> sVal, "mt" -> mtVal,
    "beta" -> betVal, "cosTheta" -> cth, "theta" -> thVal, "phi" -> phiVal,
    "mtRescaled" -> 2 mtVal/Sqrt[sVal],
    "kHat" -> kH, "rHat" -> rH, "nHat" -> nH|>
];


(* -------------------------------------------------------------------- *)
(* 4.5  DIAGNOSTIC: which of mt / bet is the free symbol?                *)
(*      Run this first if anything looks wrong.                          *)
(* -------------------------------------------------------------------- *)

psWhichFree[] := (
  Print["4.5 symbol status"];
  Print["   mt  evaluates to : ", mt];
  Print["   bet evaluates to : ", bet];
  Print["   (exactly one should be a bare symbol; the other a formula)"]);


(* -------------------------------------------------------------------- *)
(* 4.6  EVALUATE.  gs = 1, so everything is in units of gs^4.            *)
(* -------------------------------------------------------------------- *)

numAssert[expr_, label_] := If[NumericQ[expr] || (ArrayQ[expr] &&
     AllTrue[Flatten[expr], NumericQ]), expr,
   (Print["** ", label, " is NOT numeric. Surviving symbols: ",
          Variables[expr]];
    Print["   -> a replacement rule matched nothing. Run psWhichFree[]."];
    Abort[])];

psEvaluate[kin_Association] := Module[{rule, test, Rn, An, Cn, Bp, Bm},

  (* substitute BOTH; the assigned symbol ignores its rule harmlessly *)
  rule = {mt  -> kin["mtRescaled"],
          bet -> kin["beta"],
          th  -> kin["theta"],
          gs  -> 1};

  test = N[Atil /. rule];
  If[! NumericQ[test],
     Print["** substitution FAILED. Atil still contains: ", Variables[test]];
     psWhichFree[];
     Return[$Failed]];

  Rn = Chop[N[Rhel   /. rule], 10^-12];
  An = Chop[N[Atil   /. rule], 10^-12];
  Cn = Chop[N[Cij    /. rule], 10^-12];
  Bp = Chop[N[Bplus  /. rule], 10^-12];
  Bm = Chop[N[Bminus /. rule], 10^-12];

  numAssert[Rn, "R"]; numAssert[Cn, "Cij"];

  <|"R" -> Rn, "Atil" -> An, "C" -> Cn, "Bplus" -> Bp, "Bminus" -> Bm,
    "trR" -> Chop[Tr[Rn]], "M2ref" -> N[M2ref /. rule]|>
];


(* -------------------------------------------------------------------- *)
(* 4.7  CHECKS, including the azimuth-invariance argument                *)
(* -------------------------------------------------------------------- *)

psChecks[kin_, ev_] := Module[{R, rot, kH0, rH0, nH0, ph, th0, rho, eig},
  R = ev["R"]; ph = kin["phi"]; th0 = kin["theta"];

  kH0 = {Sin[th0], 0., Cos[th0]};
  rH0 = {-Cos[th0], 0., Sin[th0]};
  nH0 = {0., 1., 0.};
  rot = {{Cos[ph], -Sin[ph], 0.}, {Sin[ph], Cos[ph], 0.}, {0., 0., 1.}};

  rho = R/Tr[R];
  eig = Sort[Re[Eigenvalues[rho]], Greater];

  Print["CHECKS"];
  Print["  momentum basis vs rotated phi=0 basis (want ~0):  k ",
        Chop[Max[Abs[rot . kH0 - kin["kHat"]]]], "   r ",
        Chop[Max[Abs[rot . rH0 - kin["rHat"]]]], "   n ",
        Chop[Max[Abs[rot . nH0 - kin["nHat"]]]]];
  Print["  basis orthonormal {1,1,1,0,0,0}: ",
        Chop[{kin["kHat"] . kin["kHat"], kin["rHat"] . kin["rHat"],
              kin["nHat"] . kin["nHat"], kin["kHat"] . kin["rHat"],
              kin["kHat"] . kin["nHat"], kin["rHat"] . kin["nHat"]}]];
  Print["  R Hermitian (want 0):            ",
        Chop[Max[Abs[R - ConjugateTranspose[R]]]]];
  Print["  R real, max|Im| (want 0):        ", Chop[Max[Abs[Im[R]]]]];
  Print["  diagonal >= 0:                   ", Min[Re[Diagonal[R]]] >= -10^-12];
  Print["  Tr R - |M|^2 (want 0):           ", Chop[ev["trR"] - ev["M2ref"]]];
  Print["  4 Atil - Tr R (want 0):          ", Chop[4 ev["Atil"] - ev["trR"]]];
  Print["  B+ (want {0,0,0}):               ", Chop[ev["Bplus"]]];
  Print["  B- (want {0,0,0}):               ", Chop[ev["Bminus"]]];
  Print["  off-plane C_rn,C_nr,C_nk,C_kn:   ",
        Chop[{ev["C"][[1,2]], ev["C"][[2,1]], ev["C"][[2,3]], ev["C"][[3,2]]}]];
  Print["  C_rk - C_kr (want 0):            ", Chop[ev["C"][[1,3]] - ev["C"][[3,1]]]];
  Print["  rho eigenvalues:                 ", Chop[eig]];
  Print["  sum of eigenvalues (want 1):     ", Chop[Total[eig]]];
];


(* -------------------------------------------------------------------- *)
(* 4.8  ENTANGLEMENT.  Concurrence and tr C are BASIS INDEPENDENT, so    *)
(*      they can be quoted before the convention questions are settled.  *)
(* -------------------------------------------------------------------- *)

psEntanglement[ev_] := Module[{R, rho, sy, YY, M, lam, conc},
  R = ev["R"];
  rho = R/Tr[R];
  sy = {{0, -I}, {I, 0}};
  YY = KroneckerProduct[sy, sy];
  M = rho . YY . Conjugate[rho] . YY;
  lam = Sort[Sqrt[Abs[Re[Eigenvalues[M]]]], Greater];
  conc = Max[0, lam[[1]] - lam[[2]] - lam[[3]] - lam[[4]]];
  <|"concurrence" -> conc,
    "trC" -> Tr[ev["C"]]/ev["Atil"],
    "entangled" -> (conc > 10^-10)|>
];


(* -------------------------------------------------------------------- *)
(* 4.9  FORMATTED REPORT - this is what you send                         *)
(* -------------------------------------------------------------------- *)

psReport[mom_List] := Module[{kin, ev, ent, ax, f, v},
  ax = {"r", "n", "k"};
  f[x_]  := NumberForm[N[x], {12, 10}];
  v[x_]  := Row[{"{", Row[f /@ x, ",  "], "}"}];

  kin = psExtract[mom];
  ev  = psEvaluate[kin];
  If[ev === $Failed, Return[$Failed]];
  ent = psEntanglement[ev];

  Print["===================================================================="];
  Print[" LO SM  gg -> t tbar  spin-density matrix at the given phase-space point"];
  Print["===================================================================="];
  Print["KINEMATICS"];
  Print["   sqrt(s)        = ", f[kin["sqrtS"]], " GeV"];
  Print["   m_t            = ", f[kin["mt"]], " GeV"];
  Print["   beta           = ", f[kin["beta"]]];
  Print["   cos(theta)     = ", f[kin["cosTheta"]]];
  Print["   theta          = ", f[kin["theta"]], " rad"];
  Print["   phi (azimuth)  = ", f[kin["phi"] 180/Pi], " deg"];
  Print["   mt/E  (rescaled, = the notebook's mt) = ", f[kin["mtRescaled"]]];
  Print[""];
  Print["SPIN BASIS  (built from the actual momenta; beam along +z)"];
  Print["   kHat = ", v[kin["kHat"]], "   top direction / helicity axis"];
  Print["   rHat = ", v[kin["rHat"]], "   in production plane, perp to k"];
  Print["   nHat = ", v[kin["nHat"]], "   = rHat x kHat, normal to plane"];
  Print[""];
  Print["   CONVENTIONS - STATE THESE WHEN COMPARING:"];
  Print["     nHat = rHat x kHat  ->  (r,n,k) is LEFT-handed"];
  Print["     antitop spin labelled relative to -kHat,"];
  Print["       physical spin operator = -(sigma.a)^*  (charge conjugation)"];
  Print["     both are OPEN questions vs the group's reference code"];
  Print[""];
  Print["DENSITY MATRIX  R / gs^4    basis |++>, |+->, |-+>, |-->  (top first)"];
  Print[MatrixForm[Map[f, Re[ev["R"]], {2}]]];
  Print["   Tr R = ", f[ev["trR"]], "   ( = |M|^2 / gs^4 )"];
  Print[""];
  Print["BLOCH-FANO COEFFICIENTS   (units of gs^4)"];
  Print["   Atil = ", f[ev["Atil"]]];
  Print["   B+   = ", v[ev["Bplus"]],  "   zero at LO (parity)"];
  Print["   B-   = ", v[ev["Bminus"]]];
  Print[""];
  Print[Grid[
    Join[
      {{"", "C_i r", "C_i n", "C_i k", "", "/Atil  r", "/Atil  n", "/Atil  k"}},
      Table[
        Join[{ax[[i]]},
             Table[f[ev["C"][[i, j]]], {j, 3}],
             {""},
             Table[NumberForm[N[ev["C"][[i, j]]/ev["Atil"]], {10, 8}], {j, 3}]],
        {i, 3}]],
    Frame -> All, Alignment -> Right]];
  Print[""];
  Print["DERIVED  (basis independent - safe to quote now)"];
  Print["   tr C / Atil = ", f[ent["trC"]]];
  Print["   concurrence = ", f[ent["concurrence"]],
        If[ent["entangled"], "   -> ENTANGLED", "   -> separable"]];
  Print["===================================================================="];
  Print[""];

  psChecks[kin, ev];

  <|"kin" -> kin, "ev" -> ev, "ent" -> ent|>
];


(* -------------------------------------------------------------------- *)
(* 4.10  SELF-TEST against independently computed reference values       *)
(* -------------------------------------------------------------------- *)

psSelfTest[res_] := Module[{ev, kin, ref, got, dev},
  If[res === $Failed, Print["4.10 self-test skipped"]; Return[$Failed]];
  kin = res["kin"]; ev = res["ev"];
  ref = <|"beta" -> 0.9382345123, "cosTheta" -> -0.1667602549,
          "Atil" -> 0.0475904733, "trR" -> 0.1903618931,
          "Crr" -> 0.0362479415, "Cnn" -> -0.0295085436,
          "Ckk" -> 0.0219278671, "Crk" -> -0.0038479679|>;
  got = <|"beta" -> kin["beta"], "cosTheta" -> kin["cosTheta"],
          "Atil" -> ev["Atil"], "trR" -> ev["trR"],
          "Crr" -> ev["C"][[1,1]], "Cnn" -> ev["C"][[2,2]],
          "Ckk" -> ev["C"][[3,3]], "Crk" -> ev["C"][[1,3]]|>;
  Print["4.10 self-test vs independent reference (tolerance 1e-8)"];
  Do[
    dev = Abs[got[k] - ref[k]];
    Print["   ", StringPadRight[k, 10],
          " got ", NumberForm[N[got[k]], {12,10}],
          "   ref ", NumberForm[N[ref[k]], {12,10}],
          "   ", If[dev < 10^-8, "OK", "** DEVIATION " <> ToString[dev]]],
    {k, Keys[ref]}];
];


(* -------------------------------------------------------------------- *)
(* 4.11  RUN at Marcel's reference point                                 *)
(*                                                                       *)
(*  NOTE the azimuth is 3.66 deg, so the top is NOT in the x-z plane     *)
(*  where the notebook builds its basis.  A rotation about the beam is   *)
(*  a symmetry of the initial state, so the COEFFICIENTS are unchanged - *)
(*  but the BASIS VECTORS are rotated and must be quoted with them.      *)
(* -------------------------------------------------------------------- *)

momnum = {{500.`, 0.`, 0.`, 500.`},
          {500.`, 0.`, 0.`, -500.`},
          {500.`,  461.604563622633`,  29.53432289997799`, -78.23011321478687`},
          {500.`, -461.604563622633`, -29.53432289997799`,  78.23011321478687`}};

psWhichFree[];
result = psReport[momnum];
psSelfTest[result];

(* EXPECTED:                                                             *)
(*   sqrt(s) = 1000.000000   m_t = 173.0000000   beta = 0.9382345123     *)
(*   cos(theta) = -0.1667602549   theta = 1.7383393226   phi = 3.660901  *)
(*   Tr R = 0.1903618931          Atil = 0.0475904733                    *)
(*   B+ = B- = {0, 0, 0}                                                 *)
(*   C_rr/Atil = +0.76166382      C_nn/Atil = -0.62005148                *)
(*   C_kk/Atil = +0.46076169      C_rk/Atil = -0.08085584                *)
(*   tr C / Atil = +0.60237403    concurrence = 0.4212385 -> ENTANGLED   *)
(*   rho eigenvalues = {0.7106192488, 0.1803879332,                      *)
(*                      0.0994064930, 0.0095863249}                      *)
(*                                                                       *)
(* To evaluate elsewhere:   psReport[makePoint[1000., 173., 1.2, 0.]]    *)


(* ==================================================================== *)
(*                     ---  NLO STARTS HERE  ---                        *)
(*                                                                      *)
(*  At one loop the amplitude picks up a correction M(1), and the cross *)
(*  section gets the INTERFERENCE 2 Re( M(1) M(0)* ) - not the square   *)
(*  of the loop.  |M^(1)|^2 is order gs^8, which is NNLO.  So the tree  *)
(*  amplitudes above are still needed.                                  *)
(*                                                                      *)
(*  Each loop diagram has loop momenta in the numerator, giving TENSOR  *)
(*  integrals.  Passarino-Veltman reduction rewrites all of them in     *)
(*  terms of a small fixed basis of SCALAR integrals with 1, 2, 3 or 4  *)
(*  propagators: A0, B0, C0, D0.  That is why "list the one-loop        *)
(*  integrals" is a finite task even though there are ~34 diagrams.     *)
(*                                                                      *)
(*  THE PHYSICS TO WATCH FOR: above threshold the loop integrals        *)
(*  develop IMAGINARY parts (Im B0 = Pi Sqrt[lambda]/s).  Those give    *)
(*  relative phases between amplitudes, which switch on the T-odd       *)
(*  coefficients B+_n, B-_n, C_rn, C_nk that are identically zero at    *)
(*  tree level.  If they come back zero at one loop, something is       *)
(*  wrong.  See PART 10.                                                *)
(* ==================================================================== *)


(* ==================================================================== *)
(*  PART 5  --  NLO CONVENTIONS AND THE ANALYTIC B0                     *)
(*                                                                      *)
(*   CONVENTIONS                                                        *)
(*     B0 = Delta_eps + B0fin,  Delta_eps = 1/eps - EulerGamma + Log[4 Pi] *)
(*     propagators  1/(q^2 - m^2 + I eps)                               *)
(*     => Im B0 > 0 above threshold  (optical theorem)                  *)
(*     B0fin depends on mu^2.  ALWAYS quote mu with a number.           *)
(* ==================================================================== *)

(* -------------------------------------------------------------------- *)
(* 5.1  THE DERIVATION                                                   *)
(*                                                                       *)
(*    B0 = Delta - Int_0^1 dx Log[(p2 x^2 - (p2+m02-m12) x + m02 - Ieps)/mu2] *)
(*                                                                       *)
(*  Factor the quadratic as p2 (x-xp)(x-xm) with                         *)
(*    xpm = ( (p2+m02-m12) +- Sqrt[lambda] ) / (2 p2),                   *)
(*    lambda = (p2+m02-m12)^2 - 4 p2 m02       (Kallen function)         *)
(*  and use  Int_0^1 Log[x-a] dx = (1-a)Log[1-a] + a Log[-a] - 1         *)
(*  to get                                                               *)
(*    B0fin = 2 - Log[p2/mu2]                                            *)
(*            - Sum_i [ (1-x_i) Log[1-x_i] + x_i Log[-x_i] ]             *)
(*                                                                       *)
(*  The I eps matters.  We put -I eps on BOTH the masses and on p2.      *)
(*  Getting this wrong gives the correct magnitude with a wrong sign on  *)
(*  Im B0, or an error of 2 Pi I in the spacelike region.                *)
(* -------------------------------------------------------------------- *)

Clear[lpEps, B0fin, B0finSeries, mHval, mtval, muSq];

lpEps = 10^-20;   (* infinitesimal; small enough not to perturb the answer *)

mHval = 125.25;
mtval = 172.5;
muSq  = mtval^2;          (* <-- STATE THIS WHENEVER YOU QUOTE A NUMBER *)

B0fin[p2_?NumericQ, m02_?NumericQ, m12_?NumericQ, mu2_?NumericQ] :=
 Module[{M0, M1, bb, sq, xp, xm, res},

   M0 = m02 - I lpEps;
   M1 = m12 - I lpEps;

   (* degenerate case p^2 = 0 *)
   If[Abs[p2] < 10^-10,
      Return[
        If[Abs[m02 - m12] < 10^-10,
           -Log[M0/mu2],                                    (* equal masses *)
           1 - (M0 Log[M0/mu2] - M1 Log[M1/mu2])/(M0 - M1)   (* general      *)
        ]]
   ];

   bb = p2 + M0 - M1;
   sq = Sqrt[bb^2 - 4 p2 M0];
   xp = (bb + sq)/(2 p2);
   xm = (bb - sq)/(2 p2);

   res = 2 - Log[(p2 - I lpEps)/mu2];
   res -= (1 - xp) Log[1 - xp] + xp Log[-xp];
   res -= (1 - xm) Log[1 - xm] + xm Log[-xm];
   res
 ];

(* full B0 including the pole, as a series in eps *)
B0finSeries[p2_, m02_, m12_, mu2_, eps_] :=
   1/eps - EulerGamma + Log[4 Pi] + B0fin[p2, m02, m12, mu2];


(* -------------------------------------------------------------------- *)
(* 5.2  CLOSED-FORM SPECIAL CASES -- independent checks                  *)
(* -------------------------------------------------------------------- *)

Clear[B0massless, B0oneMass, B0zeroMom, B0equalMass, B0onShell];

(* B0(p^2, 0, 0) = 2 - Log[p2/mu2] + I Pi   for p2 > 0                   *)
(* Note the argument is -(p2) - I eps, NOT p2 - I eps: for timelike p2   *)
(* the branch point is approached from the negative side, and writing    *)
(* it the other way loses the I Pi entirely.                             *)
B0massless[p2_, mu2_] := 2 - Log[(-p2 - I lpEps)/mu2];

(* B0(p^2, 0, m^2).  Threshold here is just m^2 (one line is massless),  *)
(* so this already has an imaginary part for p2 > m^2.  The eps must be  *)
(* NEGATIVE: with +I eps you land on the wrong side and Im comes out     *)
(* with a flipped sign, i.e. off by 2x the right answer.                 *)
B0oneMass[p2_, m2_, mu2_] :=
   2 - Log[m2/mu2] + (m2/p2 - 1) Log[1 - p2/m2 - I lpEps];

(* B0(0, m0^2, m1^2) *)
B0zeroMom[m02_, m12_, mu2_] :=
   If[Abs[m02 - m12] < 10^-10, -Log[m02/mu2],
      1 - (m02 Log[m02/mu2] - m12 Log[m12/mu2])/(m02 - m12)];

(* B0(p^2, m^2, m^2) -- the top loop in the gluon self-energy, i.e.      *)
(* B0(S; mt, mt).  beta here is the velocity of the pair INSIDE the      *)
(* loop, NOT the top velocity of the outer process.                      *)
B0equalMass[p2_, m2_, mu2_] := Module[{be},
   be = Sqrt[1 - 4 m2/(p2 - I lpEps)];
   2 - Log[m2/mu2] + be Log[(be - 1)/(be + 1)]];

(* B0(m^2, 0, m^2) -- the on-shell top self-energy.  The whole           *)
(* momentum-dependent piece vanishes at p2 = m^2, leaving just a log.    *)
B0onShell[m2_, mu2_] := 2 - Log[m2/mu2];


(* -------------------------------------------------------------------- *)
(* 5.3  NUMERICAL GROUND TRUTH: direct Feynman-parameter integration     *)
(*      Use this to validate ANY loop-integral implementation,           *)
(*      including Collier.                                               *)
(* -------------------------------------------------------------------- *)

Clear[B0numeric];

B0numeric[p2_?NumericQ, m02_?NumericQ, m12_?NumericQ, mu2_?NumericQ,
          eps_: 10^-9] :=
  NIntegrate[
     -Log[(p2 x^2 - (p2 + m02 - m12) x + m02 - I eps)/mu2],
     {x, 0, 1},
     WorkingPrecision -> 20, AccuracyGoal -> 12, MaxRecursion -> 30];


(* -------------------------------------------------------------------- *)
(* 5.4  VALIDATION                                                       *)
(*      Twelve points spanning spacelike, below threshold, on threshold  *)
(*      and above threshold.  Reference values computed independently.   *)
(* -------------------------------------------------------------------- *)

Clear[lpTests, lpCheck];

lpTests = {
 {"B0[62500,  mH^2, mt^2]  <-- TASK 1",  62500.,  mHval^2, mtval^2},
 {"B0[250,    mH^2, mt^2]  <-- alt read",  250.,  mHval^2, mtval^2},
 {"B0[62500,  mt^2, mH^2]  symmetry",    62500.,  mtval^2, mHval^2},
 {"B0[100000, mH^2, mt^2]  above thr",  100000.,  mHval^2, mtval^2},
 {"B0[0,      mH^2, mt^2]",                  0.,  mHval^2, mtval^2},
 {"B0[mt^2,   0,    mt^2]",             mtval^2, 10.^-20, mtval^2},
 {"B0[0,      mt^2, mt^2]",                  0.,  mtval^2, mtval^2},
 {"B0[-50000, 0,    mt^2]  spacelike",  -50000., 10.^-20, mtval^2},
 {"B0[-10^6,  mt^2, mt^2]  spacelike",   -10.^6,  mtval^2, mtval^2},
 {"B0[10^6,   mt^2, mt^2]  s=1TeV^2",     10.^6,  mtval^2, mtval^2},
 {"B0[10^6,   0,    0   ]  massless",     10.^6, 10.^-20, 10.^-20},
 {"B0[29929,  0,    mt^2]  p2=mt^2",     29929., 10.^-20, mtval^2}
};

lpCheck[] := Module[{a, b, bad = 0},
  Print["5.4 B0 validation, mu^2 = mt^2 = ", muSq];
  Print[StringPadRight["  case", 40], StringPadLeft["closed form", 30],
        StringPadLeft["numeric", 30]];
  Do[
    a = B0fin[t[[2]], t[[3]], t[[4]], muSq];
    b = B0numeric[t[[2]], t[[3]], t[[4]], muSq];
    If[Abs[a - b] > 10^-6, bad++];
    Print[StringPadRight[If[Abs[a-b] > 10^-6, "**", "OK"] <> " " <> t[[1]], 40],
          StringPadLeft[ToString[N[a, 10]], 30],
          StringPadLeft[ToString[N[b, 10]], 30]],
    {t, lpTests}];
  Print[If[bad == 0, "  ALL AGREE", ToString[bad] <> " MISMATCHES"]];
];

lpCheck[]

(* Expected (mu^2 = 172.5^2):
     B0[62500,  mH^2, mt^2]  = +0.992443465 + 0.000000000 i   <-- TASK 1b
     B0[250,    mH^2, mt^2]  = +0.288037486 + 0.000000000 i
     B0[100000, mH^2, mt^2]  = +2.038711720 + 1.046280287 i
     B0[0,      mH^2, mt^2]  = +0.286164934 + 0.000000000 i
     B0[mt^2,   0,    mt^2]  = +2.000000000 + 0.000000000 i
     B0[0,      mt^2, mt^2]  =  0           + 0           i
     B0[-50000, 0,    mt^2]  = +0.427308974 + 0.000000000 i
     B0[-10^6,  mt^2, mt^2]  = -1.778323671 + 0.000000000 i
     B0[10^6,   mt^2, mt^2]  = -1.240391193 + 2.948707299 i
     B0[10^6,   0,    0   ]  = -1.514716085 + 3.141592654 i
     B0[29929,  0,    mt^2]  = +2.029719701 + 0.018133253 i
*)


(* -------------------------------------------------------------------- *)
(* 5.5  STRUCTURAL CHECKS that do not rely on the reference numbers      *)
(* -------------------------------------------------------------------- *)

Print["\n5.5 structural checks"];

Print["  symmetry B0(p2,m0,m1) = B0(p2,m1,m0): ",
  Chop[B0fin[62500., mHval^2, mtval^2, muSq]
     - B0fin[62500., mtval^2, mHval^2, muSq]]];

Print["  special case B0(p2,0,0)    vs general: ",
  Chop[B0massless[10.^6, muSq] - B0fin[10.^6, 10.^-20, 10.^-20, muSq]]];

Print["  special case B0(p2,0,m^2)  vs general: ",
  Chop[B0oneMass[62500., mtval^2, muSq] - B0fin[62500., 10.^-20, mtval^2, muSq]]];

Print["  special case B0(0,m0,m1)   vs general: ",
  Chop[B0zeroMom[mHval^2, mtval^2, muSq] - B0fin[0., mHval^2, mtval^2, muSq]]];

Print["  special case B0(p2,m,m)    vs general: ",
  Chop[B0equalMass[10.^6, mtval^2, muSq] - B0fin[10.^6, mtval^2, mtval^2, muSq]]];

Print["  special case B0(m^2,0,m^2) vs general: ",
  Chop[B0onShell[mtval^2, muSq] - B0fin[mtval^2, 10.^-20, mtval^2, muSq]]];

Print["  mu-independence dB0/dLog[mu^2] (want +1): ",
  Chop[(B0fin[62500., mHval^2, mtval^2, 2 muSq]
      - B0fin[62500., mHval^2, mtval^2, muSq])/Log[2]]];

Print["  Im B0 = 0 below threshold (mH+mt)^2 = ", (mHval + mtval)^2, " : ",
  Chop[Im[B0fin[62500., mHval^2, mtval^2, muSq]]]];

Print["  Im B0 > 0 above threshold: ",
  Im[B0fin[100000., mHval^2, mtval^2, muSq]] > 0];

(* the absorptive part analytically:                                     *)
(*   Im B0(s,m0^2,m1^2) = Pi Sqrt[lambda]/s   for s > (m0+m1)^2          *)
Print["  Im B0 vs Pi Sqrt[lambda]/s: ",
 Module[{s = 100000., a = mHval^2, b = mtval^2, lam},
   lam = (s - a - b)^2 - 4 a b;
   Chop[Im[B0fin[s, a, b, muSq]] - Pi Sqrt[lam]/s]]];


(* -------------------------------------------------------------------- *)
(* 5.6  COMPILED VERSION, for phase-space scans                          *)
(* -------------------------------------------------------------------- *)

Clear[B0compiled];

B0compiled = Compile[
  {{p2, _Real}, {m02, _Real}, {m12, _Real}, {mu2, _Real}},
  Module[{ep = 10.^-20, M0, M1, bb, sq, xp, xm, res},
    M0 = m02 - I ep;  M1 = m12 - I ep;
    If[Abs[p2] < 10.^-10,
       If[Abs[m02 - m12] < 10.^-10,
          -Log[M0/mu2],
          1 - (M0 Log[M0/mu2] - M1 Log[M1/mu2])/(M0 - M1)],
       bb = p2 + M0 - M1;
       sq = Sqrt[bb^2 - 4 p2 M0];
       xp = (bb + sq)/(2 p2);  xm = (bb - sq)/(2 p2);
       res = 2 - Log[(p2 - I ep)/mu2];
       res = res - ((1 - xp) Log[1 - xp] + xp Log[-xp]);
       res = res - ((1 - xm) Log[1 - xm] + xm Log[-xm]);
       res
    ]],
  CompilationTarget -> "C",          (* drop this if no C compiler *)
  RuntimeOptions -> "Speed"
];

Print["\n5.6 compiled vs interpreted: ",
  Chop[B0compiled[62500., mHval^2, mtval^2, muSq]
     - B0fin[62500., mHval^2, mtval^2, muSq]]];


(* ==================================================================== *)
(*  PART 6  --  THE COLLIER INTERFACE                                   *)
(*                                                                      *)
(*  The wrapper lives in the repository and is loaded, not duplicated    *)
(*  here, because it has to know how to launch loopdrv (an external      *)
(*  Fortran program living in WSL).                                     *)
(*                                                                      *)
(*  It provides, with memoisation:                                      *)
(*     B0cF[p10, m02, m12, mu2]  -> <|"fin", "c1", "c2"|>               *)
(*     B0c [p10, m02, m12, mu2]  -> the finite part only                *)
(*     likewise A0cF/A0c, C0cF/C0c, D0cF/D0c                            *)
(*  where  I = c2/eps^2 + c1/eps + fin.                                 *)
(*                                                                      *)
(*  ALWAYS use the cached versions.  Each scalar integral depends only  *)
(*  on the phase-space point, NOT on the helicity configuration, so the *)
(*  same number is asked for sixteen times.  Without caching each one   *)
(*  spawns a process.  This is the same lesson as the LO memoisation.   *)
(*                                                                      *)
(*  WHAT THE VALIDATION RESTS ON                                        *)
(*    B0 finite      analytic derivation (PART 5)  INDEPENDENT  15 pts  *)
(*    B0 pole        c1 = 1 by construction        INDEPENDENT  15 pts  *)
(*    C0, D0 finite  Feynman-param quadrature      INDEPENDENT  10 pts  *)
(*                   (below threshold / spacelike only - above it the   *)
(*                    integrand has a pole inside the region)           *)
(*    C0, D0 above   COLI vs DD                    internal to Collier  *)
(*                                                                      *)
(*  POLE STRUCTURE FOUND (gg -> ttbar kinematics)                       *)
(*    B0, every case           c1 = 1 exactly, c2 = 0                   *)
(*    A0(m^2)                  c1 = m^2,       c2 = 0                   *)
(*    C0 all lines massive     finite                                   *)
(*    C0 one massless line     c1 only                                  *)
(*    C0 two massless lines    c1 and c2                                *)
(*    D0 two massless gluons   c1 and c2    <- the IR-divergent box     *)
(*    D0 one massless line     FINITE                                   *)
(*  and every c2 is purely REAL: the leading soft-collinear pole        *)
(*  carries no phase.  The mass content PREDICTS the pole structure and *)
(*  it matched in every row - a real consistency check in itself.       *)
(* ==================================================================== *)

Clear[wsl, colFile, s2File];

wsl     = "\\\\wsl.localhost\\Ubuntu-22.04\\home\\vaibhav\\collier_work\\";
colFile = wsl <> "collier-interface\\mathematica\\collier_mathematica.m";
s2File  = wsl <> "formcalc_stage2_pave.mx";

(* PATH WARNING: use backslashes CONSISTENTLY.  Mixed separators pass    *)
(* FileExistsQ but fail in Get and $Path - that is what produced         *)
(* Get::noopen when FeynArts was first loaded.                           *)

Print["\nPART 6  Collier interface"];
Print["  ", If[FileExistsQ[colFile], "found   ", "MISSING "], colFile];

If[FileExistsQ[colFile], Get[colFile],
   Print["  ** not loaded - PARTS 6-10 that need numbers will not run"]];


(* ==================================================================== *)
(*  PART 7  --  THE ONE-LOOP INTEGRAL LIST          (TASK 1a - DONE)    *)
(*                                                                      *)
(*  Produced in WSL by ttbar_stages.wls.  FormCalc CANNOT run in        *)
(*  Windows Mathematica: MathLink will not launch ReadForm across the   *)
(*  Windows/WSL boundary.  Reading the .mx it writes IS fine - that is  *)
(*  a file read, not a process launch.                                  *)
(*                                                                      *)
(*  THE OPTION THAT MATTERS                                             *)
(*      CalcFeynAmp[..., FermionChains -> VA, PaVeReduce -> True]       *)
(*                                                                      *)
(*  Without PaVeReduce the same run returns 222 distinct objects with   *)
(*  38 coefficient types (bb0, bb1, bb00, cc00, dd123 ...) - tensor     *)
(*  coefficients, which would need Collier's tensor routines.  With it  *)
(*  everything reduces to bb0 / cc0 / dd0, scalars only, which the      *)
(*  existing driver already computes and the existing reference points  *)
(*  already validate.  (Marcel's correction.)                           *)
(*                                                                      *)
(*  THE DIAGRAM SET                                                     *)
(*    tree       3 diagrams (s, t, u channel)                           *)
(*    one loop  34 QCD diagrams                                         *)
(*                                                                      *)
(*  SMQCD inherits the full SM content, so the unfiltered run gives     *)
(*  112 diagrams including H, G0, G+-, W, Z, photon and b - those are   *)
(*  O(alpha), not O(alpha_s).  ExcludeParticles removes them.           *)
(*                                                                      *)
(*  The 34 include GHOST loops and LIGHT-QUARK loops.  The hand         *)
(*  catalogue had ~15 integrals and missed both entirely.               *)
(* ==================================================================== *)

Print["\nPART 7  the integral list"];
Print["  ", If[FileExistsQ[s2File], "found   ", "MISSING "], s2File];

If[FileExistsQ[s2File], Get[s2File],
   Print["  ** run ttbar_stages.wls in WSL first"]];

(* --- 7.1  how to read an entry --------------------------------------- *)
(*                                                                       *)
(*     B0i[bb0, T, 0, MT2]                                               *)
(*      |    |   |  |   `-- second internal mass^2                       *)
(*      |    |   |  `------ first internal mass^2                        *)
(*      |    |   `--------- external p^2 (here the Mandelstam t)         *)
(*      |    `------------- coefficient type; bb0 = the scalar B0        *)
(*      `------------------ FormCalc's B head                            *)
(*                                                                       *)
(*   so that is  B0(t; 0, mt^2)  ->  ./loopdrv B0 T 0 MT2 mu2            *)
(*                                                                       *)
(*  Mf2[3,Gen5] and Mf2[4,Gen5] are the light-quark masses running in    *)
(*  the gluon self-energy.  Set them to zero for massless quarks and     *)
(*  several entries collapse onto ones that already exist.               *)
(* --------------------------------------------------------------------- *)

(* The symbols in the .mx live in the FormCalc` context, not Global`, so  *)
(* a pattern like Count[ints, _A0] silently matches nothing.  Match on    *)
(* the symbol NAME instead.                                               *)

Clear[headName, intCount];
headName[x_] := SymbolName[Head[x]];
intCount[nm_] := Count[ints, _?(headName[#] === nm &)];

If[Head[ints] === List,
   Print["  ", Length[ints], " distinct scalar integrals"];
   Print[Grid[{{"A0", intCount["A0"]},   {"B0", intCount["B0i"]},
               {"C0", intCount["C0i"]},  {"D0", intCount["D0i"]}},
              Frame -> All, Alignment -> Left]];
   Print["  expected: A0 3, B0 11, C0 14, D0 6"]];


(* ==================================================================== *)
(*  PART 8  --  THE COLOUR MATRIX ON THE EXTENDED BASIS                 *)
(*                                                                      *)
(*  FormCalc returned these colour structures at one loop:              *)
(*                                                                      *)
(*     SUNT[Glu1, Glu2, Col3, Col4]   =  T^a T^b                        *)
(*     SUNT[Glu2, Glu1, Col3, Col4]   =  T^b T^a                        *)
(*     SUNT[Col3, Col4]               =  delta_ij        <- NEW         *)
(*     SUNT[Glu1, Glu2, 0, 0]         =  Tr(T^a T^b), a closed loop     *)
(*                                                                      *)
(*  THE POINT.  At tree level the basis was c1 = T^aT^b, c2 = T^bT^a,   *)
(*  c3 = f^abc T^c.  But                                                *)
(*                                                                      *)
(*      f^abc T^c = -i [T^a, T^b] = -i (c1 - c2)                        *)
(*                                                                      *)
(*  so as MATRICES the tree basis spans only a 2-dimensional space.     *)
(*                                                                      *)
(*  At one loop delta^ab * Identity appears, and that IS genuinely      *)
(*  independent.  So the natural basis becomes                          *)
(*                                                                      *)
(*     c1 = T^a T^b,   c2 = T^b T^a,   c3 = delta^ab * Identity_3       *)
(*                                                                      *)
(*  and this one IS 3-dimensional.  Everything else is expressible:     *)
(*                                                                      *)
(*     f^abc T^c = -i (c1 - c2)                                         *)
(*     d^abc T^c = c1 + c2 - (1/N) c3                                   *)
(*                                                                      *)
(*  so {c1, c2, c3} is complete for gg -> q qbar.                       *)
(* ==================================================================== *)

Clear[Nc, CF, CA, sumTrTT, lemma, CcolorNLO, colourNLO, toNewBasis];

Nc = 3;
CF = (Nc^2 - 1)/(2 Nc);          (* 4/3  - quark colour charge  *)
CA = Nc;                          (* 3    - gluon colour charge  *)

(* --- 8.1  the identities needed -------------------------------------- *)
(*                                                                       *)
(*   Tr(T^a T^b)       = (1/2) delta^ab                                  *)
(*   sum_a T^a T^a     = CF * Identity                                   *)
(*   sum_a Tr(T^a T^a) = (1/2) delta^aa = (N^2-1)/2 = 4                  *)
(*   T^a T^b T^a       = (CF - CA/2) T^b = -1/6 T^b                      *)
(*                                                                       *)
(*  That last one: move the last generator past T^b with the commutator, *)
(*  then use cyclicity and antisymmetry to get f^abc f^acd = -CA delta^bd*)
(* --------------------------------------------------------------------- *)

sumTrTT = (Nc^2 - 1)/2;           (* = 4    *)
lemma   = CF - CA/2;              (* = -1/6 *)

(* --- 8.2  the entries, derived --------------------------------------- *)
(*                                                                       *)
(*  C_mn = sum over a, b and quark colour of  c_m (c_n)^dagger           *)
(*                                                                       *)
(*  C11 = sum_ab Tr(T^aT^b T^bT^a) = CF sum_a Tr(T^aT^a) = (4/3)(4)      *)
(*  C22 = C11                               (relabel a <-> b)            *)
(*  C12 = sum_ab Tr(T^aT^b T^aT^b)                                       *)
(*      = sum_b Tr( [T^aT^bT^a] T^b ) = (CF-CA/2) sum_b Tr(T^bT^b)       *)
(*  C13 = sum_ab delta^ab Tr(T^aT^b) = sum_a Tr(T^aT^a)                  *)
(*  C23 = C13                               (same trace)                 *)
(*  C33 = sum_ab delta^ab delta^ab Tr(1) = (N^2-1) * N                   *)
(* --------------------------------------------------------------------- *)

CcolorNLO = {
  {CF sumTrTT,     lemma sumTrTT,  sumTrTT     },
  {lemma sumTrTT,  CF sumTrTT,     sumTrTT     },
  {sumTrTT,        sumTrTT,        (Nc^2-1) Nc }
};

colourNLO[a_List, b_List] := a . CcolorNLO . Conjugate[b];

Print["\nPART 8  colour matrix on the extended basis"];
Print["  basis: c1 = T^a T^b,  c2 = T^b T^a,  c3 = delta^ab * 1"];
Print[MatrixForm[CcolorNLO]];
Print["  symmetric  : ", CcolorNLO === Transpose[CcolorNLO]];
Print["  eigenvalues: ", N[Eigenvalues[CcolorNLO], 6],
      "   (all > 0, as a Gram matrix must be)"];

(* --- 8.3  consistency with the tree matrix --------------------------- *)
(*                                                                       *)
(*  The 1-1 and 1-2 entries must be UNCHANGED from the tree calculation, *)
(*  because c1 and c2 are the same objects.  Only the third row and      *)
(*  column are new.  That is a free check on the new entries.            *)
(* --------------------------------------------------------------------- *)

Print["  C11 = 16/3 ? ", CcolorNLO[[1,1]] === 16/3,
      "    C12 = -2/3 ? ", CcolorNLO[[1,2]] === -2/3];
Print["  C13 = ", CcolorNLO[[1,3]], ",  C33 = ", CcolorNLO[[3,3]],
      "   (new - no tree counterpart)"];


(* --- 8.4  CONVERT THE TREE AMPLITUDES TO THE NEW BASIS  (step N1) ---- *)
(*                                                                       *)
(*  The LO code stores A = (A_t, A_u, A_s) against c3^old = f^abc T^c.   *)
(*  Since c3^old = -i (c1 - c2),                                         *)
(*                                                                       *)
(*    A_t c1 + A_u c2 + A_s c3^old                                       *)
(*       = (A_t - i A_s) c1 + (A_u + i A_s) c2 + 0 . c3                  *)
(*                                                                       *)
(*  THE CHECK IS THE POINT.  With the new 3x3 matrix, colourNLO applied  *)
(*  to the converted vectors must reproduce the old colorProduct for     *)
(*  EVERY helicity pair - hence the same Tr R and the same five          *)
(*  coefficients.  Do this before any loop work: a sign error here       *)
(*  would silently corrupt every interference term later.                *)
(*                                                                       *)
(*  The opposite sign convention fails this check by O(1), so the test   *)
(*  genuinely discriminates.                                             *)
(* --------------------------------------------------------------------- *)

toNewBasis[A_List] := {A[[1]] - I A[[3]], A[[2]] + I A[[3]], 0};

Clear[n1Check];

n1Check[b0_, th0_] := Module[{rule, hs, devs},
  rule = {bet -> b0, th -> th0, gs -> 1};
  hs = Tuples[{helicities, helicities, helicities, helicities}];
  devs = Table[
    Module[{a = N[Avec @@ h1 /. rule], b = N[Avec @@ h2 /. rule]},
      Abs[colorProduct[a, b] - colourNLO[toNewBasis[a], toNewBasis[b]]]],
    {h1, hs}, {h2, hs}];
  Max[devs]];

Print["\n8.4 basis conversion check (want ~1e-15)"];
Print["   bet=0.86930, th=1.1 : ", n1Check[0.86930, 1.1]];
Print["   bet=0.63938, th=2.3 : ", n1Check[0.63938, 2.3]];
Print["   ** if either is O(1), the sign convention on c3 is wrong **"];


(* ==================================================================== *)
(*  PART 9  --  ASSEMBLING R AT ONE LOOP                                *)
(*                                                                      *)
(*      R_ij  ~  sum_lambda [ A0^T C A0*  +  2 Re( A1^T C A0* ) ]       *)
(*                                                                      *)
(*  NOTE this is an INTERFERENCE, not the square of the loop.           *)
(*  |M^(1)|^2 is order gs^8, which is NNLO.  Orders: R^(0) ~ gs^4,      *)
(*  R^(1) ~ gs^6.                                                       *)
(*                                                                      *)
(*  Everything downstream is ALREADY BUILT and unchanged:               *)
(*     the projectors SigR, SigN, SigK      (PART 3.3)                  *)
(*     the {r,n,k} basis                     (PART 1.2)                 *)
(*     the Bloch-Fano extraction             (PART 3.3)                 *)
(*     the phase-space machinery             (PART 4)                   *)
(*                                                                      *)
(*  You are adding ONE more amplitude vector per helicity configuration.*)
(* ==================================================================== *)

(* STILL TO DO - in order:                                              *)
(*                                                                      *)
(*  9a. UNDERSTAND THE FORMCALC OUTPUT STRUCTURE.  Not yet inspected    *)
(*      in detail.  Print Short[amp[[1]], 30], the Mat[...] objects,    *)
(*      how SUNT enters (inside Mat, in coefficients, or in Abbr[]),    *)
(*      what Den[...] and SumOver are, and what the ~331 abbreviations  *)
(*      contain.  FormCalc also has built-in colour/helicity routines   *)
(*      (look for Names["FormCalc`*ME"], e.g. ColourME, HelicityME,     *)
(*      SquaredME).  Three routes, not yet compared in practice:        *)
(*                                                                      *)
(*        (a) own evaluator - map each Mat[] Dirac chain and every      *)
(*            Pair[e[i],k[j]] contraction onto the explicit spinors,    *)
(*            polarisations and gammas of PART 1, keeping top/antitop   *)
(*            helicity labels OPEN.  Keeps the PART 1-3 conventions     *)
(*            intact, which matters because the RELATIVE PHASES of the  *)
(*            helicity amplitudes are what R is made of.                *)
(*        (b) FormCalc's own helicity machinery - but extracting        *)
(*            individual open-spin amplitudes rather than the           *)
(*            spin-summed |M|^2 needs modification, and its phase       *)
(*            conventions must be matched to ours.                      *)
(*        (c) trace projection with spin vectors instead of helicity    *)
(*            amplitudes - avoids spinors, but needs the antitop        *)
(*            projector handled with the same care as PART 3.3.         *)
(*                                                                      *)
(*      (a) fits the existing code best.                                *)
(*                                                                      *)
(*      VALIDATION GATE: the FormCalc TREE amplitude, evaluated with    *)
(*      the chosen method, must reproduce the PART 3 R matrix and all   *)
(*      five coefficients at the reference point.  If it does not,      *)
(*      STOP and fix conventions before touching the loop.              *)
(*                                                                      *)
(*  9b. FORMCALC -> COLLIER ARGUMENT TRANSLATION.  The two use          *)
(*      DIFFERENT argument conventions (FormCalc follows LoopTools;     *)
(*      Collier's C0/D0 use p_ij-type squared momentum differences).    *)
(*      A wrong permutation gives a plausible but WRONG number.  This   *)
(*      is the least-checked link in the whole chain.                   *)
(*                                                                      *)
(*      Test it by: evaluating in a limit with a known closed form;     *)
(*      checking invariance under the diagram's symmetries; and         *)
(*      checking the POLE STRUCTURE (the D0 entries with two massless   *)
(*      lines must show c2 != 0, the one-massless-line ones must be     *)
(*      finite).                                                        *)
(*                                                                      *)
(*      Invariants at the reference point:                              *)
(*        S = 10^6,  MT2 = 29929                                        *)
(*        T = MT2 - (S/2)(1 - bet Cos[th]) ~ -548301                    *)
(*        U = MT2 - (S/2)(1 + bet Cos[th]) ~ -391841                    *)
(*        check  S + T + U = 2 MT2                                      *)
(*                                                                      *)
(*  9c. RENORMALISE.  Counterterms for alpha_s, the top mass and wave   *)
(*      function, and the gluon field.  FormCalc can generate them:     *)
(*                                                                      *)
(*         insCT = InsertFields[CreateCTTopologies[1, 2 -> 2], ...]     *)
(*         ampCT = CalcFeynAmp[CreateFeynAmp[insCT]]                    *)
(*                                                                      *)
(*      The model prints "$CounterTerms are ON; 121 counterterms of     *)
(*      order 1", so the vertices exist - but the renormalisation       *)
(*      CONSTANTS still need values.  Look up their names in            *)
(*      Models/SMQCD.mod rather than guessing.                          *)
(*                                                                      *)
(*      CHECK (UV): the UV poles must cancel between ampLoop and ampCT. *)
(*      The driver returns c1 and c2 separately precisely for this.  A  *)
(*      leftover UV pole means a missing counterterm or wrong scheme.   *)
(*                                                                      *)
(*      CHECK (IR): after renormalisation the remaining 1/eps^2, 1/eps  *)
(*      poles should match the universal IR pole structure of one-loop  *)
(*      QCD amplitudes (Catani's formula, which involves               *)
(*      colour-correlated tree amplitudes).  That tests the colour      *)
(*      matrix, the boxes and the argument translation TOGETHER.        *)
(*      ** I am recalling this formula from memory - verify it in the   *)
(*         literature before relying on it. **                          *)
(*                                                                      *)
(*  9d. BUILD A^(1) per helicity configuration (16 of them) in the new  *)
(*      colour basis, keeping (fin, 1/eps, 1/eps^2) as SEPARATE numbers.*)
(*                                                                      *)
(*  9e. contract with CcolorNLO and project with SigR, SigN, SigK.      *)


(* ==================================================================== *)
(*  PART 10  --  THE PHYSICS CHECK                                      *)
(*                                                                      *)
(*  At tree level every helicity amplitude shares a common phase, which *)
(*  forces                                                              *)
(*                                                                      *)
(*      B+_i = B-_i = 0,   C_rn = C_nr = C_nk = C_kn = 0                *)
(*                                                                      *)
(*  because those are T-odd and need an absorptive part.                *)
(*                                                                      *)
(*  At one loop the integrals develop imaginary parts above threshold,  *)
(*                                                                      *)
(*      Im B0(s, m0^2, m1^2) = Pi Sqrt[lambda] / s                      *)
(*                                                                      *)
(*  (verified in PART 5.5), and those generate relative phases between  *)
(*  interfering amplitudes.  So the T-odd coefficients SWITCH ON for    *)
(*  the first time.                                                     *)
(*                                                                      *)
(*  ** IF EVERYTHING T-ODD COMES BACK EXACTLY ZERO, SOMETHING IS        *)
(*     WRONG IN THE CHAIN. **                                           *)
(*                                                                      *)
(*  Caveat: INDIVIDUAL coefficients may still vanish by symmetry at     *)
(*  special angles or for particular combinations.  Compare against the *)
(*  literature on NLO QCD spin correlations in t tbar (work by          *)
(*  Bernreuther and collaborators) rather than assuming all must be     *)
(*  nonzero.  ** That citation is from memory - check it. **            *)
(*                                                                      *)
(*  Also re-check at one loop:                                          *)
(*    - Hermiticity of R                                                *)
(*    - positivity (eigenvalues >= 0).  NLO VIRTUAL-ONLY can violate    *)
(*      this near IR-dominated regions, so interpret with care.         *)
(*    - Tr R against a known NLO virtual result, if one can be had.     *)
(*                                                                      *)
(*  Report with a STATED mu and scheme.                                 *)
(*                                                                      *)
(*  Physically these T-odd terms are the Standard Model background to   *)
(*  any CP-violating BSM signal - which is the real reason NLO matters  *)
(*  here, rather than mere precision.                                   *)
(* ==================================================================== *)

(*
Print["T-odd coefficients at one loop (expect NON-zero):"];
Print["  B+_n = ", BplusNLO[[2]]];
Print["  B-_n = ", BminusNLO[[2]]];
Print["  C_rn = ", CijNLO[[1, 2]]];
Print["  C_nk = ", CijNLO[[2, 3]]];
*)


(* ==================================================================== *)
(*  OPEN CONVENTIONS  (gate the BSM stage; invisible at LO in the SM)   *)
(*                                                                      *)
(*  1. nHat = rHat x kHat = (0,1,0) makes (r,n,k) LEFT-handed.  The     *)
(*     opposite choice flips the sign of CP-odd observables.            *)
(*  2. The antitop spin axis is labelled relative to -kHat with the     *)
(*     operator -(sigma.a)^*.                                           *)
(*                                                                      *)
(*  Settle both against the group's reference code                      *)
(*  (bernreuther_rmatrix.py) and record them in the repository.         *)
(* ==================================================================== *)

(* ==================================================================== *)
(*  OPEN QUESTIONS FOR MARCEL                                           *)
(*                                                                      *)
(*  1. EFT or the full model?  He mentioned "the 1 loop SM model file   *)
(*     in an EFT file".  A four-quark operator introduces different     *)
(*     colour structures than a propagating coloured scalar, so this    *)
(*     bears directly on PART 8.                                        *)
(*                                                                      *)
(*  2. What does "the density matrix at one loop" mean concretely?      *)
(*     The boxes carry 1/eps^2 and 1/eps IR poles and we have the       *)
(*     coefficients.  Quote them alongside the finite part, or give the *)
(*     finite part in a stated subtraction scheme?  Is real emission    *)
(*     (gg -> t tbar g, and the qg, q qbar channels) in scope?  A       *)
(*     virtual-only R is not an infrared-safe observable.               *)
(*                                                                      *)
(*  3. Which renormalisation scheme for the top mass (and alpha_s)?     *)
(*     It has to match whatever the rate-side analysis uses.            *)
(*                                                                      *)
(*  4. The two conventions above.                                       *)
(*                                                                      *)
(*  5. Spinor-helicity formalism - he asked whether you know it.        *)
(*     (Honest answer: read a bit, not used.)  May be the intended      *)
(*     route for amplitude evaluation.                                  *)
(*                                                                      *)
(*  6. Which light-quark treatment (massless n_f = 5?  4 + top?) and    *)
(*     how the top loop in the gluon self-energy is handled             *)
(*     (decoupling or not).                                             *)
(* ==================================================================== *)

(* ::Package:: *)

(* ==================================================================== *)
(*   gg -> t tbar :  spin-density matrix at LO, and the start of NLO    *)
(*                                                                      *)
(*   CONTENTS                                                           *)
(*     0 - 7    setup: kinematics, spin basis, polarisations, Dirac      *)
(*              algebra, spinors, amplitudes, colour                    *)
(*     7b       Ward identity  (gauge-invariance check)                 *)
(*     8        the trace test   Tr R = |M|^2                           *)
(*     9 - 11   the 4x4 matrix, Bloch-Fano projection, target match     *)
(*     12 - 13  numerical spot checks, conversion to physical units     *)
(*     14 - 15  NLO: one-loop integral catalogue, analytic B0           *)
(*     16       evaluate everything at a given phase-space point        *)
(*     17       Feynman diagrams (LO and NLO), drawn with Graphics      *)
(*                                                                      *)
(*   WHAT THE CALCULATION DOES, IN ONE PARAGRAPH                        *)
(*   We want the spin state of the t tbar pair, not just the rate.  So  *)
(*   instead of summing over final spins (which gives one number) we    *)
(*   keep all 16 helicity amplitudes explicitly, sum and average ONLY   *)
(*   over the initial gluon helicities and colours, and leave the two   *)
(*   top spin labels open.  Those become the row and column of a 4x4    *)
(*   matrix R.  Its trace is the ordinary |M|^2; the other fifteen      *)
(*   numbers are new information that a cross section cannot see.       *)
(*                                                                      *)
(*   DESIGN CHOICES THAT MATTER                                         *)
(*     - $Assumptions declares th, mt, gs REAL.  Without this every     *)
(*       uBar carries Conjugate[th], Conjugate[Sqrt[...]] and           *)
(*       FullSimplify wastes its time unwinding them.  gs especially:   *)
(*       colorProduct takes a Conjugate, so without gs real the trace   *)
(*       comes out as gs^2 Conjugate[gs]^2 and never matches gs^4.      *)
(*     - Everything expensive is MEMOIZED:  f[x_] := f[x] = ...         *)
(*       Building R calls the amplitude routine 128 times; without      *)
(*       caching each call re-runs FullSimplify.  This alone took the   *)
(*       notebook from "cloud timeout" to 1.6 seconds.                  *)
(*     - Energy is scaled out (Et = 1).  A 2->2 amplitude has mass      *)
(*       dimension 4 - 4 = 0, so this is EXACT, and it removes Sqrt[s]  *)
(*       from every spinor and propagator.                              *)
(*     - The TRACE is computed before the 4x4 matrix: it needs only     *)
(*       four diagonal sums and is the decisive normalisation test.     *)
(*     - R is built element by element (upper triangle only), not as    *)
(*       one giant Table plus a global FullSimplify.                    *)
(*     - The spin-projection operators are NOT sigma_1,2,3.  See the    *)
(*       derivation in section 10.                                      *)
(*                                                                      *)
(*   ORDER OF EVALUATION                                                *)
(*   Run Quit[] once, then evaluate sections 0-7 as a single block.     *)
(*   Do not leave long gaps between cells: a recycled Wolfram Cloud     *)
(*   session loses every definition, and the failure is SILENT - you    *)
(*   get plausible-looking output full of unevaluated symbols.  The     *)
(*   sanity gate after section 7 catches this in two seconds.           *)
(* ==================================================================== *)

Quit[]   (* run this first, then evaluate everything below in one go *)


(* ==================================================================== *)
(* 0.  ASSUMPTIONS                                                      *)
(*                                                                      *)
(* Mathematica treats every symbol as COMPLEX unless told otherwise.    *)
(* uBar[h] = Conjugate[uSpin[h]].gam0 therefore generates Conjugate[th],*)
(* Conjugate[Sqrt[1+mt]] and so on, and those survive into every        *)
(* amplitude.  Declaring the symbols real here is the single biggest    *)
(* speed improvement in the whole notebook.  Note gs in particular:     *)
(* colorProduct takes Conjugate[b], so without gs real the answer comes *)
(* out as gs^2 Conjugate[gs]^2 and will NOT match gs^4 |M|^2.           *)
(* ==================================================================== *)

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
(*   FullSimplify cannot see through.  The spinor residuals in section 5 *)
(*   then print as long expressions that are secretly zero.  Check them  *)
(*   numerically instead (see 5c).                                       *)

$Assumptions = assumptions;

simp[x_]  := FullSimplify[x, assumptions];
simpN[x_] := Simplify[x, assumptions];   (* cheaper; use inside hot loops *)


(* ==================================================================== *)
(* 1.  KINEMATICS  (partonic CM frame = t tbar rest frame at LO)        *)
(*                                                                      *)
(* A 2->2 tree amplitude is dimensionless, so we may set the CM energy  *)
(* per particle to 1.  That means                                       *)
(*        Et = Sqrt[s]/2 = 1   <=>   s = 4 ,   mt = Sqrt[1 - bet^2]     *)
(* This is an EXACT rescaling, not an approximation: R comes out        *)
(* identical to the version carrying s and mt explicitly.  The payoff   *)
(* is that Sqrt[s] disappears from every spinor and every propagator,   *)
(* leaving only bet and th.  Restore s at the end if you want it.       *)
(* ==================================================================== *)

Clear[Et, sMand, mt, bet, eta, mdot, k1v, k2v, ptv, ptbv];

(* --------------------------------------------------------------------- *)
(*  WHY Et = 1 AND sMand = 4   (the "1 and 4")                            *)
(*                                                                        *)
(*  In 4 dimensions an amplitude with n external legs has mass dimension  *)
(*  4 - n.  For a 2 -> 2 process n = 4, so the amplitude is EXACTLY       *)
(*  dimensionless - it cannot depend on the overall energy scale, only on *)
(*  ratios.  (Check it directly: the two spinors give Sqrt[E].Sqrt[E],    *)
(*  the propagator numerator gives E, the denominator E^2, net E^0.)      *)
(*                                                                        *)
(*  So we are free to fix the CM energy per particle to 1:                *)
(*      Et = Sqrt[s]/2 = 1   =>   s = (2 Et)^2 = 4                        *)
(*  This is EXACT, not an approximation - verified entry by entry against *)
(*  the version carrying s and mt explicitly.  The payoff is that Sqrt[s] *)
(*  disappears from every spinor normalisation and every propagator, so   *)
(*  everything depends on bet and th alone.                               *)
(*                                                                        *)
(*  In these units mt is the RATIO 2 m_t^phys / Sqrt[s^phys], in (0,1):   *)
(*      mt = Sqrt[1 - bet^2],   mt -> 1 at threshold, 0 at high energy.   *)
(*  Section 13 converts back to physical s and m_t.                       *)
(* --------------------------------------------------------------------- *)

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

(* ---- checks: on-shellness and momentum conservation ---- *)
kinCheck = simp[{
   mdot[k1v, k1v], mdot[k2v, k2v],
   mdot[ptv, ptv] - mt^2, mdot[ptbv, ptbv] - mt^2,
   k1v + k2v - ptv - ptbv
}];
Print["1. kinematics {k1^2, k2^2, pt^2-mt^2, ptb^2-mt^2, cons} = ", kinCheck];


(* ==================================================================== *)
(* 2.  SPIN BASIS  {r, n, k}                                            *)
(*                                                                      *)
(* We need three axes to write the spin correlations against.  Any      *)
(* orthonormal triad works mathematically; this one is chosen so the    *)
(* zeros are obvious.  k = top direction, r = the other in-plane        *)
(* direction, n = normal to the plane containing the beam and the top.  *)
(* At LO every coefficient with an ODD number of n indices vanishes     *)
(* (they are T-odd and need an absorptive part), so this split makes    *)
(* four free checks visible.  In a random basis they would be smeared   *)
(* across many coefficients and the checks would be lost.               *)
(*                                                                      *)
(* k = top direction, r = in the production plane perpendicular to k,   *)
(* n = normal to the production plane.                                  *)
(*                                                                      *)
(* CONVENTION WARNING.  n = r x k gives n = (0,1,0); n = k x r gives    *)
(* the opposite.  At LO in the SM this choice is INVISIBLE, because     *)
(* C_rn = C_nr = C_nk = C_kn = 0 and C_nn is even in n.  It becomes     *)
(* visible only for CP-odd observables, i.e. exactly the BSM physics    *)
(* this project is aimed at.  Settle it against bernreuther_rmatrix.py  *)
(* before the BSM stage, and record the choice in the repository.       *)
(* ==================================================================== *)

Clear[pHat, kHat, rHat, nHat];

pHat = {0, 0, 1};                          (* beam direction            *)
kHat = {Sin[th], 0, Cos[th]};              (* top direction             *)
rHat = simp[(pHat - Cos[th] kHat)/Sin[th]];   (* -> {-Cos[th],0,Sin[th]} *)
nHat = simp[Cross[rHat, kHat]];               (* -> {0,1,0}              *)

basisCheck = simp[{
   kHat . kHat, rHat . rHat, nHat . nHat,          (* want 1,1,1 *)
   kHat . rHat, kHat . nHat, rHat . nHat,          (* want 0,0,0 *)
   Cross[rHat, kHat] - nHat                  (* want {0,0,0} *)
}];
Print["2. basis {norms, orthogonality, handedness} = ", basisCheck];


(* ==================================================================== *)
(* 3.  INCOMING GLUON POLARIZATION VECTORS                              *)
(*                                                                      *)
(* Gluon 1 moves along +z, gluon 2 along -z.  eps2 is obtained from     *)
(* eps1 by a rotation of pi about the y axis (Jacob-Wick), which flips  *)
(* the x component.  Incoming particles use eps, not eps*.              *)
(* ==================================================================== *)

Clear[eps1, eps2];

eps1[+1] = -(1/Sqrt[2]) {0,  1,  I, 0};
eps1[-1] =  (1/Sqrt[2]) {0,  1, -I, 0};
eps2[+1] = -(1/Sqrt[2]) {0, -1,  I, 0};
eps2[-1] =  (1/Sqrt[2]) {0, -1, -I, 0};

(* transversality, normalisation, and the fact that summing over the    *)
(* two physical helicities reproduces the ghost-free projector used in  *)
(* the old spin-summed calculation.  This is the bridge between Stage 4 *)
(* and Stage 5.                                                          *)
pol1FromHelicities = Sum[Outer[Times, eps1[h], Conjugate[eps1[h]]], {h, {1, -1}}];
pol2FromHelicities = Sum[Outer[Times, eps2[h], Conjugate[eps2[h]]], {h, {1, -1}}];
pol1Physical = -eta + (Outer[Times, k1v, k2v] + Outer[Times, k2v, k1v])/mdot[k1v, k2v];
pol2Physical = -eta + (Outer[Times, k2v, k1v] + Outer[Times, k1v, k2v])/mdot[k2v, k1v];

polCheck = simp[{
   Table[mdot[k1v, eps1[h]], {h, {1, -1}}],          (* want 0,0  *)
   Table[mdot[k2v, eps2[h]], {h, {1, -1}}],          (* want 0,0  *)
   Table[mdot[eps1[h], Conjugate[eps1[h]]], {h, {1, -1}}],   (* want -1,-1 *)
   Flatten[pol1FromHelicities - pol1Physical],       (* want all 0 *)
   Flatten[pol2FromHelicities - pol2Physical]        (* want all 0 *)
}];
Print["3. polarizations {k.e, k.e, e.e*, completeness1, completeness2} = ", polCheck];


(* ==================================================================== *)
(* 4.  DIRAC ALGEBRA  (Dirac representation, metric (+,-,-,-))          *)
(*                                                                      *)
(* The only property of the gamma matrices that matters is the Clifford *)
(* algebra {gam^m, gam^n} = 2 g^mn; everything else follows from it, so *)
(* that is what we check.  slash4[p] = gam^mu p_mu takes CONTRAVARIANT  *)
(* components and the metric supplies the minus signs.                  *)
(* ==================================================================== *)

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

(* one-line definitions on purpose: a line break after a syntactically  *)
(* complete expression silently terminates the input in Mathematica     *)
slash4[p_List] := gam0 p[[1]] - gam1 p[[2]] - gam2 p[[3]] - gam3 p[[4]];
sigDot[v_List] := v[[1]] sig1 + v[[2]] sig2 + v[[3]] sig3;
Sigma4[v_List] := ArrayFlatten[{{sigDot[v], zero2}, {zero2, sigDot[v]}}];
lower4[v_List] := eta . v;

cliffordCheck = simp[Union[Flatten[
   Table[gammaList[[m]] . gammaList[[n]] + gammaList[[n]] . gammaList[[m]]
         - 2 eta[[m, n]] id4, {m, 4}, {n, 4}]]]];
Print["4. Clifford {gam^m,gam^n}=2g^mn residual (want {0}) = ", cliffordCheck];


(* ==================================================================== *)
(* 5.  TOP AND ANTITOP HELICITY SPINORS                                 *)
(*                                                                      *)
(* chiT[h]  : eigenspinor of sigma.kHat  with eigenvalue h              *)
(* chiTb[h] : eigenspinor of sigma.(-kHat) with eigenvalue h            *)
(*            (the antitop moves along -kHat)                           *)
(*                                                                      *)
(* u_h = ( a chi , h b chi ) , v_h = ( h b chibar , a chibar )          *)
(* with a = Sqrt[Et+mt], b = Sqrt[Et-mt].  Note b = |p|/a, which is why *)
(* these solve the Dirac equation.                                      *)
(*                                                                      *)
(* IMPORTANT: h on vSpin is the label carried by chiTb, NOT the physical *)
(* antitop spin.  See section 10 for the correct projection operator.   *)
(* ==================================================================== *)

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

(* ---- the full spinor validation checklist ---- *)
spinorCheck = simp[Union[Flatten[Table[{
   (slash4[ptv]  - mt id4) . uSpin[h],        (* Dirac eqn, want 0     *)
   (slash4[ptbv] + mt id4) . vSpin[h],        (* Dirac eqn, want 0     *)
   Sigma4[kHat]  . uSpin[h] - h uSpin[h],     (* helicity eigenvalue   *)
   Sigma4[-kHat] . vSpin[h] - h vSpin[h]      (* helicity eigenvalue   *)
}, {h, {1, -1}}]]]];
normCheck = simp[Table[
   {uBar[h] . uSpin[h], vBar[h] . vSpin[h], Conjugate[uSpin[h]] . uSpin[h]},
   {h, {1, -1}}]];
Print["5a. spinor residuals (want {0}) = ", spinorCheck];
Print["5b. {ubar u, vbar v, u^dag u} (want {2mt, -2mt, 2Et}) = ", normCheck];

(* --------------------------------------------------------------------- *)
(* 5c.  NUMERICAL spinor check                                            *)
(*                                                                        *)
(* With bet as the free symbol the residuals in 5a are proportional to    *)
(*   aSpin (1 - mt) - bet bSpin ,                                         *)
(* which vanishes because 1 - mt = bSpin^2 and bet = aSpin bSpin.  But    *)
(* FullSimplify cannot prove it through the nested Sqrt[1 +- Sqrt[1-b^2]],*)
(* so 5a prints long expressions that are secretly zero.  Check them      *)
(* numerically instead - this is conclusive.                              *)
(* --------------------------------------------------------------------- *)

spinorCheckNum = Chop[Max[Abs[Flatten[
   N[spinorCheck /. {bet -> 0.8693, th -> 1.1}]]]]];
Print["5c. spinor residuals, NUMERICAL (want 0) = ", spinorCheckNum];



(* ==================================================================== *)
(* 6.  FIXED-HELICITY AMPLITUDES                                        *)
(*                                                                      *)
(* Colour is stripped off and handled separately in a 3-dim basis       *)
(*     c1 = T^a T^b  (t-channel), c2 = T^b T^a  (u-channel),            *)
(*     c3 = f^{abc} T^c  (s-channel)                                    *)
(* so At, Au, As below are the COLOUR-STRIPPED coefficients.            *)
(*                                                                      *)
(* Internal momenta follow the fermion arrow:  qT = pt - k1 satisfies   *)
(* qT + k1 = pt at the top-side vertex.  Same for qU = pt - k2.         *)
(*                                                                      *)
(* Every definition is memoized (f[x_] := f[x] = ...).  Without this,   *)
(* building R calls Avec 128 times and re-runs FullSimplify each time.  *)
(* ==================================================================== *)

Clear[qTv, qUv, denT, denU, propT, propU,
      epsSlash1, epsSlash2, At, Au, As, Avec,
      p3v, V3Tensor, v3Contract, qCurrent];

qTv = ptv - k1v;
qUv = ptv - k2v;

denT = simp[mdot[qTv, qTv] - mt^2];    (* = t - mt^2 = -(s/2)(1 - b Cos) *)
denU = simp[mdot[qUv, qUv] - mt^2];    (* = u - mt^2 = -(s/2)(1 + b Cos) *)
Print["6a. {denT, denU} = ", {denT, denU}];

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

(* contract the two external gluon indices; lower4 supplies eps_mu so   *)
(* the result X^sigma is contravariant                                  *)
v3Contract[h1_, h2_] := v3Contract[h1, h2] = simp[Table[
   Sum[lower4[eps1[h1]][[m]] lower4[eps2[h2]][[n]] V3Tensor[m, n, sg],
       {m, 4}, {n, 4}], {sg, 4}]];

(* the fixed-helicity top current J^rho = ubar gam^rho v.               *)
(* simp[] here matters: without it Conjugate[] wrappers propagate into  *)
(* all 16 s-channel amplitudes.                                         *)
qCurrent[ht_, hb_] := qCurrent[ht, hb] = simp[
   Table[uBar[ht] . gammaList[[r]] . vSpin[hb], {r, 4}]];

(* quark vertex (I gs) times gluon propagator (-I g_rs / s)             *)
As[h1_, h2_, ht_, hb_] := As[h1, h2, ht, hb] = simp[
   (I gs) (-I)/sMand mdot[qCurrent[ht, hb], v3Contract[h1, h2]]];

Avec[h1_, h2_, ht_, hb_] := Avec[h1, h2, ht, hb] =
   {At[h1, h2, ht, hb], Au[h1, h2, ht, hb], As[h1, h2, ht, hb]};

helicities  = {+1, -1};
finalStates = {{+1, +1}, {+1, -1}, {-1, +1}, {-1, -1}};

(* check the three-gluon contraction against the closed form that       *)
(* follows from eps1.k2 = eps2.k1 = 0 in this frame                     *)
v3Check = simp[Union[Flatten[Table[
   v3Contract[h1, h2] - gs mdot[eps1[h1], eps2[h2]] (k1v - k2v),
   {h1, helicities}, {h2, helicities}]]]];
Print["6b. three-gluon contraction residual (want {0}) = ", v3Check];

(* build and cache all 16 amplitude vectors ONCE, up front *)
Print["6c. caching 16 amplitudes, timing: ",
   AbsoluteTiming[Do[Avec[h1, h2, ht, hb],
      {h1, helicities}, {h2, helicities},
      {ht, helicities}, {hb, helicities}];][[1]], " s"];
Print["6d. Avec[+,+,+,+] = ", Avec[1, 1, 1, 1]];


(* ==================================================================== *)
(* 7.  COLOUR                                                           *)
(*                                                                      *)
(* Colour is kept completely separate from the Dirac algebra.  Putting  *)
(* SU(3) generators inside a Dirac chain means two different            *)
(* non-commutative structures in one product, which is what broke the   *)
(* original FeynCalc attempt ("failed to isolate the colour algebraic   *)
(* expression").  Here c1, c2, c3 are ordinary COMMUTING symbols, so    *)
(* each amplitude is just a 3-vector and the Dirac algebra never sees   *)
(* them.  The three diagrams are never added together: the bilinear     *)
(* form A_i^T . C . A_j^* automatically produces all nine diagram-pair  *)
(* interference terms with the right colour weights.                    *)
(*                                                                      *)
(* C_mn = Sum_colours c_m c_n^dagger with the quark indices traced:     *)
(*   C11 = C22 = 16/3 , C12 = C21 = -2/3 , C33 = 12 ,                   *)
(*   C13 = -C31 = 6i , C23 = -C32 = -6i                                 *)
(* using Tr(T^a T^b) = d^ab/2, T^a T^b T^a = -T^b/6, f^abc f^abc = 24,  *)
(* Tr(T^a T^b T^c) = (d^abc + i f^abc)/4.                               *)
(*                                                                      *)
(* Then  Sum_colours M_i M_j^*  =  A_i^T . C . A_j^*                    *)
(* ==================================================================== *)

Clear[Ccolor, colorProduct];

Ccolor = {
   {16/3, -2/3,  6 I},
   {-2/3, 16/3, -6 I},
   {-6 I,  6 I,   12}
};

(* NOTE the _List patterns.  If Avec is ever undefined, its head is     *)
(* Avec rather than List, this pattern silently fails to match, and the *)
(* whole calculation returns unevaluated symbols with NO error message. *)
(* The guard below turns that silent failure into a loud one.           *)
colorProduct[a_List, b_List] := a . Ccolor . Conjugate[b];
colorProduct::badarg = "colorProduct received non-list arguments: `1`, `2`. \
Definitions are missing - the kernel has probably been restarted.";
colorProduct[a_, b_] := (Message[colorProduct::badarg, a, b]; Abort[]);

Print["7a. Ccolor Hermitian: ", ConjugateTranspose[Ccolor] === Ccolor];
Print["7b. colorProduct[++;++, ++;++] (expect 0, c1-c2 = i c3 cancellation) = ",
      simp[colorProduct[Avec[1, 1, 1, 1], Avec[1, 1, 1, 1]]]];


(* ==================================================================== *)
(* 7b.  WARD IDENTITY  --  the strongest single check in the notebook   *)
(*                                                                      *)
(* Gauge invariance says: replace an external gluon's polarisation      *)
(* vector by its own momentum, and the amplitude must vanish.           *)
(*                                                                      *)
(*      eps1 -> k1     =>     k1_mu M^{mu nu} = 0                       *)
(*                                                                      *)
(* WHY THIS IS WORTH DOING: it only holds for the SUM of the three      *)
(* diagrams, so one number simultaneously tests                         *)
(*    - the relative signs of the t, u and s channels                   *)
(*    - the three-gluon vertex convention (signs and momentum routing)  *)
(*    - the internal propagator momenta q_T = pt - k1 etc.              *)
(*    - the external polarisation treatment                             *)
(* Nothing else in the notebook probes those four things together.      *)
(*                                                                      *)
(* NOTE ON THE FORM: because c1 - c2 = i c3 (the SU(3) commutator), the *)
(* three colour-stripped pieces do NOT each vanish separately.  The     *)
(* gauge-invariant combination is  i (At - Au) + As = 0.                *)
(* ==================================================================== *)

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

Print["7c. Ward identity  i(At-Au)+As  (want {0}) = ", wardCheck];


(* ==================================================================== *)
(* --- SANITY GATE: run this before anything expensive ---              *)
(* If any line below is wrong, STOP.  The kernel has lost its state and *)
(* everything downstream will be meaningless unevaluated symbols.       *)
(* ==================================================================== *)

Print["GATE simp defined       (want True) : ", Head[simp[1]] =!= simp];
Print["GATE Avec head          (want List) : ", Head[Avec[1, 1, 1, 1]]];
Print["GATE helicities         (want {1,-1}): ", helicities];
Print["GATE colorProduct works (want 0)    : ",
      simp[colorProduct[Avec[1, 1, 1, 1], Avec[1, 1, 1, 1]]]];

(* Save the state so a kernel restart costs one second, not one hour.   *)
DumpSave["stage5.mx", {simp, simpN, assumptions, mt, bet, eta, mdot,
   k1v, k2v, ptv, ptbv, kHat, rHat, nHat, eps1, eps2,
   gam0, gam1, gam2, gam3, gammaList, id4, sig1, sig2, sig3,
   slash4, sigDot, Sigma4, lower4, chiT, chiTb, uSpin, vSpin, uBar, vBar,
   At, Au, As, Avec, qCurrent, v3Contract, Ccolor, colorProduct,
   helicities, finalStates}];
Print["state saved to stage5.mx  (restore later with  << \"stage5.mx\" )"];


(* ==================================================================== *)
(* 8.  THE TRACE  -- the decisive test, and the cheap one               *)
(*                                                                      *)
(*   Tr R = (1/256) Sum_{h1,h2,ht,htb} A^T C A^*                        *)
(*                                                                      *)
(* By construction this is exactly the old spin-summed, initial-state   *)
(* averaged |M|^2.  Prove this BEFORE building any off-diagonal entry:  *)
(* it needs only the 4 diagonal sums, not the full 4x4.                 *)
(*                                                                      *)
(* 1/4  = two physical helicities for each of two incoming gluons       *)
(* 1/64 = 8 x 8 incoming adjoint colours                                *)
(* ==================================================================== *)

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

Print["8a. Tr R                  = ", traceR];
Print["8b. Tr R - |M|^2_standard (want 0) = ", simp[traceR - M2ref]];
Print["8c. Tr R - 4 gs^4 Atilde  (want 0) = ", simp[traceR - 4 gs^4 Atilde]];


(* ==================================================================== *)
(* 9.  THE 4x4 DENSITY MATRIX, ELEMENT BY ELEMENT                       *)
(*                                                                      *)
(* R is Hermitian, so only the 10 upper-triangle entries are            *)
(* independent; at LO it is also real, and Bose symmetry of the two     *)
(* gluons plus parity cuts it further to just FIVE independent numbers. *)
(* Computing one element at a time keeps each FullSimplify small - far  *)
(* better than one big Table followed by a global simplification.       *)
(*                                                                      *)
(*   R_ij = (1/256) Sum_{h1,h2} A_i^T C A_j^*                           *)
(*                                                                      *)
(* Basis order  |++>, |+->, |-+>, |-->  with the first label the top.   *)
(* Initial helicities and colours are summed and averaged; the final    *)
(* top spins are kept OPEN.  That is what turns a scalar rate into a    *)
(* 4x4 spin matrix.                                                     *)
(*                                                                      *)
(* Hermiticity means only the 10 upper-triangle entries are independent.*)
(* Computing them one at a time keeps each FullSimplify small.          *)
(* ==================================================================== *)

Clear[Rel, Rhel];

Rel[i_, j_] := Rel[i, j] = simp[(1/(4*64)) Sum[
   colorProduct[Avec[h1, h2, Sequence @@ finalStates[[i]]],
                Avec[h1, h2, Sequence @@ finalStates[[j]]]],
   {h1, helicities}, {h2, helicities}]];

Do[Print["9a. R[", i, ",", j, "] = ", Rel[i, j]], {i, 4}, {j, i, 4}];

Rhel = Table[If[j >= i, Rel[i, j], Conjugate[Rel[j, i]]], {i, 4}, {j, 4}];

Print["9b. Hermiticity residual (want {0}) = ",
      simp[Union[Flatten[Rhel - ConjugateTranspose[Rhel]]]]];
Print["9c. trace consistency    (want 0)   = ", simp[Tr[Rhel] - traceR]];

(* --- closed form to compare against, derived from the Appendix-B      *)
(* --- coefficients and verified numerically at four (bet,th) points.   *)
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
Print["9d. Rhel - RTarget (want {0}) = ",
      simp[Union[Flatten[Rhel - RTarget]]]];


(* ==================================================================== *)
(* 10.  BLOCH-FANO PROJECTION                                           *)
(*                                                                      *)
(*   R = Atil (1 x 1) + B+_i (sig_i x 1) + B-_i (1 x sig_i)             *)
(*                     + C_ij (sig_i x sig_j)                           *)
(*                                                                      *)
(* Using Tr sig_i = 0 and Tr(sig_i sig_j) = 2 delta_ij:                 *)
(*   Atil = Tr R /4 ,  C_ij = Tr[R (Sig_i x Sig_j)]/4 , etc.            *)
(*                                                                      *)
(* THE OPERATORS ARE NOT sig1, sig2, sig3.  Two separate effects:       *)
(*                                                                      *)
(*  (a) The spinor basis {chiT[+1], chiT[-1]} is the helicity basis,    *)
(*      so k is the quantisation axis.  In that basis                   *)
(*         sigma.kHat -> sig3 ,  sigma.nHat -> sig2 ,                   *)
(*         sigma.rHat -> -sig1                                          *)
(*      The minus sign on r is because rHat maps to -x when kHat is     *)
(*      rotated onto z.  Assuming sig_r = sig1 flips C_rk and C_kr.     *)
(*                                                                      *)
(*  (b) For the antitop the label carried by vSpin is related to the    *)
(*      physical spin by charge conjugation, so the correct operator is *)
(*      -(sigma.a)^* , not sigma.a and not -(sigma.a).                  *)
(*      Evaluating -(sigma.a)^* in the chiTb basis happens to give the  *)
(*      SAME three matrices as for the top, which is why both slots     *)
(*      below use SigR, SigN, SigK.                                     *)
(*                                                                      *)
(* Verified numerically to machine precision at four (bet,th) points.   *)
(* Getting either effect wrong still leaves Tr R correct, so the trace  *)
(* test alone does NOT validate the off-diagonal coefficients.          *)
(* ==================================================================== *)

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

Print["10a. Atil = ", Atil];
Print["10b. B+ (want {0,0,0}) = ", Bplus];
Print["10c. B- (want {0,0,0}) = ", Bminus];
Do[Print["10d. C_", axisName[[i]], axisName[[j]], " = ", Cij[[i, j]]],
   {i, 3}, {j, 3}];


(* ==================================================================== *)
(* 11.  MATCH THE APPENDIX-B TARGETS                                    *)
(*                                                                      *)
(* Five independent structures at LO in the SM:                         *)
(*   Atil, C_rr, C_nn, C_kk, C_rk = C_kr                                *)
(* plus the zeros B+_i = B-_i = C_rn = C_nr = C_nk = C_kn = 0.          *)
(* ==================================================================== *)

Clear[CrrT, CnnT, CkkT, CrkT];

CrrT = -Fg (1 - bet^2 (2 - bet^2)(1 + Sin[th]^4));
CnnT = -Fg (1 - 2 bet^2 + bet^4 (1 + Sin[th]^4));
CkkT = -Fg (1 - (bet^2/2) Sin[2 th]^2 - bet^4 (1 + Sin[th]^4));
CrkT =  Fg Sqrt[1 - bet^2] bet^2 Sin[2 th] Sin[th]^2;

Print["11a. Atil - gs^4 Atilde (want 0) = ", simp[Atil - gs^4 Atilde]];
Print["11b. Crr  - gs^4 CrrT   (want 0) = ", simp[Cij[[1, 1]] - gs^4 CrrT]];
Print["11c. Cnn  - gs^4 CnnT   (want 0) = ", simp[Cij[[2, 2]] - gs^4 CnnT]];
Print["11d. Ckk  - gs^4 CkkT   (want 0) = ", simp[Cij[[3, 3]] - gs^4 CkkT]];
Print["11e. Crk  - gs^4 CrkT   (want 0) = ", simp[Cij[[1, 3]] - gs^4 CrkT]];
Print["11f. Ckr  - gs^4 CrkT   (want 0) = ", simp[Cij[[3, 1]] - gs^4 CrkT]];
Print["11g. off-plane zeros (want {0}) = ",
      simp[Union[{Cij[[1,2]], Cij[[2,1]], Cij[[2,3]], Cij[[3,2]]}]]];


(* ==================================================================== *)
(* 12.  NUMERICAL SPOT CHECKS                                           *)
(*                                                                      *)
(* Purpose: separate a PHYSICS error from a FullSimplify timeout.  If   *)
(* the numbers here are right but the symbolic check above is not zero, *)
(* the algebra is fine and the simplifier just gave up.                 *)
(*                                                                      *)
(* Run these FIRST whenever something is slow or suspicious.  They      *)
(* separate a genuine physics error from a FullSimplify timeout.        *)
(* Reference values (gs = 1), computed independently:                   *)
(*   bet = 0.86930, th = 1.1  ->  Tr R = 0.3113740409                   *)
(*   bet = 0.63938, th = 2.3  ->  Tr R = 0.3317768301                   *)
(* ==================================================================== *)

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

Print["12a. ", numCheck[0.86930, 1.1]];
Print["12b. ", numCheck[0.63938, 2.3]];
Print["12c. ", numCheck[0.40000, 0.7]];


(* ==================================================================== *)
(* 13.  RESTORING PHYSICAL UNITS                                        *)
(*                                                                      *)
(* Everything above used Et = 1.  R and all Bloch-Fano coefficients are *)
(* dimensionless, so nothing needs rescaling.  To express results in    *)
(* s and mt instead of bet, substitute                                  *)
(*     bet -> Sqrt[1 - 4 mtPhys^2/sPhys]                                *)
(* and, if you want Mandelstam variables,                               *)
(*     t -> mtPhys^2 - sPhys/2 + sPhys bet Cos[th]/2                    *)
(*     u -> mtPhys^2 - sPhys/2 - sPhys bet Cos[th]/2                    *)
(* ==================================================================== *)

toPhysical = {bet -> Sqrt[1 - 4 mtPhys^2/sPhys]};
Print["13. Atil in (s, mt): ", Simplify[Atil /. toPhysical,
      {sPhys > 4 mtPhys^2, mtPhys > 0, 0 < th < Pi, gs > 0}]];


(* ==================================================================== *)
(*                     ---  NLO STARTS HERE  ---                        *)
(*                                                                      *)
(* At one loop the amplitude picks up a correction M(1), and the cross  *)
(* section gets the INTERFERENCE 2 Re( M(1) M(0)* ) - not the square of *)
(* the loop.  So the tree amplitudes above are still needed.            *)
(*                                                                      *)
(* Each loop diagram has loop momenta in the numerator, giving TENSOR   *)
(* integrals.  Passarino-Veltman reduction rewrites all of them in      *)
(* terms of a small fixed basis of SCALAR integrals with 1, 2, 3 or 4   *)
(* propagators: A0, B0, C0, D0.  That is why "list the one-loop         *)
(* integrals" is a finite task even though there are ~30-40 diagrams.   *)
(*                                                                      *)
(* THE PHYSICS TO WATCH FOR: above threshold the loop integrals develop *)
(* IMAGINARY parts (Im B0 = Pi Sqrt[lambda]/s).  Those give relative    *)
(* phases between amplitudes, which switch on the T-odd coefficients    *)
(* B+_n, B-_n, C_rn, C_nk that are identically zero at tree level.  If  *)
(* they come back zero at one loop, something is wrong.                 *)
(* ==================================================================== *)


(* ==================================================================== *)
(* 14.  CATALOGUE OF ONE-LOOP SCALAR INTEGRALS                          *)
(*                                                                      *)
(*  External legs: 2 massless (gluons, p^2=0), 2 massive (tops, p^2=mt^2) *)
(*  Internal lines: gluons (massless) and tops (mass mt)                *)
(*  2 -> 2 process, so nothing beyond boxes.                            *)
(*                                                                      *)
(*  NOTE: this hand list is for UNDERSTANDING.  For production, generate *)
(*  the list automatically (FeynArts+FormCalc, or Package-X) and use    *)
(*  this to sanity-check the output.  The exact set depends on the      *)
(*  diagram set and gauge, and some entries are related by crossing.    *)
(* ==================================================================== *)


Clear[lpCatalogue];

lpCatalogue = <|

 "A0" -> {
   {"A0[mt^2]",        "massive tadpole, from top self-energy"},
   {"A0[0] = 0",       "scaleless, vanishes in dim reg"}
 },

 "B0" -> {
   {"B0[s, 0, 0]",        "gluon self-energy insertion, massless loop"},
   {"B0[s, mt^2, mt^2]",  "gluon self-energy, top loop"},
   {"B0[t, 0, mt^2]",     "t-channel bubble"},
   {"B0[u, 0, mt^2]",     "u-channel bubble"},
   {"B0[mt^2, 0, mt^2]",  "top self-energy, on shell -> needs care, IR"},
   {"B0[0, mt^2, mt^2]",  "zero-momentum bubble"},
   {"B0[0, 0, 0] = 0",    "scaleless, vanishes"}
 },

 "C0" -> {
   {"C0[0,0,s; 0,0,0]",           "three-gluon vertex, massless loop"},
   {"C0[0,0,s; mt,mt,mt]",        "gg -> g* through a top loop"},
   {"C0[mt^2,mt^2,s; mt,0,mt]",   "g t tbar vertex correction"},
   {"C0[mt^2,mt^2,s; 0,mt,0]",    "g t tbar vertex, gluon exchange"},
   {"C0[0,mt^2,t; 0,0,mt]",       "mixed triangle, t-channel"},
   {"C0[0,mt^2,u; 0,0,mt]",       "mixed triangle, u-channel"},
   {"(+ mass permutations of the above)", ""}
 },

 "D0" -> {
   {"D0[0,0,mt^2,mt^2; s,t; 0,0,0,mt]",    "box, mostly massless internal"},
   {"D0[0,0,mt^2,mt^2; s,t; 0,0,mt,mt]",   "box"},
   {"D0[0,0,mt^2,mt^2; s,t; 0,mt,mt,mt]",  "box"},
   {"(+ u-channel partners of each)", ""}
 }
|>;

lpShowCatalogue[] := Module[{rows},
  rows = Join @@ Table[
     Join[{{Style[k, Bold], ""}}, lpCatalogue[k]],
     {k, Keys[lpCatalogue]}];
  Grid[rows, Alignment -> Left, Frame -> All,
       Background -> {None, {{None, LightYellow}}}]
];

Print["14. one-loop scalar integral catalogue:"];
lpShowCatalogue[]


(* ==================================================================== *)
(* 15.  ANALYTIC B0                                                     *)
(*                                                                      *)
(*  Derivation:                                                          *)
(*    B0 = Delta - Integrate[ Log[(p2 x^2 - (p2+m02-m12) x + m02 - I eps)/mu2], *)
(*                            {x,0,1} ]                                  *)
(*  Factor the quadratic as p2 (x-xp)(x-xm) with                         *)
(*    xpm = ( (p2+m02-m12) +- Sqrt[lambda] ) / (2 p2),                   *)
(*    lambda = (p2+m02-m12)^2 - 4 p2 m02       (Kallen function)         *)
(*  and use  Integrate[Log[x-a],{x,0,1}] = (1-a)Log[1-a] + a Log[-a] - 1 *)
(*  to get                                                               *)
(*    B0fin = 2 - Log[p2/mu2]                                            *)
(*            - Sum_i [ (1-x_i) Log[1-x_i] + x_i Log[-x_i] ]             *)
(*                                                                       *)
(*  The I eps matters.  We put -I eps on BOTH the masses and on p2.      *)
(*  Getting this wrong gives the correct magnitude with a wrong sign on  *)
(*  Im B0, or an error of 2 Pi I in the spacelike region.                *)
(* ==================================================================== *)


Clear[lpEps, B0fin, B0finSeries];

lpEps = 10^-20;   (* infinitesimal; small enough not to perturb the answer *)

B0fin[p2_?NumericQ, m02_?NumericQ, m12_?NumericQ, mu2_?NumericQ] :=
 Module[{M0, M1, bb, sq, xp, xm, res},

   M0 = m02 - I lpEps;
   M1 = m12 - I lpEps;

   (* --- degenerate case p^2 = 0 --- *)
   If[Abs[p2] < 10^-10,
      Return[
        If[Abs[m02 - m12] < 10^-10,
           -Log[M0/mu2],                                   (* equal masses *)
           1 - (M0 Log[M0/mu2] - M1 Log[M1/mu2])/(M0 - M1)  (* general      *)
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
(* 15b.  Closed-form special cases -- useful as independent checks       *)
(* -------------------------------------------------------------------- *)

Clear[B0massless, B0oneMass, B0zeroMom];

(* B0(p^2, 0, 0) = 2 - Log[p2/mu2] + I Pi   for p2 > 0                *)
(* Note the argument is -(p2) - I eps, NOT p2 - I eps: for timelike p2 *)
(* the branch point is approached from the negative side, and writing  *)
(* it the other way loses the I Pi entirely.                           *)
B0massless[p2_, mu2_] := 2 - Log[(-p2 - I lpEps)/mu2];

(* B0(p^2, 0, m^2) *)
(* threshold here is just m^2 (one line is massless), so this already  *)
(* has an imaginary part for p2 > m^2.  The eps must be NEGATIVE: with  *)
(* +I eps you land on the wrong side and Im comes out with a flipped    *)
(* sign, i.e. off by 2 x the right answer.                              *)
B0oneMass[p2_, m2_, mu2_] :=
   2 - Log[m2/mu2] + (m2/p2 - 1) Log[1 - p2/m2 - I lpEps];

(* B0(0, m0^2, m1^2) *)
B0zeroMom[m02_, m12_, mu2_] :=
   If[Abs[m02 - m12] < 10^-10, -Log[m02/mu2],
      1 - (m02 Log[m02/mu2] - m12 Log[m12/mu2])/(m02 - m12)];


(* --- 15b-extra.  Two more closed forms, both NEEDED for gg -> t tbar - *)

Clear[B0equalMass, B0onShell];

(* B0(p^2, m^2, m^2)  -- this is the top loop in the gluon self-energy,  *)
(* i.e. B0(s; mt, mt) in the catalogue.  beta here is the velocity of    *)
(* the pair inside the loop, NOT the top velocity of the outer process.  *)
B0equalMass[p2_, m2_, mu2_] := Module[{be},
   be = Sqrt[1 - 4 m2/(p2 - I lpEps)];
   2 - Log[m2/mu2] + be Log[(be - 1)/(be + 1)]];

(* B0(m^2, 0, m^2) -- the on-shell top self-energy.  The whole           *)
(* momentum-dependent piece vanishes at p2 = m^2, leaving just a log.    *)
B0onShell[m2_, mu2_] := 2 - Log[m2/mu2];

Print["  equal-mass form vs general (want 0):  ",
  Chop[B0equalMass[10.^6, mtval^2, muSq]
     - B0fin[10.^6, mtval^2, mtval^2, muSq]]];
Print["  on-shell  form vs general (want 0):  ",
  Chop[B0onShell[mtval^2, muSq]
     - B0fin[mtval^2, 10.^-20, mtval^2, muSq]]];


(* -------------------------------------------------------------------- *)
(* 15c.  Numerical ground truth: direct Feynman-parameter integration    *)
(*       Use this to validate ANY loop-integral implementation,          *)
(*       including Collier once it is built.                             *)
(* -------------------------------------------------------------------- *)

Clear[B0numeric];

B0numeric[p2_?NumericQ, m02_?NumericQ, m12_?NumericQ, mu2_?NumericQ,
          eps_: 10^-9] :=
  NIntegrate[
     -Log[(p2 x^2 - (p2 + m02 - m12) x + m02 - I eps)/mu2],
     {x, 0, 1},
     WorkingPrecision -> 20, AccuracyGoal -> 12, MaxRecursion -> 30];


(* -------------------------------------------------------------------- *)
(* 15d.  VALIDATION                                                      *)
(*       Twelve points spanning spacelike, below threshold, on threshold *)
(*       and above threshold.  Reference values computed independently.  *)
(* -------------------------------------------------------------------- *)

Clear[mHval, mtval, muSq, lpTests, lpCheck];

mHval = 125.25;
mtval = 172.5;
muSq  = mtval^2;          (* <-- STATE THIS WHENEVER YOU QUOTE A NUMBER *)

lpTests = {
 {"B0[62500,  mH^2, mt^2]  <-- TASK 1",  62500.,  mHval^2, mtval^2},
 {"B0[250,    mH^2, mt^2]  <-- alt read",  250.,  mHval^2, mtval^2},
 {"B0[62500,  mt^2, mH^2]  symmetry",    62500.,  mtval^2, mHval^2},
 {"B0[100000, mH^2, mt^2]  above thr",  100000.,  mHval^2, mtval^2},
 {"B0[0,      mH^2, mt^2]",                  0.,  mHval^2, mtval^2},
 {"B0[mt^2,   0,    mt^2]",           mtval^2, 10.^-20, mtval^2},
 {"B0[0,      mt^2, mt^2]",                  0.,  mtval^2, mtval^2},
 {"B0[-50000, 0,    mt^2]  spacelike",  -50000., 10.^-20, mtval^2},
 {"B0[-10^6,  mt^2, mt^2]  spacelike",   -10.^6,  mtval^2, mtval^2},
 {"B0[10^6,   mt^2, mt^2]  s=1TeV^2",     10.^6,  mtval^2, mtval^2},
 {"B0[10^6,   0,    0   ]  massless",     10.^6, 10.^-20, 10.^-20},
 {"B0[29929,  0,    mt^2]  p2=mt^2",     29929., 10.^-20, mtval^2}
};

lpCheck[] := Module[{a, b, bad = 0},
  Print["15d. B0 validation, mu^2 = mt^2 = ", muSq];
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

(* Expected output (mu^2 = 172.5^2):
     B0[62500,  mH^2, mt^2]  = +0.992443465 + 0.000000000 i     <-- TASK 1
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
(* 15e.  Structural checks that do not rely on the reference numbers     *)
(* -------------------------------------------------------------------- *)

Print["\n15e. structural checks"];

Print["  symmetry B0(p2,m0,m1) = B0(p2,m1,m0): ",
  Chop[B0fin[62500., mHval^2, mtval^2, muSq]
     - B0fin[62500., mtval^2, mHval^2, muSq]]];

Print["  special case B0(p2,0,0)   vs general: ",
  Chop[B0massless[10.^6, muSq] - B0fin[10.^6, 10.^-20, 10.^-20, muSq]]];

Print["  special case B0(p2,0,m^2) vs general: ",
  Chop[B0oneMass[62500., mtval^2, muSq] - B0fin[62500., 10.^-20, mtval^2, muSq]]];

Print["  special case B0(0,m0,m1)  vs general: ",
  Chop[B0zeroMom[mHval^2, mtval^2, muSq] - B0fin[0., mHval^2, mtval^2, muSq]]];

Print["  mu-independence of dB0/dLog[mu^2] (want +1): ",
  Chop[(B0fin[62500., mHval^2, mtval^2, 2 muSq]
      - B0fin[62500., mHval^2, mtval^2, muSq])/Log[2]]];

Print["  Im B0 = 0 below threshold (mH+mt)^2 = ", (mHval + mtval)^2, " : ",
  Chop[Im[B0fin[62500., mHval^2, mtval^2, muSq]]]];

Print["  Im B0 > 0 above threshold: ",
  Im[B0fin[100000., mHval^2, mtval^2, muSq]] > 0];

(* Analytic check of the absorptive part:
   Im B0(s,m0^2,m1^2) = Pi Sqrt[lambda]/s   for s > (m0+m1)^2          *)
Print["  Im B0 vs Pi Sqrt[lambda]/s: ",
 Module[{s = 100000., a = mHval^2, b = mtval^2, lam},
   lam = (s - a - b)^2 - 4 a b;
   Chop[Im[B0fin[s, a, b, muSq]] - Pi Sqrt[lam]/s]]];


(* -------------------------------------------------------------------- *)
(* 15f.  Compiled version, for phase-space scans                         *)
(*       (Marcel's protip: compile the inner loop)                       *)
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

Print["\n15f. compiled vs interpreted: ",
  Chop[B0compiled[62500., mHval^2, mtval^2, muSq]
     - B0fin[62500., mHval^2, mtval^2, muSq]]];

Print["     timing, 10^4 evaluations:"];
Print["       interpreted: ",
  First@AbsoluteTiming[Do[B0fin[62500. + n, mHval^2, mtval^2, muSq], {n, 10^4}]], " s"];
Print["       compiled:    ",
  First@AbsoluteTiming[Do[B0compiled[62500. + n, mHval^2, mtval^2, muSq], {n, 10^4}]], " s"];


(* ==================================================================== *)
(* 16.  PHASE-SPACE POINT EVALUATION                                    *)
(* ==================================================================== *)

Clear[psWhichFree, psExtract, psEvaluate, psChecks, psEntanglement,
      psReport, psSelfTest, numAssert, momnum, result];


(* -------------------------------------------------------------------- *)
(* 16a.  Diagnostic: which of mt / bet is the free symbol?               *)
(*       Run this first if anything looks wrong.                         *)
(* -------------------------------------------------------------------- *)

psWhichFree[] := Module[{mtFree, betFree},
  mtFree  = (Head[mt]  === Symbol) && (mt  === HoldForm@mt  || ! NumericQ[mt]);
  betFree = (Head[bet] === Symbol) && (bet === HoldForm@bet || ! NumericQ[bet]);
  Print["16a. symbol status"];
  Print["   mt  evaluates to : ", mt];
  Print["   bet evaluates to : ", bet];
  Print["   (exactly one should be a bare symbol; the other a formula)"];
];


(* -------------------------------------------------------------------- *)
(* 16b.  Extract invariants and spin basis from four momenta             *)
(*       Input: {k1, k2, pt, ptbar}, each {E, px, py, pz}                *)
(* -------------------------------------------------------------------- *)

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
(* 16c.  Evaluate.  gs = 1, so everything is in units of gs^4.           *)
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
(* 16d.  Checks, including the azimuth-invariance argument               *)
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
(* 16e.  Entanglement.  Concurrence and tr C are BASIS INDEPENDENT, so   *)
(*       they can be quoted before the convention questions are settled. *)
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
(* 16f.  Formatted report - this is what you send                        *)
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
  Print["   phi (azimuth)  = ", f[kin["phi"] 180/Pi],
        " deg   <-- top NOT in the x-z plane"];
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
  Print["   B+   = ", v[ev["Bplus"]],  "   zero at LO (parity, no absorptive part)"];
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
(* 16g.  Self-test against independently computed reference values       *)
(* -------------------------------------------------------------------- *)

psSelfTest[res_] := Module[{ev, kin, ref, got, dev},
  If[res === $Failed, Print["16g. self-test skipped (evaluation failed)"];
     Return[$Failed]];
  kin = res["kin"]; ev = res["ev"];
  ref = <|"beta" -> 0.9382345123, "cosTheta" -> -0.1667602549,
          "Atil" -> 0.0475904733, "trR" -> 0.1903618931,
          "Crr" -> 0.0362479415, "Cnn" -> -0.0295085436,
          "Ckk" -> 0.0219278671, "Crk" -> -0.0038479679|>;
  got = <|"beta" -> kin["beta"], "cosTheta" -> kin["cosTheta"],
          "Atil" -> ev["Atil"], "trR" -> ev["trR"],
          "Crr" -> ev["C"][[1,1]], "Cnn" -> ev["C"][[2,2]],
          "Ckk" -> ev["C"][[3,3]], "Crk" -> ev["C"][[1,3]]|>;
  Print["16g. self-test vs independent reference (tolerance 1e-8)"];
  Do[
    dev = Abs[got[k] - ref[k]];
    Print["   ", StringPadRight[k, 10],
          " got ", NumberForm[N[got[k]], {12,10}],
          "   ref ", NumberForm[N[ref[k]], {12,10}],
          "   ", If[dev < 10^-8, "OK", "** DEVIATION " <> ToString[dev]]],
    {k, Keys[ref]}];
];


(* -------------------------------------------------------------------- *)
(* 16h.  RUN                                                             *)
(* -------------------------------------------------------------------- *)

psWhichFree[];
Print[""];

momnum = {{500.`, 0.`, 0.`, 500.`},
          {500.`, 0.`, 0.`, -500.`},
          {500.`,  461.604563622633`,  29.53432289997799`, -78.23011321478687`},
          {500.`, -461.604563622633`, -29.53432289997799`,  78.23011321478687`}};

result = psReport[momnum];
Print[""];
psSelfTest[result];


(* ==================================================================== *)
(*  EXPECTED VALUES                                                     *)
(*                                                                      *)
(*   sqrt(s) = 1000.000000   m_t = 173.0000000   beta = 0.9382345123    *)
(*   cos(theta) = -0.1667602549   theta = 1.7383393226                  *)
(*   phi = 3.660901 deg          mt (rescaled) = 0.3460000000           *)
(*                                                                      *)
(*   kHat = { 0.9839854697,  0.0629572298, -0.1667602549}               *)
(*   rHat = { 0.1664199678,  0.0106478606,  0.9859974733}               *)
(*   nHat = {-0.0638513095,  0.9979594232,  0.0000000000}               *)
(*                                                                      *)
(*   R / gs^4 =                                                         *)
(*     0.0695183404   0.0038479679   0.0038479679   0.0657564851        *)
(*     0.0038479679   0.0256626062   0.0067393980  -0.0038479679        *)
(*     0.0038479679   0.0067393980   0.0256626062  -0.0038479679        *)
(*     0.0657564851  -0.0038479679  -0.0038479679   0.0695183404        *)
(*                                                                      *)
(*   Tr R = 0.1903618931        Atil = 0.0475904733                     *)
(*   B+ = B- = {0, 0, 0}                                                *)
(*                                                                      *)
(*   C_rr = +0.0362479415    /Atil = +0.76166382                        *)
(*   C_nn = -0.0295085436    /Atil = -0.62005148                        *)
(*   C_kk = +0.0219278671    /Atil = +0.46076169                        *)
(*   C_rk = C_kr = -0.0038479679   /Atil = -0.08085584                  *)
(*   C_rn = C_nr = C_nk = C_kn = 0                                      *)
(*                                                                      *)
(*   tr C / Atil = +0.60237403                                          *)
(*   concurrence = 0.4212385005   -> ENTANGLED                          *)
(*   rho eigenvalues = {0.7106192488, 0.1803879332,                     *)
(*                      0.0994064930, 0.0095863249}                     *)
(* ==================================================================== *)


(* ==================================================================== *)
(* 17.  FEYNMAN DIAGRAMS for gg -> t tbar                               *)
(*  Section 1: helpers + LO      Section 2: NLO                         *)
(*  Pure Graphics, no packages needed.                                  *)
(* ==================================================================== *)


(* ==================================================================== *)
(*  1.  HELPERS                                                         *)
(* ==================================================================== *)

Clear[gluon, fermion, gluonLoop, fermionLoop, lbl, frame];

(* gluon = coiled spring; kAmp > 1 makes the loops close *)
gluon[p1_, p2_, n_: 7, kAmp_: 1.45] := Module[{d, len, u, v, a},
  d = N[p2 - p1]; len = Norm[d]; u = d/len; v = {-u[[2]], u[[1]]};
  a = kAmp len/(2 Pi n);
  {Thickness[0.005],
   Line[Table[p1 + u (s len - a Sin[2 Pi n s]) + v (a Cos[2 Pi n s] - a),
              {s, 0, 1, 0.004}]]}];

fermion[p1_, p2_] :=
  {Thickness[0.005], Arrowheads[{{0.045, 0.55}}], Arrow[{N[p1], N[p2]}]};

gluonLoop[c_, r_] := {Thickness[0.005],
  Line[Table[N[c] + (r + 0.06 Sin[12 t]) {Cos[t], Sin[t]}, {t, 0, 2 Pi, 0.01}]]};

fermionLoop[c_, r_] := Table[{Thickness[0.005], Arrowheads[{{0.045, 0.6}}],
   Arrow[Table[N[c] + r {Cos[t], Sin[t]}, {t, t0, t0 + Pi, Pi/60}]]},
   {t0, {0, Pi}}];

lbl[txt_, pos_, sz_: 15] := Text[Style[txt, sz, Black], N[pos]];

frame[body_, title_, sub_] := Graphics[
  {body, lbl[Style[title, Bold, 14], {-0.3, -2.3}],
         lbl[Style[sub, Italic, 11], {-0.3, -2.7}]},
  PlotRange -> {{-2.6, 2.1}, {-3.0, 1.8}}, ImageSize -> 240];


(* ==================================================================== *)
(*  2.  LEADING ORDER  --  three tree diagrams                          *)
(* ==================================================================== *)

loT = frame[{
   gluon[{-2.1, 1.1}, {-0.6, 1.1}], gluon[{-2.1, -1.1}, {-0.6, -1.1}],
   fermion[{-0.6, -1.1}, {-0.6, 1.1}],
   fermion[{-0.6, 1.1}, {1.0, 1.1}], fermion[{1.0, -1.1}, {-0.6, -1.1}],
   lbl[Subscript["g", 1], {-2.35, 1.1}], lbl[Subscript["g", 2], {-2.35, -1.1}],
   lbl["t", {1.25, 1.1}], lbl[OverBar["t"], {1.25, -1.1}],
   lbl[Subscript["q", "T"], {-0.25, 0}, 13]},
  "t-channel", "c1 = T^a T^b"];

loU = frame[{
   fermion[{-0.6, -1.1}, {-0.6, 1.1}],
   fermion[{-0.6, 1.1}, {1.0, 1.1}], fermion[{1.0, -1.1}, {-0.6, -1.1}],
   gluon[{-2.1, 1.1}, {-0.6, -1.1}], gluon[{-2.1, -1.1}, {-0.6, 1.1}],
   lbl[Subscript["g", 1], {-2.35, 1.1}], lbl[Subscript["g", 2], {-2.35, -1.1}],
   lbl["t", {1.25, 1.1}], lbl[OverBar["t"], {1.25, -1.1}],
   lbl[Subscript["q", "U"], {-0.28, 0}, 13]},
  "u-channel", "c2 = T^b T^a   (gluons crossed)"];

loS = frame[{
   gluon[{-2.1, 1.1}, {-1.0, 0}], gluon[{-2.1, -1.1}, {-1.0, 0}],
   gluon[{-1.0, 0}, {0.3, 0}],
   fermion[{0.3, 0}, {1.4, 1.1}], fermion[{1.4, -1.1}, {0.3, 0}],
   lbl[Subscript["g", 1], {-2.35, 1.1}], lbl[Subscript["g", 2], {-2.35, -1.1}],
   lbl["t", {1.65, 1.1}], lbl[OverBar["t"], {1.65, -1.1}],
   lbl[Superscript["g", "*"], {-0.35, 0.35}, 13]},
  "s-channel", "c3 = f^abc T^c"];

Print[Style["LO   gg -> t tbar", Bold, 16]];
Print[Row[{loT, loU, loS}]];


(* ==================================================================== *)
(*  3.  NLO  --  virtual topology classes                               *)
(* ==================================================================== *)

nloSE = frame[{
   gluon[{-2.2, 0}, {-1.0, 0}], gluonLoop[{-0.35, 0}, 0.62],
   gluon[{0.3, 0}, {1.5, 0}]},
  "self-energy", "gluon or top line"];

nloVtx = frame[{
   gluon[{-2.2, 0}, {-0.9, 0}],
   fermion[{-0.9, 0}, {0.7, 1.1}], fermion[{0.7, -1.1}, {-0.9, 0}],
   gluon[{0.02, 0.62}, {0.02, -0.62}, 5],
   lbl["t", {0.95, 1.1}], lbl[OverBar["t"], {0.95, -1.1}]},
  "vertex correction", "g t tbar, or three-gluon"];

nloBox = frame[{
   gluon[{-2.2, 1.5}, {-0.55, 0.55}], gluon[{-2.2, -1.5}, {-0.55, -0.55}],
   fermion[{-0.55, -0.55}, {-0.55, 0.55}], fermion[{-0.55, 0.55}, {0.55, 0.55}],
   fermion[{0.55, 0.55}, {0.55, -0.55}], fermion[{0.55, -0.55}, {-0.55, -0.55}],
   fermion[{0.55, 0.55}, {1.6, 1.4}], fermion[{1.6, -1.4}, {0.55, -0.55}],
   lbl["t", {1.85, 1.4}], lbl[OverBar["t"], {1.85, -1.4}]},
  "box", "four propagators in the loop"];

Print[Style["NLO   virtual: the three classes", Bold, 16]];
Print[Row[{nloSE, nloVtx, nloBox}]];


(* ==================================================================== *)
(*  4.  NLO  --  examples, with the integral each one produces          *)
(* ==================================================================== *)

nloGSE = frame[{
   gluon[{-2.2, 1.0}, {-1.3, 0}], gluon[{-2.2, -1.0}, {-1.3, 0}],
   gluon[{-1.3, 0}, {-0.75, 0}], fermionLoop[{-0.15, 0}, 0.6],
   gluon[{0.45, 0}, {1.0, 0}],
   fermion[{1.0, 0}, {1.9, 1.0}], fermion[{1.9, -1.0}, {1.0, 0}]},
  "gluon self-energy", "B0(s; mt, mt)"];

nloTSE = frame[{
   gluon[{-2.2, 1.1}, {-0.7, 1.1}], gluon[{-2.2, -1.1}, {-0.7, -1.1}],
   fermion[{-0.7, -1.1}, {-0.7, -0.3}], gluonLoop[{-0.7, 0.05}, 0.34],
   fermion[{-0.7, 0.42}, {-0.7, 1.1}],
   fermion[{-0.7, 1.1}, {1.0, 1.1}], fermion[{1.0, -1.1}, {-0.7, -1.1}]},
  "top self-energy", "B0(mt^2; 0, mt)"];

nloVC = frame[{
   gluon[{-2.2, 1.1}, {-0.7, 1.1}], gluon[{-2.2, -1.1}, {-0.7, -1.1}],
   fermion[{-0.7, -1.1}, {-0.7, 1.1}],
   fermion[{-0.7, 1.1}, {1.0, 1.1}], fermion[{1.0, -1.1}, {-0.7, -1.1}],
   gluon[{-0.7, 0.55}, {0.35, 1.1}, 5]},
  "g t tbar vertex", "C0(mt^2, mt^2, s; mt, 0, mt)"];

nloXB = frame[{
   gluon[{-2.2, 1.5}, {-0.55, -0.55}], gluon[{-2.2, -1.5}, {-0.55, 0.55}],
   fermion[{-0.55, -0.55}, {-0.55, 0.55}], fermion[{-0.55, 0.55}, {0.55, 0.55}],
   fermion[{0.55, 0.55}, {0.55, -0.55}], fermion[{0.55, -0.55}, {-0.55, -0.55}],
   fermion[{0.55, 0.55}, {1.6, 1.4}], fermion[{1.6, -1.4}, {0.55, -0.55}],
   lbl["t", {1.85, 1.4}], lbl[OverBar["t"], {1.85, -1.4}]},
  "crossed box", "D0(...; s, u; ...)"];

Print[Style["NLO   representative diagrams", Bold, 16]];
Print[Row[{nloGSE, nloTSE}]];
Print[Row[{nloVC, nloXB}]];


(* ==================================================================== *)
(*  5.  NLO  --  scalar integral basis after PV reduction               *)
(* ==================================================================== *)

box4[t_, s_, n_] := Graphics[{Thickness[0.005], t,
   lbl[Style[s, Bold, 17], {0, -1.3}], lbl[n, {0, -1.75}, 11]},
  PlotRange -> {{-1.6, 1.6}, {-2.0, 1.7}}, ImageSize -> 180];

intA0 = box4[{Circle[{0, 0.75}, 0.55], Line[{{0, 0.2}, {0, -0.7}}]},
             Subscript["A", 0], "1 propagator"];
intB0 = box4[{Circle[{0, 0.3}, 0.55], Line[{{-1.35, 0.3}, {-0.55, 0.3}}],
              Line[{{0.55, 0.3}, {1.35, 0.3}}]},
             Subscript["B", 0], "2 propagators"];
intC0 = box4[{Line[{{-0.55, -0.1}, {0.55, -0.1}, {0, 0.85}, {-0.55, -0.1}}],
              Line[{{-0.55, -0.1}, {-1.25, -0.75}}],
              Line[{{0.55, -0.1}, {1.25, -0.75}}],
              Line[{{0, 0.85}, {0, 1.55}}]},
             Subscript["C", 0], "3 propagators"];
intD0 = box4[{Line[{{-0.45, -0.15}, {-0.45, 0.75}, {0.45, 0.75},
                    {0.45, -0.15}, {-0.45, -0.15}}],
              Line[{{-0.45, 0.75}, {-1.15, 1.45}}],
              Line[{{0.45, 0.75}, {1.15, 1.45}}],
              Line[{{-0.45, -0.15}, {-1.15, -0.85}}],
              Line[{{0.45, -0.15}, {1.15, -0.85}}]},
             Subscript["D", 0], "4 propagators"];

Print[Style["Scalar basis after Passarino-Veltman reduction", Bold, 16]];
Print[Row[{intA0, intB0, intC0, intD0}]];
Print[Style["No pentagons: four external legs cannot make a 5-point function.",
            12]];


(* ==================================================================== *)
(*  6.  NLO  --  real emission (2 -> 3, needed for IR finiteness)       *)
(* ==================================================================== *)

realEm = frame[{
   gluon[{-2.2, 1.1}, {-0.7, 1.1}], gluon[{-2.2, -1.1}, {-0.7, -1.1}],
   fermion[{-0.7, -1.1}, {-0.7, 1.1}],
   fermion[{-0.7, 1.1}, {1.0, 1.1}], fermion[{1.0, -1.1}, {-0.7, -1.1}],
   gluon[{0.25, 1.1}, {1.3, 0.25}, 5],
   lbl["t", {1.25, 1.35}], lbl[OverBar["t"], {1.25, -1.1}],
   lbl["g", {1.55, 0.15}]},
  "real emission", "gg -> t tbar g"];

Print[Style["NLO   real emission", Bold, 16]];
Print[realEm];


(* --- optional export ------------------------------------------------- *)
(*
Export["LO_diagrams.pdf", Row[{loT, loU, loS}]];
Export["NLO_classes.pdf", Row[{nloSE, nloVtx, nloBox}]];
Export["NLO_basis.pdf",   Row[{intA0, intB0, intC0, intD0}]];
*)


(* ==================================================================== *)
(*   TASK 1 : towards NLO                                               *)
(*     (a) catalogue of the one-loop scalar integrals for gg -> t tbar  *)
(*     (b) analytic B0(p^2, m0^2, m1^2), validated                      *)
(*                                                                      *)
(*   Append this AFTER the LO sections of LOQCDGOAL1_rewrite.m.         *)
(*   It is self-contained: it does not touch any LO definition, and it  *)
(*   uses its own symbols (prefix "lp") so nothing collides.            *)
(*                                                                      *)
(*   CONVENTIONS                                                        *)
(*     B0 = Delta_eps + B0fin,   Delta_eps = 1/eps - EulerGamma + Log[4 Pi] *)
(*     propagators  1/(q^2 - m^2 + I eps)                                *)
(*     => Im B0 > 0 above threshold  (optical theorem)                  *)
(*     B0fin depends on mu^2.  ALWAYS quote mu with a number.           *)
(* ==================================================================== *)

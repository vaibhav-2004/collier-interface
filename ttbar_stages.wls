#!/usr/bin/env wolframscript
(* ====================================================================
   gg -> t tbar : STAGES 1 AND 2, run in WSL

   Stage 1   FeynArts    enumerate the diagrams
   Stage 2   FormCalc    Dirac + colour algebra, PV reduction
                         (calls FORM and ReadForm, external programs)

   WHY THIS IS A SEPARATE FILE FROM THE NOTEBOOK
   FormCalc cannot run in Windows Mathematica: MathLink will not launch
   ReadForm across the Windows/WSL boundary.  It fails either with
   ReadForm::notcompiled or by hanging at load, and the nasty part is
   that a failed load leaves CreateFeynAmp working while CalcFeynAmp
   silently does nothing - you get raw uncontracted SUNF[Glu1,Glu2,Glu5]
   and an EMPTY integral list, with no error.  Hence the explicit load
   check below.

   Writing a .mx here and reading it from Windows Mathematica is fine:
   that is a file read, not a process launch.

   RUN:
       cd ~/collier_work
       ./ttbar_stages.wls 2>&1 | tee ttbar_stages_output.txt

   PREREQUISITES (one-time, already done):
     - Wolfram Engine installed in WSL and activated
     - sudo ln -s /usr/local/bin/WolframKernel /usr/local/bin/math
       (FormCalc's ./compile looks for a command called 'math')
     - cd FormCalc-9.10 && ./compile    builds native Linux ReadForm
     - FORM ships with FormCalc (bin/Linux-x86-64/form); no apt needed
   ==================================================================== *)

base = "/home/vaibhav/collier_work/";

(* LOAD ORDER MATTERS: FeynArts first, then FormCalc.  Setting
   $FeynArtsDir in the Global` context beforehand causes a shadowing
   abort - FormCalc finds FeynArts on its own.                        *)

Get[base <> "FeynArts-3.11/FeynArts.m"];
Get[base <> "FormCalc-9.10/FormCalc.m"];

(* the check that catches a silent failure *)
faOK = Length[DownValues[CreateFeynAmp]] > 0;
fcOK = Length[DownValues[CalcFeynAmp]]  > 0;
Print["FeynArts loaded : ", faOK];
Print["FormCalc loaded : ", fcOK];
If[! (faOK && fcOK),
   Print["\n** ABORTING: a package did not load.  Everything downstream"];
   Print["   would silently return nothing.  Check the paths above."];
   Exit[1]];


(* ====================================================================
   STAGE 1 -- FeynArts
   ==================================================================== *)

(* Fields:  V[5]       the gluon
            F[3,{3}]   third-generation up-type quark = top
           -F[3,{3}]   the antitop                                    *)

Print["\n--- stage 1: tree (expect 3 Classes insertions) ---"];

topsTree = CreateTopologies[0, 2 -> 2];
insTree  = InsertFields[topsTree,
             {V[5], V[5]} -> {F[3, {3}], -F[3, {3}]},
             Model -> "SMQCD", InsertionLevel -> {Classes}];

(* SMQCD inherits the FULL SM content, so without this filter you get
   112 diagrams including H, G0, G+-, W, Z, photon and b.  Those are
   O(alpha), not O(alpha_s).  Kept: gluons V[5], gluon ghosts U[5],
   quarks F[3] and F[4] - the light quarks are needed because they run
   in the gluon self-energy loop.                                      *)

qcdOnly = {S[1], S[2], S[3],    (* Higgs, G0, G+-     *)
           V[1], V[2], V[3],    (* photon, Z, W       *)
           F[1], F[2]};         (* neutrinos, leptons *)

(* Tadpoles vanish in dim reg for massless fields and are absorbed into
   renormalisation.  External-leg (wave function) corrections go into
   the field renormalisation constants, so drawing them would double
   count.                                                              *)

Print["\n--- stage 1: one loop, QCD only (expect 34) ---"];

topsLoop = CreateTopologies[1, 2 -> 2,
             ExcludeTopologies -> {Tadpoles, WFCorrections}];
insLoop  = InsertFields[topsLoop,
             {V[5], V[5]} -> {F[3, {3}], -F[3, {3}]},
             Model -> "SMQCD", InsertionLevel -> {Classes},
             ExcludeParticles -> qcdOnly];

(* Read the counts off the "in total: N Classes insertions" lines that
   InsertFields prints.  Do NOT use DiagramExtract[ins, All] - it
   returns 0 in FeynArts 3.11.
   The 34 include GHOST loops (u_g) and LIGHT-QUARK loops (u_m, d_m);
   the hand catalogue of ~15 integrals had neither.
   No pentagons: four external legs cannot make a 5-point function.    *)

DumpSave[base <> "feynarts_stage1.mx", {insTree, insLoop}];
Print["\nsaved ", base, "feynarts_stage1.mx"];

(* To LOOK at the diagrams, Get that .mx in Windows Mathematica and
   Paint[insLoop, ColumnsXRows -> {4,5}, PaintLevel -> {Classes}].
   There is no front end under wolframscript.                          *)


(* ====================================================================
   STAGE 2 -- FormCalc

   THE OPTIONS THAT MATTER (Marcel's correction):

     FermionChains -> VA    vector/axial basis rather than chiral
                            projectors; fewer, simpler structures

     PaVeReduce -> True     reduces the TENSOR coefficients all the way
                            down to scalars.  Without it the same run
                            gives 222 objects with 38 coefficient types
                            (bb0, bb1, bb00, cc00, dd123 ...) and the
                            Collier driver would need its tensor
                            routines B_cll / C_cll / D_cll.  With it,
                            only bb0 / cc0 / dd0 survive - which the
                            existing scalar driver already computes and
                            the existing reference points validate.
   ==================================================================== *)

Print["\n--- stage 2: tree amplitude ---"];
ampTree = CalcFeynAmp[CreateFeynAmp[insTree],
            FermionChains -> VA, PaVeReduce -> True];

Print["\n--- stage 2: one-loop amplitude (slow) ---"];
ampLoop = CalcFeynAmp[CreateFeynAmp[insLoop],
            FermionChains -> VA, PaVeReduce -> True];

(* Abbr[] is NOT preserved by DumpSave - it reads FormCalc's internal
   tables, so query it HERE, in the session that computed the
   amplitude.  This is why an earlier attempt that reloaded the .mx and
   then called Abbr[] got 0 and an empty integral list.                *)

abbr = Abbr[];
Print["\nAbbr[] entries: ", Length[abbr]];

(* FormCalc's heads are B0i, C0i, D0i - NOT B0, C0, D0.  Searching for
   the latter finds nothing.  The integrals live in the abbreviation
   table as well as in the amplitude, so search BOTH.                  *)

ints = Union[Cases[{ampLoop, abbr},
                   (A0 | A00 | B0i | C0i | D0i)[__], Infinity]];

types = Union[Cases[ints, (B0i | C0i | D0i)[t_, ___] :> t, Infinity]];

colours = Union[Cases[{ampLoop, abbr}, SUNT[__] | SUNF[__], Infinity]];


(* ====================================================================
   REPORT  --  this is Task 1a
   ==================================================================== *)

Print["\n===================================================="];
Print[" TASK 1a : ", Length[ints], " distinct scalar integrals"];
Print["===================================================="];
Scan[Print["  ", #] &, ints];

Print["\n--- coefficient types (want only bb0, cc0, dd0) ---"];
Print[types];

Print["\n--- colour structures ---"];
Scan[Print["  ", #] &, colours];
Print["
  SUNT[Glu1, Glu2, Col3, Col4] = T^a T^b    = c1
  SUNT[Glu2, Glu1, Col3, Col4] = T^b T^a    = c2
  SUNT[Col3, Col4]             = delta_ij   <- NOT in the tree basis
  SUNT[Glu1, Glu2, 0, 0]       = Tr(T^a T^b), a closed quark loop

  At tree level the basis was c1, c2, f^abc T^c.  But
      f^abc T^c = -i (c1 - c2)
  so the tree basis spans only 2 dimensions as matrices.  delta^ab * 1
  IS independent, so the one-loop basis is {c1, c2, delta^ab 1} and the
  colour matrix must be recomputed - see PART 8 of the notebook."];

DumpSave[base <> "formcalc_stage2_pave.mx",
         {ampTree, ampLoop, abbr, ints, types, colours}];
Print["\nsaved ", base, "formcalc_stage2_pave.mx"];
Print["Now open the notebook in Windows Mathematica - PART 7 reads it."];

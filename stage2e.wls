#!/usr/bin/env wolframscript

base = "/home/vaibhav/collier_work/";
Get[base <> "FeynArts-3.11/FeynArts.m"];
Get[base <> "FormCalc-9.10/FormCalc.m"];
Get[base <> "feynarts_stage1.mx"];

Print["computing one loop with PaVeReduce..."];
amp = CalcFeynAmp[CreateFeynAmp[insLoop],
        FermionChains -> VA, PaVeReduce -> True];

abbr = Abbr[];
ints = Union[Cases[{amp, abbr},
          (A0|A00|B0|B0i|C0|C0i|D0|D0i)[__], Infinity]];

Print["\n=== integrals: ", Length[ints], " ==="];
Scan[Print["  ", #] &, ints];

Print["\n=== heads appearing ==="];
Print[Union[Cases[ints, h_[___] :> h, Infinity]]];

Print["\n=== any tensor coefficients left? ==="];
Print[Union[Cases[ints, (B0i|C0i|D0i)[t_, ___] :> t, Infinity]]];

Print["\n=== colour structures ==="];
Print[Union[Cases[{amp, abbr}, SUNT[__]|SUNF[__], Infinity]]];

DumpSave[base <> "formcalc_stage2_pave.mx", {amp, abbr, ints}];

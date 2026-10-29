#!/usr/bin/env wolframscript

Get["/home/vaibhav/collier_work/FeynArts-3.11/FeynArts.m"];
Get["/home/vaibhav/collier_work/FormCalc-9.10/FormCalc.m"];
Get["/home/vaibhav/collier_work/feynarts_stage1.mx"];

amp = CalcFeynAmp[CreateFeynAmp[insLoop]];

ints = Union[Cases[{amp, Abbr[]}, (A0|A00|B0i|C0i|D0i)[__], Infinity]];

Print["\n=== TASK 1a: ", Length[ints], " distinct loop integrals ==="];
Scan[Print["  ", #] &, ints];

Print["\n=== which coefficient types appear ==="];
Print[Union[Cases[ints, (B0i|C0i|D0i)[t_, ___] :> t, Infinity]]];

Print["\n=== colour structures ==="];
Print[Union[Cases[{amp, Abbr[]}, SUNT[__]|SUNF[__], Infinity]]];

DumpSave["/home/vaibhav/collier_work/formcalc_stage2.mx",
         {amp, ints}];

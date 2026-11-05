#!/usr/bin/env wolframscript

base = "/home/vaibhav/collier_work/";

(* do NOT set $FeynArtsDir beforehand - a script that did so aborted *)
Get[base <> "FeynArts-3.11/FeynArts.m"];

(* tree level: expect 3 (s, t, u channel) *)
topsTree = CreateTopologies[0, 2 -> 2];
insTree  = InsertFields[topsTree,
             {V[5], V[5]} -> {F[3, {3}], -F[3, {3}]},
             Model -> "SMQCD", InsertionLevel -> {Classes}];

(* one loop, QCD only: SMQCD carries the full SM content, so without
   ExcludeParticles you get 112 diagrams with H, G0, G+-, W, Z, photon, b.
   Tadpoles and external-leg corrections are excluded: tadpoles are
   absorbed by renormalisation, external legs go into field
   renormalisation constants (drawing them would double count). *)
qcdOnly = {S[1], S[2], S[3], V[1], V[2], V[3], F[1], F[2]};

topsLoop = CreateTopologies[1, 2 -> 2,
             ExcludeTopologies -> {Tadpoles, WFCorrections}];
insLoop  = InsertFields[topsLoop,
             {V[5], V[5]} -> {F[3, {3}], -F[3, {3}]},
             Model -> "SMQCD", InsertionLevel -> {Classes},
             ExcludeParticles -> qcdOnly];

DumpSave[base <> "feynarts_stage1.mx", {insTree, insLoop}];
Print["\nsaved ", base, "feynarts_stage1.mx"];

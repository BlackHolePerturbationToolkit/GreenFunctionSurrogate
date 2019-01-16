(* ::Package:: *)

BeginPackage["GreenFunctionSurrogate`EvaluateSurrogate`"];


SurrogateData::usage =
 "SurrogateData[\!\(\*
StyleBox[\"sur\",\nFontSlant->\"Italic\"]\)] represents a surrogate with properties and data stored in \!\(\*
StyleBox[\"sur\",\nFontSlant->\"Italic\"]\).";
ReadSurrogate::usage =
 "ReadSurrogate[\!\(\*
StyleBox[\"file\",\nFontSlant->\"Italic\"]\)] creates a SurrogateData object representing the surrogate stored in \!\(\*
StyleBox[\"file\",\nFontSlant->\"Italic\"]\).";
GreenFunctionSurrogate::usage =
 "GreenFunctionSurrogate[\!\(\*
StyleBox[\"s\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"rstar\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"rstar0\",\nFontSlant->\"Italic\"]\)] gives the Green function for the spin \!\(\*
StyleBox[\"s\",\nFontSlant->\"Italic\"]\) Regge-Wheeler "<>
 "with eigenvalue \!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\) and at field point (\!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"rstar\",\nFontSlant->\"Italic\"]\)) for a base point at (0, \!\(\*
StyleBox[\"rstar0\",\nFontSlant->\"Italic\"]\)).";


Begin["`Private`"];


(* If the h5mma is not found, then just use Mathematica's built-in HDF5 support *)
$h5mma = If[Quiet[Get["h5mma`"], {Needs::nocont, Get::noopen}]===$Failed, False, True];
If[$h5mma, SetOptions[ImportHDF5, Turbo->True]];

ReadHDF5[file_String, opts_:"Datasets"] :=
  If[$h5mma, ImportHDF5[file, opts], Import[file, opts]];


surData[s_Integer, file_] := surData[s, file] = ReadSurrogate[file];


sur[s_Integer] := surData[s, $SurrogateFile[s]];


$SurrogateFile[0] = FileNameJoin[{FileNameDrop[FindFile["GreenFunctionSurrogate`"], -2], "SurrogateData", "G_s0", "G_sur_s0_sigma0.1_h0.01_t240_hyp.h5"}];
$SurrogateFile[2] = FileNameJoin[{FileNameDrop[FindFile["GreenFunctionSurrogate`"], -2], "SurrogateData", "G_s2", , "G_sur_s2_sigma0.1_h0.01_t240_hyp.h5"}];


(* Read in data for surrogate *)
ReadSurrogate[file_String] :=
 Module[{eim, B, s, lmax, \[Sigma], rstar, rstar0, times},
  (* Load lightweight data *)
  s = ReadHDF5[file, {"Datasets", "s"}];
  lmax = ReadHDF5[file, {"Datasets", "lMax"}];
  \[Sigma] = ReadHDF5[file, {"Datasets", "sigma"}];
  rstar = ReadHDF5[file, {"Datasets", "fieldPoints"}];
  rstar0 = ReadHDF5[file, {"Datasets", "basePoints"}];
  times = ReadHDF5[file, {"Datasets", "times"}];

  (* Create functions which load the surrogate data on demand (with caching) *)
  eim[l_] := eim[l] = ReadHDF5[file,{ "Datasets", "l="<>ToString[l]<>"/eimData"}];
  B[l_] := B[l]= Transpose[ReadHDF5[file, {"Datasets", "l="<>ToString[l]<>"/B"}]];

  (* Return a SurrogateData object representing the surrogate model *)
  SurrogateData[<|"DataFile"->file,
    "s"->s, "lmax"->lmax, "sigma"->\[Sigma], "rstar"->rstar, "rstar0"->rstar0, "times"->times,
    "B"->B, "eim"->eim|>]
];


Format[SurrogateData[assoc_Association]] := "SurrogateData["<>assoc["DataFile"]<>",<<>>]";

SurrogateData[assoc_Association][y_String] := assoc[y];

partPattern = (_Integer | {__Integer} | Span[(_Integer|All)] | Span[(_Integer|All), (_Integer|All)] |
  Span[(_Integer|All), (_Integer|All), (_Integer|All)] | All);

SurrogateData[assoc_Association][l_Integer, ti_, rstari_, rstar0i_] :=
 Module[{},
  If[!MatchQ[ti, partPattern] || !MatchQ[rstari, partPattern] || !MatchQ[rstar0i, partPattern],
    Return[$Failed];
  ];
  
  If[l>assoc["lmax"] || l<assoc["s"],
    Return[$Failed];
  ];

  assoc["B"][l][[ti]].assoc["eim"][l][[All, rstari,rstar0i]]
];


(* FIXME: Add support for higher order interpolation *)
(* FIXME: Use off-centred stencil when within one grid cell of the r=r' kink *)
(* FIXME: Is $MachineEpsilon too stringent here? *)
findPatch[x_, grid_] :=
 Module[{nearest, index, patch, patchcoords, nearestWithPadding},
  nearest = First[Nearest[grid -> All, x]];
  If[nearest["Distance"] < $MachineEpsilon,
    patch = index = nearest["Index"];
    patchcoords = grid[[index]];,
    (* Always choose the nearest point to the left of the specified point *)
    If[x > nearest["Element"],
      index = nearest["Index"];,
      index = nearest["Index"]-1;
    ];
    nearestWithPadding = Min[Max[index, 2], Length[grid]-2];
    patch = Span[nearestWithPadding-1, nearestWithPadding+2];
    patchcoords = grid[[List@@patch]];,
    $Failed
  ];
  {index, patch, patchcoords}
];


GreenFunctionSurrogate[s_, l_,
       t : (_?NumericQ | {_?NumericQ ..}),
   rstar : (_?NumericQ | {_?NumericQ ..}),
  rstar0 : (_?NumericQ | {_?NumericQ ..})] :=
 Module[{tRet, ti, tp, trange, rstari, rstarp, rstarrange, rstar0i, rstar0p, rstar0range, T, Rstar, Rstar0, interp},
  (* Account for time retardation *)
  tRet = t - Abs[rstar - rstar0];

  (* Green function is zero for points outside the lightcone *)
  If[tRet < 0, Return[0.]];

  {ti, tp, trange} = findPatch[tRet, sur[s]["times"]];
  {rstari, rstarp, rstarrange} = findPatch[rstar, sur[s]["rstar"]];
  {rstar0i, rstar0p, rstar0range} = findPatch[rstar0, sur[s]["rstar0"]];

  (* Handle the cases where one or more coordinates coincides with a grid point *)
  If[Head[tp]=!=Span,
    trange = Nothing;
    T = Nothing;,
    T = tRet;
  ];

  If[Head[rstarp]=!=Span,
    rstarrange = Nothing;
    Rstar = Nothing;,
    Rstar=rstar;
  ];

  If[Head[rstar0p]=!=Span,
    rstar0range = Nothing;
    Rstar0 = Nothing;,
    Rstar0 = rstar0;
  ];

  If[SameQ[Rstar0, Rstar, T, Nothing],
    sur[s][l, tp, rstarp, rstar0p],
    interp = ListInterpolation[sur[s][l, tp, rstarp, rstar0p], {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0}
  ]
];


End[];
EndPackage[];


(* Add h5mma to $ContextPath since Get[] inside `Private` does not do so. *)
If[SimulationTools`ReadHDF5`Private`$h5mma,
  If[!MemberQ[$ContextPath, "h5mma`"], AppendTo[$ContextPath, "h5mma`"]];
];

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
$SurrogateFile[2] = FileNameJoin[{FileNameDrop[FindFile["GreenFunctionSurrogate`"], -2], "SurrogateData", "G_s2", "G_sur_s2_sigma0.1_h0.01_t240_hyp.h5"}];


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
  If[Block[{Internal`$SameQTolerance = 4}, SameQ[nearest["Element"], x]],
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
 Module[{tRet, ti, tp, trange, rstari, rstarp, rstarrange, rstar0i, rstar0p, rstar0range, T, Rstar, Rstar0, data, interp},
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

  Which[
  (* We have this grid point in our surrogate data, no need to interpolate *)
  SameQ[Rstar0, Rstar, T, Nothing],
    Sow[{"Exact: ", {tp, rstarp, rstar0p}}];
    sur[s][l, tp, rstarp, rstar0p],

  (* Either we're far enough away from the r=r' kink to use a standard interpolation
     stencil, or we don't need to interpolate in rstar or rstar0 as they are both
     on our surrogate data grid. *)
  Abs[rstari - rstar0i] >= 3 || SameQ[Rstar0, Rstar, Nothing],
    Sow[{"Interp: ", {tp, rstarp, rstar0p}, {trange, rstarrange, rstar0range}, {T, Rstar, Rstar0}}];
    interp = ListInterpolation[sur[s][l, tp, rstarp, rstar0p], {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0},

  (* We're close to the r=r' kink; use a skewed stencil to avoid the discontinuity.
     There are 5 subcases depending on which edge we're nearest to: *)

  (* Case 1: r>r0 and r << rmax - skew grid in +r direction *)
  (0 < rstar - rstar0) && (rstari - rstar0i < 3) && rstari < Length[sur[s]["rstar"]]/2,
    (* NB: in the following being on an r grid point doesn't help,
       as it won't be on a skewed gridpoint. Conversely, it's possible
       that we're on a skewed grid point and miss it, but that's an optimisation
       that can be added later. *)
    Rstar = rstar;
    Which[
    SameQ[T, Rstar0, Nothing], (* Have to interpolate in r only *)
      data = Table[sur[s][l, ti, j, rstar0i], {j, rstar0i, rstar0i + 3}];
      rstarrange = sur[s]["rstar"][[{rstar0i, rstar0i + 3}]];,
    SameQ[Rstar0, Nothing], (* Have to interpolate in r and t *)
      data = Table[sur[s][l, i, j, rstar0i], {i, tp[[1]], tp[[2]]}, {j, rstar0i, rstar0i + 3}];
      rstarrange = sur[s]["rstar"][[{rstar0i, rstar0i + 3}]];,
    SameQ[T, Nothing] || SameQ[T, Rstar, Nothing], (* Have to interpolate in r and r0 *)
      data = Table[sur[s][l, ti, j, k], {j, k, k+3}, {k, rstar0p[[1]], rstar0p[[2]]}];
      rstarrange = {rstar0, rstar0 + (sur[s]["rstar"][[rstari+3]]-sur[s]["rstar"][[rstari]])};,
    SameQ[Rstar, Nothing] || True, (* Have to interpolate in t, r and r0 *)
      data = Table[sur[s][l, i, j, k], {i, tp[[1]], tp[[2]]}, {j, k, k+3}, {k, rstar0p[[1]], rstar0p[[2]]}];
      rstarrange = {rstar0, rstar0 + (sur[s]["rstar"][[rstari+3]]-sur[s]["rstar"][[rstari]])};
    ];
    Sow[{"Interp skew 1: ", data, {trange, rstarrange, rstar0range}, {T, Rstar, Rstar0}}];
    interp = ListInterpolation[data, {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0},

  (* Case 2: r>r0 and r >> 0 - skew in -r0 direction *)
  (0 < rstar - rstar0) && (rstari - rstar0i < 3) && rstari >= Length[sur[s]["rstar"]]/2,
    (* NB: in the following being on an r0 grid point doesn't help,
       as it won't be on a skewed gridpoint. Conversely, it's possible
       that we're on a skewed grid point and miss it, but that's an optimisation
       that can be added later. *)
    Rstar0 = rstar0;
    Which[
    SameQ[T, Rstar, Nothing], (* Have to interpolate in r0 only *)
      data = Table[sur[s][l, ti, rstari, j], {j, rstari - 3, rstari}];
      rstar0range = sur[s]["rstar0"][[{rstari - 3, rstari}]];,
    SameQ[Rstar, Nothing], (* Have to interpolate in r0 and t *)
      data = Table[sur[s][l, i, rstari, j], {i, tp[[1]], tp[[2]]}, {j, rstari - 3, rstari}];
      rstar0range = sur[s]["rstar0"][[{rstari - 3, rstari}]];,
    SameQ[T, Nothing] || SameQ[T, Rstar0, Nothing], (* Have to interpolate in r and r0 *)
      data = Table[sur[s][l, ti, j, k], {j, rstarp[[1]], rstarp[[2]]}, {k, j - 3, j}];
      rstar0range = {rstar - (sur[s]["rstar0"][[rstar0i]] - sur[s]["rstar0"][[rstar0i - 3]]), rstar};,
    SameQ[Rstar0, Nothing] || True, (* Have to interpolate in t, r and r0 *)
      data = Table[sur[s][l, i, j, k], {i, tp[[1]], tp[[2]]}, {j, rstarp[[1]], rstarp[[2]]}, {k, j - 3, j}];
      rstar0range = {rstar - (sur[s]["rstar0"][[rstar0i]] - sur[s]["rstar0"][[rstar0i - 3]]), rstar};
    ];
    Sow[{"Interp skew 2: ", data, {trange, rstarrange, rstar0range}, {T, Rstar, Rstar0}}];
    interp = ListInterpolation[data, {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0},

  (* Case 3: r<r0 and r << rmax - skew in +r0 direction*)
  (-3 < rstari - rstar0i) && (rstar - rstar0 < 0) && rstari < Length[sur[s]["rstar"]]/2,
    (* NB: in the following being on an r0 grid point doesn't help,
       as it won't be on a skewed gridpoint. Conversely, it's possible
       that we're on a skewed grid point and miss it, but that's an optimisation
       that can be added later. *)
    Rstar0 = rstar0;
    Which[
    SameQ[T, Rstar, Nothing], (* Have to interpolate in r0 only *)
      data = Table[sur[s][l, ti, rstari, j], {j, rstari, rstari + 3}];
      rstar0range = sur[s]["rstar0"][[{rstari, rstari + 3}]];,
    SameQ[Rstar, Nothing], (* Have to interpolate in r0 and t *)
      data = Table[sur[s][l, i, rstari, j], {i, tp[[1]], tp[[2]]}, {j, rstari, rstari + 3}];
      rstar0range = sur[s]["rstar0"][[{rstari, rstari + 3}]];,
    SameQ[T, Nothing] || SameQ[T, Rstar0, Nothing], (* Have to interpolate in r and r0 *)
      data = Table[sur[s][l, ti, j, k], {j, rstarp[[1]], rstarp[[2]]}, {k, j, j + 3}];
      rstar0range = {rstar, rstar + (sur[s]["rstar0"][[rstar0i + 3]] - sur[s]["rstar0"][[rstar0i]])};,
    SameQ[Rstar0, Nothing] || True, (* Have to interpolate in t, r and r0 *)
      data = Table[sur[s][l, i, j, k], {i, tp[[1]], tp[[2]]}, {j, rstarp[[1]], rstarp[[2]]}, {k, j, j + 3}];
      rstar0range = {rstar, rstar + (sur[s]["rstar0"][[rstar0i + 3]] - sur[s]["rstar0"][[rstar0i]])};
    ];
    Sow[{"Interp skew 3: ", data, {trange, rstarrange, rstar0range}, {T, Rstar, Rstar0}}];
    interp = ListInterpolation[data, {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0},

  (* Case 4: r<r0 and r >> 0 - skew in -r direction *)
  (-3 < rstari - rstar0i) && (rstar - rstar0 < 0) && rstari >= Length[sur[s]["rstar"]]/2,
    (* NB: in the following being on an r grid point doesn't help,
       as it won't be on a skewed gridpoint. Conversely, it's possible
       that we're on a skewed grid point and miss it, but that's an optimisation
       that can be added later. *)
    Rstar = rstar;
    Which[
    SameQ[T, Rstar0, Nothing], (* Have to interpolate in r only *)
      data = Table[sur[s][l, ti, j, rstar0i], {j, rstar0i - 3, rstar0i}];
      rstarrange = sur[s]["rstar"][[{rstar0i - 3, rstar0i}]];,
    SameQ[Rstar0, Nothing], (* Have to interpolate in r and t *)
      data = Table[sur[s][l, i, j, rstar0i], {i, tp[[1]], tp[[2]]}, {j, rstar0i - 3, rstar0i}];
      rstarrange = sur[s]["rstar"][[{rstar0i - 3, rstar0i}]];,
    SameQ[T, Nothing] || SameQ[T, Rstar, Nothing], (* Have to interpolate in r and r0 *)
      data = Table[sur[s][l, ti, j, k], {j, k-3, k}, {k, rstar0p[[1]], rstar0p[[2]]}];
      rstarrange = {rstar0 - (sur[s]["rstar"][[rstari]] - sur[s]["rstar"][[rstari - 3]]), rstar0};,
    SameQ[Rstar, Nothing] || True, (* Have to interpolate in t, r and r0 *)
      data = Table[sur[s][l, i, j, k], {i, tp[[1]], tp[[2]]}, {j, k-3, k}, {k, rstar0p[[1]], rstar0p[[2]]}];
      rstarrange = {rstar0 - (sur[s]["rstar"][[rstari]] - sur[s]["rstar"][[rstari - 3]]), rstar0};
    ];
    Sow[{"Interp skew 4: ", data, {trange, rstarrange, rstar0range}, {T, Rstar, Rstar0}}];
    interp = ListInterpolation[data, {trange, rstarrange, rstar0range}];
    interp @@ {T, Rstar, Rstar0},

  (* Case 5: r=r0 *)
  rstar == rstar0,
    If[SameQ[T, Nothing],
      data = Table[sur[s][l, ti, j, j], {j, rstar0p[[1]], rstar0p[[2]]}];
      Sow[{"Interp skew 5a: ", data, rstar0range, {T, Rstar, Rstar0}}];
      interp = ListInterpolation[data, rstar0range];
      interp @@ {Rstar0}
    ,
      data = Table[sur[s][l, i, j, j], {i, tp[[1]], tp[[2]]}, {j, rstar0p[[1]], rstar0p[[2]]}];
      Sow[{"Interp skew 5b: ", data, {trange, rstar0range}, {T, Rstar, Rstar0}}];
      interp = ListInterpolation[data, {trange, rstar0range}];
      interp @@ {T, Rstar0}
    ],

  True, $Failed (* We should never reach here *)
  ]
];


End[];
EndPackage[];


(* Add h5mma to $ContextPath since Get[] inside `Private` does not do so. *)
If[SimulationTools`ReadHDF5`Private`$h5mma,
  If[!MemberQ[$ContextPath, "h5mma`"], AppendTo[$ContextPath, "h5mma`"]];
];

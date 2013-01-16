(* ::Package:: *)

(* ::Title:: *)
(*AWGN channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in AWGN channels.*)
(*Copyright (C) 2012 Donagh Horgan.*)
(*Email: donaghh@rennes.ucc.ie.*)
(**)
(*This program is free software : you can redistribute it and/or modify*)
(*it under the terms of the GNU General Public License as published by*)
(*the Free Software Foundation, either version 3 of the License, or*)
(*(at your option) any later version.*)
(**)
(*This program is distributed in the hope that it will be useful,*)
(*but WITHOUT ANY WARRANTY; without even the implied warranty of*)
(*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See *)
(*COPYING for more details.*)
(**)
(*You should have received a copy of the GNU General Public License*)
(*along with this program. If not, see http://www.gnu.org/licenses.*)


(* ::Subsection:: *)
(*Version information*)


(* ::Text:: *)
(*14/12/2012*)
(*1.32*)


(* ::Subsection:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.32: Various bug fixes.*)
(*Version 1.31: Removed LowSNR option in favour of new LowSNR approximate method and amalgamated Method and Algorithm options.*)
(*Version 1.30: Added diversity reception support to the sample complexity function.*)
(*Version 1.24: Added SEC support, and removed SSC support.*)
(*Version 1.23: Moved LowSNRErrorBound to the Nakagami package.*)
(*Version 1.22: Moved help functions to Network package.*)
(*Version 1.21: Added generic help functions for consistent documentation.*)
(*Version 1.2: Recoded ProbabilityOfFalseAlarm, ProbabilityOfDetection and \[Lambda] functions, so they are easier to read.*)
(*Version 1.12: Added EGC support.*)
(*Version 1.11: Added protection for symbols.*)
(*Version 1.1: Added diversity types to functions.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["AWGN`"];


(* ::Subsection::Closed:: *)
(*PDF of the recieved energy*)


AWGNPDF;


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm*)


AWGNProbabilityOfFalseAlarm;


(* ::Subsubsection::Closed:: *)
(*Probability of detection*)


AWGNProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Threshold*)


\[Lambda];


(* ::Subsection::Closed:: *)
(*Sample complexity*)


AWGNSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<QFunction`;
<<Extras`;
<<Network`;


(* ::Subsection::Closed:: *)
(*PDF of the recieved energy*)


Options[AWGNPDF] = {Method->OptionValue[ProbabilityOfDetection,Method]};
AWGNPDF::usage="AWGNPDF[M, \[Gamma], x] evaluates the probability density function of the recieved energy at a single energy detector operating on an AWGN channel at x.
AWGNPDF[M, \[Gamma], x, n] evaluates the probability density function of the recieved energy at the fusion center of a cooperative network operating on an AWGN fading channel at x.\n\n"<>MethodHelp[AWGNPDF, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\""}];
AWGNPDF[M_,\[Gamma]_,x_,OptionsPattern[]]:=Module[{n = 1},AWGNPDF[M,\[Gamma],x,n,Method->OptionValue[Method]]]
AWGNPDF[M_,\[Gamma]_,x_,n_,OptionsPattern[]]:=Module[{method, mn},
	{method, mn} = ProcessMethod[OptionValue[Method]];

	Which[
		method == "ExactNumerical",
			PDF[NoncentralChiSquareDistribution[M n, M n \[Gamma]], x],
		method == "ApproximateNumerical",
			PDF[NormalDistribution[M n(1+\[Gamma]),Sqrt[2M n(1+2\[Gamma])]],x],
		method == "ApproximateNumericalLowSNR",
			PDF[NormalDistribution[M n(1+\[Gamma]),Sqrt[2M n]],x],
		True,
			Undefined
	]
]


(* ::Subsection:: *)
(*Decision probabilities*)


(* ::Subsubsection::Closed:: *)
(*Probability of false alarm*)


Options[AWGNProbabilityOfFalseAlarm]={Method->OptionValue[ProbabilityOfFalseAlarm,Method], DiversityType->OptionValue[ProbabilityOfFalseAlarm,DiversityType]};
AWGNProbabilityOfFalseAlarm::usage="AWGNProbabilityOfFalseAlarm[M, \[Lambda]] calculates the probability of false alarm for a single energy detector operating on an AWGN channel.
AWGNProbabilityOfFalseAlarm[M, \[Lambda], n] calculates the probability of false alarm for energy detection with diversity reception in an AWGN channel.\n\n"<>MethodHelp[AWGNProbabilityOfFalseAlarm, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\""}]<>"\n\n"<>DiversityTypeHelp[AWGNProbabilityOfFalseAlarm];
AWGNProbabilityOfFalseAlarm[M_,\[Lambda]_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AWGNProbabilityOfFalseAlarm][[All,1]]],Options[target][[All,1]]];
	AWGNProbabilityOfFalseAlarm[M,\[Lambda],n,#/.(DiversityType/.#)->"None"&[RelevantOptions[AWGNProbabilityOfFalseAlarm]]]
]
AWGNProbabilityOfFalseAlarm[M_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{diversityType = OptionValue[DiversityType], \[Gamma]t, method, mn, g},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	g[a_] := Which[
		method == "Exact" || StringTake[method, 5] == "Exact",
			GammaRegularized[M a / 2, \[Lambda] / 2],
		method == "Approximate" || StringTake[method, 11] == "Approximate" || method == "NGaussian",
			Q[(\[Lambda] - M a) / Sqrt[2M a]],
		True,
			Undefined
	];

	Which[
		diversityType == "None" || diversityType == "MRC" || diversityType == "EGC" || diversityType == "SC" || diversityType == "SEC",
			g[1],
		diversityType == "SLC",
			g[n],
		diversityType == "SLS",
			1 - (1 - g[1])^n,
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Probability of detection*)


Options[AWGNProbabilityOfDetection]={Method->OptionValue[ProbabilityOfDetection,Method], DiversityType->OptionValue[ProbabilityOfDetection,DiversityType]};
AWGNProbabilityOfDetection::usage="AWGNProbabilityOfDetection[M, \[Gamma], \[Lambda]] calculates the approximate probability of detection for a single energy detector operating on an AWGN channel.
AWGNProbabilityOfDetection[M, \[Gamma], \[Lambda], n] calculates the approximate probability of detection for energy detection with diversity reception in an AWGN channel.\n\n"<>MethodHelp[AWGNProbabilityOfDetection, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\""}]<>"\n\n"<>DiversityTypeHelp[AWGNProbabilityOfDetection];
AWGNProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AWGNProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	AWGNProbabilityOfDetection[M,\[Gamma],\[Lambda],n,#/.(DiversityType/.#)->"None"&[RelevantOptions[AWGNProbabilityOfDetection]]]
]
AWGNProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,n_,OptionsPattern[]]:=Module[{RelevantOptions, diversityType = OptionValue[DiversityType], \[Gamma]t, \[Gamma]0, method, mn, g},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[a_,b_] := Which[
		method == "ExactNumerical",
			MarcumQ[M a / 2, Sqrt[M b \[Gamma]0], Sqrt[\[Lambda]]],
		method == "ApproximateNumerical",
			Q[(\[Lambda] - M (a + b \[Gamma]0)) / Sqrt[2M (a + b 2 \[Gamma]0)]],
		method == "ApproximateNumericalLowSNR",
			Q[(\[Lambda] - M (a + b \[Gamma]0)) / Sqrt[2M a]],
		True,
			Undefined
	];

	Which[
		diversityType == "None" || diversityType == "SC",
			g[1, 1],
		diversityType == "MRC" || diversityType == "EGC",
			g[1, n],
		diversityType == "SEC" && ListQ[\[Gamma]] && (Length[\[Gamma]] == 2),
			(* No \[Gamma]0 here - this is a special case *)
			If[\[Gamma][[1]] >= \[Gamma]t,
				AWGNProbabilityOfDetection[M, \[Gamma][[1]], \[Lambda], DiversityType->"None", Method->OptionValue[Method]],
				AWGNProbabilityOfDetection[M, \[Gamma][[2]], \[Lambda], DiversityType->"None", Method->OptionValue[Method]]
			],
		diversityType == "SLC",
			g[n, n],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - AWGNProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], DiversityType->"None", Method->OptionValue[Method]],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - AWGNProbabilityOfDetection[M, \[Gamma]0, \[Lambda], DiversityType->"None", Method->OptionValue[Method]])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Threshold*)


Options[\[Lambda]] = {Method->"Exact", DiversityType->OptionValue[AWGNProbabilityOfFalseAlarm,DiversityType]};
\[Lambda]::usage="\[Lambda][M, Pf] calculates a threshold suitable for use in the calculation of the decision probabilities for a single energy detector.
\[Lambda][M, Pf, n] calculates a threshold suitable for use in the calculation of the decision probabilities for energy detection with diversity reception.\n\n"<>MethodHelp[\[Lambda], {"\"Exact\"", "\"Approximate\""}]<>"\n\n"<>DiversityTypeHelp[\[Lambda]];
\[Lambda][M_,Pf_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[\[Lambda]][[All,1]]],Options[target][[All,1]]];
	\[Lambda][M,Pf,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[\[Lambda]]]]
]
\[Lambda][M_,Pf_,n_,OptionsPattern[]]:=Module[{method, mn, diversityType = OptionValue[DiversityType], \[Gamma]t, g},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	g[a_,b_] := Which[
		(* Check to see if method is "Exact" or starts with "Exact" *)
		method == "Exact" || StringTake[method, 5] == "Exact",
			2 InverseGammaRegularized[M a / 2, b],
		(* Keep NGaussian here for legacy purposes *)
		method == "Approximate" || StringTake[method, 11] == "Approximate" || method == "NGaussian",
			Sqrt[2M a] InverseQ[b] + M a,
		True,
			Undefined
	];

	Which[
		diversityType == "None" || diversityType == "MRC" || diversityType == "EGC" || diversityType == "SC" || diversityType == "SEC",
			g[1, Pf],
		diversityType == "SLC",
			g[n, Pf],
		diversityType == "SLS",
			g[1, 1 - (1 - Pf)^(1 / n)],
		True,
			Undefined
	]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[AWGNSampleComplexity] = {Method->OptionValue[SampleComplexity,Method], DiversityType->OptionValue[SampleComplexity,DiversityType]}
AWGNSampleComplexity::usage="AWGNSampleComplexity[\[Gamma], Pf, Pd, n] calculates the number of samples required for the specified decision probabilities and signal to noise ratio in an AWGN channel.\n\n"<>MethodHelp[AWGNSampleComplexity, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\""}]<>"\n\n"<>DiversityTypeHelp[AWGNSampleComplexity];
AWGNSampleComplexity[\[Gamma]_,Pf_,Pd_,n_:1,OptionsPattern[]] := Module[{method, mn, diversityType = OptionValue[DiversityType], \[Gamma]t, f, g, M, RelevantOptions},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AWGNSampleComplexity][[All,1]]],Options[target][[All,1]]];

	g[a_,b_] := Which[
		method == "ExactNumerical",
			f[x_?NumericQ] := AWGNProbabilityOfDetection[x, \[Gamma], \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], n, RelevantOptions[AWGNProbabilityOfDetection]];
			M/.FindRoot[f[M] == Pd, {M, AWGNSampleComplexity[\[Gamma], Pf, Pd, n, DiversityType->diversityType, Method->"ApproximateNumerical"], 1, \[Infinity]}],
		method == "ApproximateNumerical",
			2 ((Sqrt[a] InverseQ[Pf] - Sqrt[a + 2 b \[Gamma]] InverseQ[Pd]) / (b \[Gamma]))^2,
		method == "ApproximateNumericalLowSNR" || method == "ApproximateAsymptotic",
			2a ((InverseQ[Pf] - InverseQ[Pd]) / (b \[Gamma]))^2,
		True,
			Undefined
	];

	Quiet[Which[
		diversityType == "None",
			g[1, 1],
		diversityType == "MRC" || diversityType == "EGC",
			g[1, n],
		diversityType == "SLC",
			g[n, n],
		True,
			Undefined
	]]
]


End[]


EndPackage[];

(* ::Package:: *)

(* ::Title:: *)
(*Rice channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in Rice channels.*)
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
(*Version 1.32: Amalgamated algorithms and methods, removed LowSNR option completely.*)
(*Version 1.31: Created NumericalLowSNR method instead of using LowSNR option.*)
(*Version 1.30: Major updates to RiceSampleComplexity, added diversity support, exact and approximate methods.*)
(*Version 1.25: Minor bug fixes for limit functions.*)
(*Version 1.24: Added error bound functions and renamed SmallK method to IntegerKN.*)
(*Version 1.23: Added SEC and limited SC (numerical methods only) support.*)
(*Version 1.22: Finished MRC and SLC implementations.*)
(*Version 1.21: Finished SLS implementations.*)
(*Version 1.2: Moved timing functions to RiceProbabilityOfDetection and added RiceLimit function for public access to truncation points. More minor bug fixes.*)
(*Version 1.1: Major clean up of code, added approximations and numerical methods for no diversity Rice channels.*)
(*Version 1.02: Moved database logging functions to the Network package.*)
(*Version 1.01: Added sample complexity function.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Rice`"]; 


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


RicePDF;


(* ::Subsection:: *)
(*Probabiity of detection*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


RiceProbabilityOfDetection;


RiceLimit;


(* ::Subsubsection::Closed:: *)
(*Annamalai's method*)


AnnamalaiRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Integer-KN method*)


IntegerKNRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Small Pf method*)


SmallPfRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Asymptotic method*)


AsymptoticRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate)*)


NGaussianRiceProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate, low SNR)*)


NGaussianLowSNRRiceProbabilityOfDetection;


(* ::Subsection:: *)
(*Detection probability error bounds*)


(* ::Subsubsection::Closed:: *)
(*Low SNR assumption error*)


LowSNRAssumptionErrorRice;


(* ::Subsubsection::Closed:: *)
(*Asymptotic approximation error*)


AsymptoticErrorRice;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


RiceSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<Network`;
<<AWGN`;
<<ErfApprox`;
<<QFunction`;
<<InverseMarcumQ`;
<<Extras`;


(* ::Subsection::Closed:: *)
(*Help generation*)


GenerateMethodHelp[fName_, methodName_] := ToString[fName] <> "[M, \[Gamma], \[Lambda], K] calculates the probability of detection for a single energy detector operating in a Rice fading channel using the " <> methodName <> " method.
" <> ToString[fName] <> "[M, \[Gamma], \[Lambda], K, n] calculates the probability of detection for energy detection with diversity reception in a Rice fading channel using the " <> methodName <> " method."<>"\n\n"<>DiversityTypeHelp[fName]<>"\n\n"<>TimingHelp[fName];


GenerateTruncationHelp[fName_,methodName_] := ToString[fName] <> "[M, \[Gamma], \[Lambda], K] calculates the truncation point for use in the " <> methodName <> " method for a single energy detector operating on a Rice channel.
" <> ToString[fName] <> "[M, \[Gamma], \[Lambda], K, n] calculates the truncation point for use in the " <> methodName <> " method for energy detection with diversity reception in a Rice channel."<>"\n\n"<>DiversityTypeHelp[fName]<>"\n\n"<>ToleranceHelp[fName];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


Options[RicePDF] = {Method->"Exact", DiversityType->OptionValue[ProbabilityOfDetection,DiversityType]};
RicePDF::usage="RicePDF[\[Gamma], K, x] evaluates the probability density function of the instantaneous signal to noise ratio at a single energy detector operating on a Rice fading channel at x.
RicePDF[\[Gamma], K, x, n] evaluates the probability density function of the average instantaneous signal to noise ratio for energy detection with diversity reception in a Rice fading channel.\n\n"<>MethodHelp[RicePDF, {"\"Exact\"", "\"Approximate\""}]<>"\n\n"<>DiversityTypeHelp[RicePDF];
RicePDF[\[Gamma]_,K_,x_,OptionsPattern[]]:=Module[{n = 1}, RicePDF[\[Gamma], K, x, n, Method->OptionValue[Method], DiversityType->"None"]]
RicePDF[\[Gamma]_,K_,x_,n_,OptionsPattern[]]:=Module[{method, mn, diversityType = OptionValue[DiversityType], \[Gamma]t, \[Gamma]0, g},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		method == "Exact" || StringTake[method, 5] == "Exact",
			Which[
				diversityType == "None",
					((2 (K + 1)) / \[Gamma]0) PDF[NoncentralChiSquareDistribution[2, 2K], (2 (K + 1) x) / \[Gamma]0],
				diversityType == "MRC",
					((2 (K + 1)) / \[Gamma]0) PDF[NoncentralChiSquareDistribution[2n, 2K n], (2 (K + 1) x) / \[Gamma]0],
				diversityType == "EGC",
					Undefined,
				diversityType == "SC",
					(1 - MarcumQ[1, Sqrt[2K],Sqrt[2(K + 1) x / \[Gamma]0]])^(n - 1) n ((2 (K + 1)) / \[Gamma]0) PDF[NoncentralChiSquareDistribution[2, 2K], (2 (K + 1) x) / \[Gamma]0],
				diversityType == "SEC",
					Piecewise[{{(1 - MarcumQ[1, Sqrt[2K], Sqrt[2 (K + 1) \[Gamma]t / \[Gamma]0]])^(n - 1) RicePDF[\[Gamma]0, K, x, Method->method], x < \[Gamma]t}, {Sum[(1 - MarcumQ[1, Sqrt[2K], Sqrt[2 (K + 1) \[Gamma]t / \[Gamma]0]])^j, {j, 0, n - 1}] RicePDF[\[Gamma]0, K, x, Method->method], x >= \[Gamma]t}}],
				diversityType == "SLC",
					((2 (K + 1)) / \[Gamma]0) PDF[NoncentralChiSquareDistribution[2n, 2K n], (2 (K + 1) x) / \[Gamma]0],
				diversityType == "SLS",
					Which[
						ListQ[\[Gamma]0],
							Product[RicePDF[\[Gamma]0[[i]], K, x, DiversityType->"None", Method->"Exact"], {i, n}],
						!ListQ[\[Gamma]0],
							RicePDF[\[Gamma]0, K, x, DiversityType->"None", Method->"Exact"]^n,
						True,
							Undefined
					],
				True,
					Undefined
			],
		method == "Approximate" || StringTake[method, 11] == "Approximate",
			Which[
				diversityType == "None",
					PDF[NormalDistribution[\[Gamma]0, Sqrt[2K + 1] (\[Gamma]0 / (K + 1))], x],
				diversityType == "MRC",
					PDF[NormalDistribution[n \[Gamma]0, Sqrt[(1 + 2K) n] (\[Gamma]0 / (K + 1))], x],
				diversityType == "EGC",
					Undefined,
				diversityType == "SC",
					(1 - MarcumQ[1, Sqrt[2K],Sqrt[2(K + 1) x / \[Gamma]0]])^(n - 1) n PDF[NormalDistribution[\[Gamma]0, Sqrt[2K + 1] (\[Gamma]0 / (K + 1))], x],
				diversityType == "SEC",
					Undefined,
				diversityType == "SLC",
					PDF[NormalDistribution[n \[Gamma]0, Sqrt[(1 + 2K) n] (\[Gamma]0 / (K + 1))], x],
				diversityType == "SLS",
					Which[
						ListQ[\[Gamma]0],
							Product[RicePDF[\[Gamma]0[[i]], K, x, DiversityType->"None", Method->"Approximate"], {i, n}],
						!ListQ[\[Gamma]0],
							RicePDF[\[Gamma]0, K, x, DiversityType->"None", Method->"Approximate"]^n,
						True,
							Undefined
					],
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsection:: *)
(*Probabiity of detection*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


Options[RiceProbabilityOfDetection]={Method->OptionValue[ProbabilityOfDetection,Method],DiversityType->OptionValue[ProbabilityOfDetection,DiversityType],Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
RiceProbabilityOfDetection::usage="RiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K] calculates the probability of detection for a single energy detector operating on a Rice fading channel.
RiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n] calculates the probability of detection for energy detection with diversity reception in a Rice fading channel.\n\n"<>MethodHelp[RiceProbabilityOfDetection, {"\"ExactNumerical\"", "\"ExactAnnamalai\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateIntegerKN\"", "\"ApproximateAsymptotic\""}]<>"\n\n"<>DiversityTypeHelp[RiceProbabilityOfDetection]<>"\n\n"<>TimingHelp[RiceProbabilityOfDetection];
RiceProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,K_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	RiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, RelevantOptions[RiceProbabilityOfDetection]]
]
RiceProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,K_,n_,OptionsPattern[]]:=Module[{RelevantOptions, method, mn, limit, f, totaltime = 0, iterations = 0, time, result},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	limit = RiceLimit[M, \[Gamma], \[Lambda], K, n, Method->method, RelevantOptions[RiceLimit]];

	f := Which[
		(* Catch extreme values - they can cause errors *)
		\[Lambda] == -\[Infinity] || M == \[Infinity] || \[Gamma] == \[Infinity] || n == \[Infinity],
			1,
		method == "ExactNumerical",
			NumericalRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, RelevantOptions[NumericalRiceProbabilityOfDetection]],
		method == "ExactAnnamalai",
			AnnamalaiRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, Limit->limit, RelevantOptions[AnnamalaiRiceProbabilityOfDetection]],
		method == "ApproximateNumerical",
			NGaussianRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, RelevantOptions[NGaussianRiceProbabilityOfDetection]],
		method == "ApproximateIntegerKN",
			IntegerKNRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, Limit->limit, RelevantOptions[IntegerKNRiceProbabilityOfDetection]],
		method == "ApproximateSmallPf",
			SmallPfRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, Limit->limit, RelevantOptions[SmallPfRiceProbabilityOfDetection]],
		method == "ApproximateAsymptotic",
			AsymptoticRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, RelevantOptions[AsymptoticRiceProbabilityOfDetection]],
		method == "ApproximateNumericalLowSNR",
			NGaussianLowSNRRiceProbabilityOfDetection[M, \[Gamma], \[Lambda], K, n, RelevantOptions[NGaussianLowSNRRiceProbabilityOfDetection]],
		True,
			Undefined
	];

	If[OptionValue[Timed],
		(* Evaluate result until MaxTime seconds of CPU time have been used or MaxIterations have been performed, whichever comes first *)
		While[totaltime < OptionValue[MaxTime] && iterations < OptionValue[MaxIterations],
			ClearSystemCache[];
			{time, result} = TimeConstrained[Timing[f],OptionValue[MaxTime],{OptionValue[MaxTime],Null}];
			totaltime += time;
			iterations++;
		];
		{result,totaltime/iterations},
		f
	]
]


Options[RiceLimit] = {Method->"ExactAnnamalai", DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType], Tolerance->10^-6};
RiceLimit::usage = "RiceLimit[M, \[Gamma], \[Lambda], K] calculates the truncation point for use in the specified method for a single energy detector operating on a Rice channel.
RiceLimit[M, \[Gamma], \[Lambda], K, n] calculates the truncation point for use in the specified method for energy detection with diversity reception in a Rice channel.\n\n" <> MethodHelp[RiceLimit, {"\"ExactNumerical\"", "\"ExactAnnamalai\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateIntegerKN\"", "\"ApproximateAsymptotic\""}] <> "\n\n" <> DiversityTypeHelp[RiceLimit];
RiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{n = 1}, RiceLimit[M, \[Gamma], \[Lambda], K, n, DiversityType->"None", Method->OptionValue[Method], Tolerance->OptionValue[Tolerance]]]
RiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, \[Gamma]t, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType], method, mn},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RiceLimit][[All,1]]],Options[target][[All,1]]];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	Which[
		(* Catch extreme values - they can cause errors *)
		\[Lambda] == -\[Infinity] || M == \[Infinity] || \[Gamma] == \[Infinity] || n == \[Infinity],
			Undefined,
		method == "ExactAnnamalai",
			AnnamalaiRiceLimit[M, \[Gamma], \[Lambda], K, n, RelevantOptions[AnnamalaiRiceLimit]],
		method == "ApproximateIntegerKN",
			IntegerKNRiceLimit[M, \[Gamma], \[Lambda], K, n, RelevantOptions[IntegerKNRiceLimit]],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Annamalai's method*)


Options[AnnamalaiRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType], Limit->Null};
AnnamalaiRiceProbabilityOfDetection::usage = GenerateMethodHelp[AnnamalaiRiceProbabilityOfDetection, "Annamalai"];
AnnamalaiRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AnnamalaiRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	AnnamalaiRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[AnnamalaiRiceProbabilityOfDetection]]
]
AnnamalaiRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], f, \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = RiceLimit[M, \[Gamma], \[Lambda], K, n, Method->"ExactAnnamalai", DiversityType->diversityType]];

	Which[
		diversityType == "None",
			1 - (2 (1 + K) / (2 + 2K + M \[Gamma]0)) Exp[-K M \[Gamma]0 / (2 + 2K + M \[Gamma]0)] (1 - GammaRegularized[M/2, \[Lambda] / 2]) - Total[Table[(k!) (2 (1 + K) / (2 + 2K + M \[Gamma]0)) (((M \[Gamma]0) / (2 + 2K + M \[Gamma]0))^k) Exp[-K M \[Gamma]0 / (2 + 2K + M \[Gamma]0)] (1 / k! + Total[Table[((2 K (1 + K) / (2 + 2K + M \[Gamma]0))^i) / (i! i! (k - i)!),{i , 1, k}]]) (1 - GammaRegularized[M / 2 + k, \[Lambda] / 2]),{k, 1, limit}]],
		diversityType == "MRC",
			1 - ((1 - GammaRegularized[M / 2, \[Lambda] / 2]) (2 (1 + K) / (2 + 2K + M \[Gamma]0))^n Exp[-K M n \[Gamma]0 / (2 + 2K + M \[Gamma]0)]) - Total[Table[(1 - GammaRegularized[M / 2 + k, \[Lambda] / 2]) ((n + k - 1)!) (((M \[Gamma]0) / (2 + 2K + M \[Gamma]0))^k) ((2 (1 + K) / (2 + 2K + M \[Gamma]0))^n) Exp[-K M n \[Gamma]0 / (2 + 2K + M \[Gamma]0)] Total[Table[(1 / (i! (n + i - 1)! (k - i)!)) (2K n (K + 1) / (2 + 2K + M \[Gamma]))^i,{i, 0, k}]], {k, 1, limit}]],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			1 - ((2 (1 + K) / (2 + 2K + M \[Gamma]0))^n) Exp[-K n M \[Gamma]0 / (2 + 2K + M \[Gamma]0)] (1 - GammaRegularized[M n/2, \[Lambda] / 2]) - Total[Table[((n + k - 1)!) ((2 (1 + K) / (2 + 2K + M \[Gamma]0))^n) (((M \[Gamma]0) / (2 + 2K + M \[Gamma]0))^k) Exp[-K n M \[Gamma]0 / (2 + 2K + M \[Gamma]0)] (1 / ((n - 1)! k!) + Total[Table[((2 K n (1 + K) / (2 + 2K + M \[Gamma]0))^i) / (i! (n + i - 1)! (k - i)!),{i , 1, k}]]) (1 - GammaRegularized[M n / 2 + k, \[Lambda] / 2]),{k, 1, limit}]],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - AnnamalaiRiceProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], K, Limit->limit[[i]], DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - AnnamalaiRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, Limit->limit, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


Options[AnnamalaiRiceLimit] = {DiversityType->OptionValue[RiceLimit,DiversityType], Tolerance->OptionValue[RiceLimit,Tolerance]};
AnnamalaiRiceLimit::usage = GenerateTruncationHelp[AnnamalaiRiceLimit, "Annamalai"];
AnnamalaiRiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{n = 1},AnnamalaiRiceLimit[M,\[Gamma],\[Lambda],K,n,DiversityType->"None",Tolerance->OptionValue[Tolerance]]]
AnnamalaiRiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	Quiet[
		Which[
			diversityType == "None",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol];
				j/.FindRoot[1 - GammaRegularized[(M / 2) + j + 1, \[Lambda] / 2] == tol,{j, j0, 1, \[Infinity]}],
			diversityType == "MRC",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol];
				j/.FindRoot[1 - GammaRegularized[(M / 2) + j + 1, \[Lambda] / 2] == tol,{j, j0, 1, \[Infinity]}],
			diversityType == "EGC",
				Undefined,
			diversityType == "SC",
				Undefined,
			diversityType == "SEC",
				Undefined,
			diversityType == "SLC",
				j0 = (\[Lambda] / 2) - (M n / 2) - Sqrt[M n / 2] InverseQ[1 - tol];
				j/.FindRoot[1 - GammaRegularized[(M / 2) n + j + 1, \[Lambda] / 2] == tol,{j, j0, 1, \[Infinity]}],
			diversityType == "SLS",
				Which[
					ListQ[\[Gamma]],
						Table[AnnamalaiRiceLimit[M,\[Gamma][[i]],\[Lambda],K,DiversityType->"None",Tolerance->tol], {i, Length[\[Gamma]]}],
					!ListQ[\[Gamma]],
						AnnamalaiRiceLimit[M,\[Gamma],\[Lambda],K,DiversityType->"None",Tolerance->tol],
					True,
						Undefined
				],
			True,
				Undefined
		]//N//Ceiling
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method*)


Options[NumericalRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType]};
NumericalRiceProbabilityOfDetection::usage = GenerateMethodHelp[NumericalRiceProbabilityOfDetection, "Numerical"];
NumericalRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NumericalRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NumericalRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[NumericalRiceProbabilityOfDetection]]
]
NumericalRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{lim, f, g, \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType], \[ScriptCapitalD], x},
	\[ScriptCapitalD] = NumericalRiceDistribution[M, \[Gamma], K, n, diversityType];

	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None" || diversityType == "MRC" || diversityType == "EGC" || diversityType == "SC" || diversityType == "SEC" || diversityType == "SLC",
			1 - CDF[\[ScriptCapitalD], \[Lambda]],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - NumericalRiceProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], K, DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NumericalRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]//N
]


NumericalRiceDistribution[M_?NumericQ, \[Gamma]_, K_?NumericQ, n_?IntegerQ, diversityType0_] := NumericalRiceDistribution[M, \[Gamma], K, n, diversityType0] = Module[{\[Gamma]t, \[Gamma]0, g, numberOfPoints = 1000000, x, diversityType, y, \[ScriptCapitalD]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType0];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	\[ScriptCapitalD] = ProbabilityDistribution[RicePDF[\[Gamma]0, K, y, n, DiversityType->diversityType0, Method->"Exact"], {y, 0, \[Infinity]}];

	g[a_] := SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M a, M x], x \[Distributed] \[ScriptCapitalD]], numberOfPoints]];

	Which[
		diversityType == "None",
			g[1],
		diversityType == "MRC",
			g[1],
		diversityType == "EGC",
			x = Table[Unique[], {n}];
			SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M, M Sum[Sqrt[x[[i]]], {i, n}]^2/n], Table[x[[i]] \[Distributed] ProbabilityDistribution[RicePDF[\[Gamma]0, K, y, n, DiversityType->"None", Method->"Exact"], {y, 0, \[Infinity]}], {i, n}]], numberOfPoints]],
		diversityType == "SC",
			g[1],
		diversityType == "SEC",
			g[1],
		diversityType == "SLC",
			g[n],
		diversityType == "SLS",
			(* We don't need the SLS PDF here, the master function will use the no diversity PDF instead *)
			Null,
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Integer-KN method*)


Options[IntegerKNRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType], Limit->Null};
IntegerKNRiceProbabilityOfDetection::usage = GenerateMethodHelp[IntegerKNRiceProbabilityOfDetection, "IntegerKN"];
IntegerKNRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[IntegerKNRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	IntegerKNRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[IntegerKNRiceProbabilityOfDetection]]
];
IntegerKNRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Nu], \[Gamma]0, \[Gamma]t, f, J, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = RiceLimit[M, \[Gamma], \[Lambda], K, n, Method->"ApproximateIntegerKN", DiversityType->diversityType]];

	(* Precision is set to 100 here - required for FaddeevaDerivative to be stable for large orders *)
	J[k_, a_, b_, c_] := (Exp[-a^2] / 2) ((I c / (2 b))^k / k!) FaddeevaDerivative[k, N[-I (a + (c / (2b))), 100]] // Re;

	Which[
		diversityType == "None",
			With[{a = (\[Lambda] - M) / (2 Sqrt[M]), b = - Sqrt[M] / 2, c = (K + 1) / \[Gamma]0},
				AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType] + Exp[-K] (J[0, a, b, c] + Total[Table[K^k / k! Total[Table[J[p, a, b, c], {p, 0, k}]], {k, 1, limit}]])
			],
		diversityType == "MRC",
			With[{a = (\[Lambda] - M) / (2 Sqrt[M]), b = - Sqrt[M] / 2, c = (K + 1) / \[Gamma]0},
				AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType] + Exp[-K n] (Total[Table[J[p, a, b, c], {p, 0, n - 1}]] + Total[Table[(K n)^k / k! Total[Table[J[p, a, b, c], {p, 0, n + k - 1}]], {k, 1, limit}]])
			],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			With[{a = (\[Lambda] - M) / (2 Sqrt[M]), b = - Sqrt[M] / 2, c = (K + 1) / \[Gamma]0},
				(1 - MarcumQ[1, Sqrt[2K], Sqrt[2 c \[Gamma]t]])^(n - 1) IntegerKNRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, Limit->limit[[1]]] + Total[Table[(1 - MarcumQ[1, Sqrt[2K], Sqrt[2 c \[Gamma]t]])^j, {j, 0, n - 2}]] (((1 / 2) (1 - Erf[a + b \[Gamma]t]) MarcumQ[1, Sqrt[2K], Sqrt[2 c \[Gamma]t]]) + (Exp[-K] (GammaRegularized[1, c \[Gamma]t] J[0, a + b \[Gamma]t, b, c] + Total[Table[K^k / k! Total[Table[GammaRegularized[k + 1 - p, c \[Gamma]t] J[p, a + b \[Gamma]t, b, c], {p, 0, k}]], {k, 1, limit[[2]]}]])))
			],
		diversityType == "SLC",
			With[{a = (\[Lambda] - M n) / (2 Sqrt[M n]), b = - Sqrt[M / n] / 2, c = (K + 1) / \[Gamma]0},
				AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType] + Exp[-K n] (Total[Table[J[p, a, b, c], {p, 0, n - 1}]] + Total[Table[(K n)^k / k! Total[Table[J[p, a, b, c], {p, 0, n + k - 1}]], {k, 1, limit}]])
			],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - IntegerKNRiceProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], K, Limit->limit[[i]], DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - IntegerKNRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, Limit->limit, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


Options[IntegerKNRiceLimit] = {DiversityType->OptionValue[RiceLimit,DiversityType], Tolerance->OptionValue[RiceLimit,Tolerance]};
IntegerKNRiceLimit::usage = GenerateTruncationHelp[IntegerKNRiceLimit, "IntegerKN"];
IntegerKNRiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{n = 1},IntegerKNRiceLimit[M,\[Gamma],\[Lambda],K,n,DiversityType->"None",Tolerance->OptionValue[Tolerance]]]
IntegerKNRiceLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			j0 = 1;
			j/.FindRoot[(1 - GammaRegularized[j + 1, K]) (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType]) == tol, {j, j0, 1, \[Infinity]}],
		diversityType == "MRC",
			j0 = 1;
			j/.FindRoot[(1 - GammaRegularized[j + 1, K n]) (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType]) == tol, {j, j0, 1, \[Infinity]}],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			j0 = 1;
			With[{a = (\[Lambda] - M) / (2 Sqrt[M]), b = - Sqrt[M] / 2, c = (K + 1) / \[Gamma]0},
				{IntegerKNRiceLimit[M, \[Gamma], \[Lambda], K, Tolerance->tol], j/.FindRoot[(1 - MarcumQ[1, Sqrt[2K], Sqrt[2 c \[Gamma]t]])^j (1 - GammaRegularized[j + 1, K]) (1 / 2) (1 + Erf[a + b \[Gamma]t])  == tol, {j, j0, 1, \[Infinity]}]}
			],
		diversityType == "SLC",
			j0 = 1;
			j/.FindRoot[(1 - GammaRegularized[j + 1, K n]) (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, DiversityType->diversityType]) == tol, {j, j0, 1, \[Infinity]}],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]],
					Table[IntegerKNRiceLimit[M,\[Gamma][[i]],\[Lambda],K,DiversityType->"None",Tolerance->tol], {i, Length[\[Gamma]]}],
				!ListQ[\[Gamma]],
					IntegerKNRiceLimit[M,\[Gamma],\[Lambda],K,DiversityType->"None",Tolerance->tol],
				True,
					Undefined
			],
		True,
			Undefined
	]//Ceiling//Quiet
]


(* ::Subsubsection:: *)
(*Small Pf method*)


Options[SmallPfRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType], Limit->Null};
SmallPfRiceProbabilityOfDetection::usage = GenerateMethodHelp[SmallPfRiceProbabilityOfDetection, "SmallPf"];
SmallPfRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[SmallPfRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	SmallPfRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[SmallPfRiceProbabilityOfDetection]]
]
SmallPfRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Nu], \[Gamma]0, \[Gamma]t, f, g, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = RiceLimit[M, \[Gamma], \[Lambda], K, n, Method->"ApproximateSmallPf", DiversityType->diversityType]];

	(*g[a_, b_, c_] := AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType] + (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType]) Exp[-K] Total[Table[K^k / k! GammaRegularized[k + 1, - a c / b], {k, 0, limit}]];*)
	g[a_, b_, c_] := AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType] + (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType]) MarcumQ[1, Sqrt[2K], Sqrt[-2 a c / b]];

	Which[
		diversityType == "None",
			g[(\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, (K + 1) / \[Gamma]0],
		diversityType == "MRC",
			g[(\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, (K + 1) / \[Gamma]0],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			Undefined,
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - SmallPfRiceProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], K, Limit->limit[[i]], DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - SmallPfRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, Limit->limit, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Asymptotic method*)


Options[AsymptoticRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType]};
AsymptoticRiceProbabilityOfDetection::usage = GenerateMethodHelp[AsymptoticRiceProbabilityOfDetection, "Asymptotic"];
AsymptoticRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AsymptoticRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	AsymptoticRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[AsymptoticRiceProbabilityOfDetection]]
]
AsymptoticRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{lim, f, g, \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[a_, b_, \[Mu]_, \[Sigma]_] := Q[(a + b \[Mu]) / Sqrt[b^2 \[Sigma]^2 + 1]];

	Which[
		diversityType == "None",
			g[(\[Lambda] - M) / Sqrt[2 M], - Sqrt[M / 2], \[Gamma]0, Sqrt[2K + 1] (\[Gamma]0 / (K + 1))],
		diversityType == "MRC",
			g[(\[Lambda] - M) / Sqrt[2 M], - Sqrt[M / 2], n \[Gamma]0, Sqrt[n (2K + 1)] (\[Gamma]0 / (K + 1))],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			g[(\[Lambda] - M n) / Sqrt[2 M n], - Sqrt[M / (2 n)], n \[Gamma]0, Sqrt[n (2K + 1)] (\[Gamma]0 / (K + 1))],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - AsymptoticRiceProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], K, DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - AsymptoticRiceProbabilityOfDetection[M, \[Gamma]0, \[Lambda], K, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate)*)


Options[NGaussianRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType]};
NGaussianRiceProbabilityOfDetection::usage = GenerateMethodHelp[NGaussianRiceProbabilityOfDetection, "NGaussian"];
NGaussianRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NGaussianRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[NGaussianRiceProbabilityOfDetection]]
]
NGaussianRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, diversityType = OptionValue[DiversityType], method = "ApproximateNumerical"},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "MRC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "EGC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SEC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->{diversityType, \[Gamma]t}],{x,0,\[Infinity]}],
		diversityType == "SLC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - NGaussianRiceProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],K,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NGaussianRiceProbabilityOfDetection[M,\[Gamma]0,\[Lambda],K,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate, low SNR)*)


Options[NGaussianLowSNRRiceProbabilityOfDetection] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType]};
NGaussianLowSNRRiceProbabilityOfDetection::usage = GenerateMethodHelp[NGaussianLowSNRRiceProbabilityOfDetection, "NGaussian"];
NGaussianLowSNRRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianLowSNRRiceProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NGaussianLowSNRRiceProbabilityOfDetection[M,\[Gamma],\[Lambda],K,n,RelevantOptions[NGaussianLowSNRRiceProbabilityOfDetection]]
]
NGaussianLowSNRRiceProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, diversityType = OptionValue[DiversityType], method = "ApproximateNumericalLowSNR"},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "MRC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "EGC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SEC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->{diversityType, \[Gamma]t}],{x,0,\[Infinity]}],
		diversityType == "SLC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]RicePDF[\[Gamma]0,K,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - NGaussianRiceProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],K,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NGaussianRiceProbabilityOfDetection[M,\[Gamma]0,\[Lambda],K,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsection:: *)
(*Detection probability error bounds*)


(* ::Subsubsection::Closed:: *)
(*Low SNR assumption error*)


Options[LowSNRAssumptionErrorRice] = {DiversityType->OptionValue[RiceProbabilityOfDetection,DiversityType]};
LowSNRAssumptionErrorRice::usage="LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K, n] calculates the upper bound for the low SNR approximation error.\n\n"<>DiversityTypeHelp[LowSNRAssumptionErrorRice];
LowSNRAssumptionErrorRice[M_,\[Gamma]_,\[Lambda]_,K_] := LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K] = Module[{n = 1},
	LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K, n, DiversityType->"None"]
]
LowSNRAssumptionErrorRice[M_,\[Gamma]_,\[Lambda]_,K_,n_,OptionsPattern[]] := Module[{diversityType = OptionValue[DiversityType], \[Gamma]t, g, \[Epsilon], \[Epsilon]max = 10},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	g[\[Epsilon]_?NumericQ] := Which[
		diversityType == "None" || diversityType == "SLC",
			Abs[AWGNProbabilityOfDetection[M, \[Epsilon], \[Lambda], n, DiversityType->diversityType, Method->ApproximateNumerical] - AWGNProbabilityOfDetection[M, \[Epsilon], \[Lambda], n, DiversityType->diversityType, Method->ApproximateNumericalLowSNR]],
		diversityType == "MRC" || diversityType == "EGC" || diversityType == "SC" || diversityType == "SEC",
			LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K],
		True,
			Undefined
	];

	If[diversityType == "SLS",
		(* This bound is VERY loose *)
		With[{Pm = (1 - ProbabilityOfDetection[M, \[Gamma], \[Lambda], ChannelType -> {"Rice", K}, Method -> "ApproximateNumerical", DiversityType -> "None"])},
			Min[Abs[{Pm^n - (Pm - LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K])^n, Pm^n - (Pm + LowSNRAssumptionErrorRice[M, \[Gamma], \[Lambda], K])^n}]]
		],
		NMaximize[{g[\[Epsilon]], 0 <= \[Epsilon] <= \[Epsilon]max}, {\[Epsilon], 0, \[Epsilon]max}][[1]]
	]
]


(* ::Subsubsection::Closed:: *)
(*Asymptotic approximation error*)


AsymptoticErrorRice::usage="AsymptoticErrorRice[Pf, K, n] gives the upper bound for the error of the asymptotic method for the specified parameters.";
AsymptoticErrorRice[Pf_,K_,n_] := AsymptoticErrorRice[Pf, K, n] = Module[{f, z},
	f[z_] := (Pf / 2) Erfc[(K + 1) n / Sqrt[2 (2K + 1) n]] + (1 - Pf) (MarcumQ[n, Sqrt[2K n], Sqrt[2 z]] - (1 / 2) Erfc[(z - (K + 1) n)/Sqrt[2 (2K + 1) n]]);

	NMaximize[{Abs[f[z]], z >= 0}, {z, 0, (K + 1) n}][[1]]
]


(* ::Subsection:: *)
(*Sample complexity*)


Options[RiceSampleComplexity] = {DiversityType->OptionValue[SampleComplexity,DiversityType], Method->OptionValue[SampleComplexity,Method]};
RiceSampleComplexity::usage="RiceSampleComplexity[\[Gamma], Pf, Pd, K, n] calculates the number of samples required for the specified decision probabilities and signal to noise ratio in a Rice channel.\n\n"<>MethodHelp[RiceSampleComplexity, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateAsymptotic\""}]<>"\n\n"<>DiversityTypeHelp[RiceSampleComplexity];
RiceSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,K_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RiceSampleComplexity][[All,1]]],Options[target][[All,1]]];
	RiceSampleComplexity[\[Gamma], Pf, Pd, K, n, RelevantOptions[RiceSampleComplexity]]
]
RiceSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,K_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, RelevantOptions, diversityType = OptionValue[DiversityType], method, mn, g, M, limit},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RiceSampleComplexity][[All,1]]],Options[target][[All,1]]];

	Quiet[Which[
		method == "ExactNumerical",
			f[x_?NumericQ] := RiceProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], K, n, RelevantOptions[RiceProbabilityOfDetection]/."ExactNumerical"->"ExactAnnamalai"];
			M/.FindRoot[f[M] == Pd, {M, RiceSampleComplexity[\[Gamma]0, Pf, Pd, K, n, Method->"ApproximateNumerical", DiversityType->diversityType], 1, \[Infinity]}],
		method == "ApproximateNumerical",
			f[x_?NumericQ] := RiceProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], K, n, RelevantOptions[RiceProbabilityOfDetection]];
			M/.FindRoot[f[M] == Pd, {M, AWGNSampleComplexity[\[Gamma]0, Pf, Pd, n, RelevantOptions[RiceSampleComplexity]], 1, \[Infinity]}],
		method == "ApproximateNumericalLowSNR",
			f[x_?NumericQ] := RiceProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], K, n, RelevantOptions[RiceProbabilityOfDetection]];
			M/.FindRoot[f[M] == Pd, {M, AWGNSampleComplexity[\[Gamma]0, Pf, Pd, n, RelevantOptions[RiceSampleComplexity]], 1, \[Infinity]}],
		method == "ApproximateSmallPf",
			Which[
				diversityType == "None",
					2 (2 (K + 1) InverseQ[Pf] / (\[Gamma]0 InverseMarcumQ[1, Sqrt[2K], (Pd - Pf) / (1 - Pf)]^2))^2,
				diversityType == "MRC",
					Undefined,
				diversityType == "SLC",
					Undefined,
				True,
					Undefined
			],
		method == "ApproximateAsymptotic",
			Which[
				diversityType == "None",
					(2 (1+K)^2 (-(1+2 K) InverseQ[Pd]^4+InverseQ[Pd]^2 ((1+K)^2+(1+2 K) InverseQ[Pf]^2)+(1+K) ((1+K) InverseQ[Pf]^2+2 \[Sqrt](InverseQ[Pd]^2 InverseQ[Pf]^2 ((1+K)^2-(1+2 K) InverseQ[Pd]^2+(1+2 K) InverseQ[Pf]^2)))))/(\[Gamma]0^2 ((1+K)^2-(1+2 K) InverseQ[Pd]^2)^2),
				diversityType == "MRC",
					(2 (1+K)^2 (-(1+2 K) InverseQ[Pd]^4+InverseQ[Pd]^2 ((1+K)^2 n+(1+2 K) InverseQ[Pf]^2)+(1+K) ((1+K) n InverseQ[Pf]^2+2 \[Sqrt](n InverseQ[Pd]^2 InverseQ[Pf]^2 ((1+K)^2 n-(1+2 K) InverseQ[Pd]^2+(1+2 K) InverseQ[Pf]^2)))))/(n ((1+K)^2 n \[Gamma]0-(\[Gamma]0+2 K \[Gamma]0) InverseQ[Pd]^2)^2),
				diversityType == "SLC",
					(2 (1+K)^2 (-(1+2 K) InverseQ[Pd]^4+InverseQ[Pd]^2 ((1+K)^2 n+(1+2 K) InverseQ[Pf]^2)+(1+K) ((1+K) n InverseQ[Pf]^2+2 \[Sqrt](n InverseQ[Pd]^2 InverseQ[Pf]^2 ((1+K)^2 n-(1+2 K) InverseQ[Pd]^2+(1+2 K) InverseQ[Pf]^2)))))/((1+K)^2 n \[Gamma]0-(\[Gamma]0+2 K \[Gamma]0) InverseQ[Pd]^2)^2,
				True,
					Undefined
			],
		method == "ApproximateAsymptoticUpperBound",
			\[Epsilon] = AsymptoticErrorRice[Pf, K, n];
			If[Pd + \[Epsilon] >= 1, \[Infinity], RiceSampleComplexity[\[Gamma], Pf, Pd + \[Epsilon], K, n, Method->"ApproximateAsymptotic", DiversityType->diversityType]],
		method == "ApproximateAsymptoticLowerBound",
			\[Epsilon] = AsymptoticErrorRice[Pf, K, n];
			If[Pd - \[Epsilon] <= 0, 0, RiceSampleComplexity[\[Gamma], Pf, Pd - \[Epsilon], K, n, Method->"ApproximateAsymptotic", DiversityType->diversityType]],
		True,
			Undefined
	]]
]


End[];


EndPackage[];

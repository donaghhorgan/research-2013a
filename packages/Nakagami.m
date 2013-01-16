(* ::Package:: *)

(* ::Title:: *)
(*Nakagami channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in Nakagami channels.*)
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
(*09/01/2012*)
(*1.71*)


(* ::Subsection:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.71: Improved the computation time of the AnnamalaiNakagamiLimit function.*)
(*Version 1.70: Added no diversity, MRC and EGC methods to the LargeSNRNakagamiProbabilityOfDetection function. Added LargeSNR sample complexity algorithms for no diversity, MRC, EGC and SLC.*)
(*Version 1.63: Added asymptotic error bounds to the NakagamiSampleComplexity function.*)
(*Version 1.62: Amalgamated algorithms and methods.*)
(*Version 1.61: Created NumericalLowSNR method instead of using LowSNR option.*)
(*Version 1.60: Major updates for NakagamiSampleComplexity: added diversity support, exact and approximate methods.*)
(*Version 1.59: Fixed minor bugs in limit functions.*)
(*Version 1.58: Moved timing functions to the main function, and minor bug fixes.*)
(*Version 1.57: Moved FaddeevaDerivative function to the Extras package.*)
(*Version 1.56: Added approximaton error for Nakagami's PDF for EGC diversity.*)
(*Version 1.55: Fixed a bug in the asymptotic error bound function.*)
(*Version 1.54: Added asymptotic error bound functions for EGC and MRC.*)
(*Version 1.53: Added SEC support to IntegerMN, NGaussian and Numerical methods. Herath's method supports SEC with n = 2, i.e. SSC.*)
(*Version 1.52: Added new EGC methods to IntegerMN and Asymptotic algorithms, based on Nakagami's approximation for the sum of iid Nakagami rvs. Dharmawansa's exact PDF-based method is no longer enabled.*)
(*Version 1.51: Improved exact numerical method with smooth kernel distributions and caching.*)
(*Version 1.5: Finished polishing up error bound functions. SLS bound is very loose, but works.*)
(*Version 1.49: Moved error bounds to separate section, added LowSNRAssumptionErrorNakagami from AWGN package.*)
(*Version 1.48: Updated all function help definitions and fixed bug where LargeSNR algorithm was not publicly accessible.*)
(*Version 1.47: Moved ProcessDiversityType and ProcessSNR functions to the Extras package.*)
(*Version 1.46: Finished updating IntegerMN method with FaddeevaDerivative-based algorithms.*)
(*Version 1.45: Updated Asymptotic approximation for no diversity, MRC, SLC and SLS diversity cases.*)
(*Version 1.44: Recoded Numerical, Sun, Herath, Digham and Annamalai's methods to use new ProcessDiversityType and ProcessSNR functions. All methods are much easier to read and understand now.*)
(*Version 1.43: Retitled LargeMN method to Lopez-Benitez method, and recoded for speed. Also added Asymptotic method.*)
(*Version 1.42: Changed IntegerMN method to support faster computation with FaddeevaDerivative function.*)
(*Version 1.41: Added large SNR approximation method for SLC and no diversity types.*)
(*Version 1.4: Added an exact numeric method, minor bug fixes and retitled some methods.*)
(*Version 1.32: Added SC and SSC diversity to the IntegerMN method.*)
(*Version 1.31: Filled out all known diversity types for each exact method.*)
(*Version 1.3: Added MRC, EGC, SLC, SSC, SC and SLS diversity cases to all methods.*)
(*Version 1.2: Split Horgan's approximation into separate IntegerMN and LargeMN functions, so that full functionality can be accessed through the NakagamiProbabilityOfDetection interface.*)
(*Version 1.11: Moved database logging functions to the Network package.*)
(*Version 1.1: Introduced RelevantOptions function and changed function definitions, so that child options are inherited from parents. The Gaussian approximation method is now called Horgan's approximation to avoid confusion with the numerical Gaussian method.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Nakagami`"]; 


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


NakagamiPDF;


(* ::Subsection:: *)
(*Detection probability*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


NakagamiProbabilityOfDetection;


NakagamiLimit;


(* ::Subsubsection::Closed:: *)
(*Annamalai' s method*)


AnnamalaiNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Digham' s method*)


DighamNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Herath' s method*)


HerathNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Sun' s method*)


SunNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical method*)


NumericalNakagamiProbabilityOfDetection;


NumericalDistribution;


(* ::Subsubsection::Closed:: *)
(*Integer mn method*)


IntegerMNNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Asymptotic method*)


AsymptoticNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Large SNR method*)


LargeSNRNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Lopez-Benitez asymptotic method*)


LopezBenitezNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate)*)


NGaussianNakagamiProbabilityOfDetection;


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate, low SNR)*)


NGaussianLowSNRNakagamiProbabilityOfDetection;


(* ::Subsection::Closed:: *)
(*Detection probability error bounds*)


LowSNRAssumptionErrorNakagami;


NakagamiEGCPDFApproximationError;


AsymptoticErrorNakagami;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


NakagamiSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<Network`;
<<AWGN`;
<<ErfApprox`;
<<QFunction`;
<<Extras`;


(* ::Subsection::Closed:: *)
(*Help generation*)


GenerateMethodHelp[fName_, methodName_] := ToString[fName] <> "[M, \[Gamma], \[Lambda], m] calculates the probability of detection for a single energy detector operating in a Nakagami-m fading channel using the " <> methodName <> " method.
" <> ToString[fName] <> "[M, \[Gamma], \[Lambda], m, n] calculates the probability of detection for energy detection with diversity reception in a Nakagami-m fading channel using the " <> methodName <> " method."<>"\n\n"<>DiversityTypeHelp[fName]<>"\n\n"<>TimingHelp[fName];


GenerateTruncationHelp[fName_,methodName_] := ToString[fName] <> "[M, \[Gamma], \[Lambda], m] calculates the truncation point for use in the " <> methodName <> " method for a single energy detector operating on a Nakagami-m channel.
" <> ToString[fName] <> "[M, \[Gamma], \[Lambda], m, n] calculates the truncation point for use in the " <> methodName <> " method for energy detection with diversity reception in a Nakagami-m channel."<>"\n\n"<>DiversityTypeHelp[fName]<>"\n\n"<>ToleranceHelp[fName];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


Options[NakagamiPDF] = {Method->"Exact", DiversityType->OptionValue[ProbabilityOfDetection,DiversityType]};
NakagamiPDF::usage = "NakagamiPDF[\[Gamma], m, x] evaluates the probability density function of the instantaneous signal to noise ratio at a single energy detector operating on a Nakagami-m fading channel at x.
NakagamiPDF[\[Gamma], m, x, n] evaluates the probability density function of the average instantaneous signal to noise ratio for energy detection with diversity reception in a Nakagami-m fading channel.\n\n"<>MethodHelp[NakagamiPDF, {"\"Exact\"", "\"Approximate\""}]<>"\n\n"<>DiversityTypeHelp[NakagamiPDF];
NakagamiPDF[\[Gamma]_,m_,x_,OptionsPattern[]]:=Module[{n = 1},NakagamiPDF[\[Gamma],m,x,n,Method->OptionValue[Method],DiversityType->OptionValue[DiversityType]]]
NakagamiPDF[\[Gamma]_,m_,x_,n_,OptionsPattern[]]:=Module[{method, mn, diversityType = OptionValue[DiversityType], \[Gamma]t, \[Gamma]0, g},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[a_, b_] := Which[
		method == "Exact" || StringTake[method, 5] == "Exact",
			PDF[GammaDistribution[m a, b / m], x],
		method == "Approximate" || StringTake[method, 11] == "Approximate",
			PDF[NormalDistribution[a b, Sqrt[a b^2 / m]], x],
		True,
			Undefined
	];

	Which[
		method == "Exact" || StringTake[method, 5] == "Exact",
			Which[
				diversityType == "None",
					g[1, \[Gamma]0],
				diversityType == "MRC",
					g[n, \[Gamma]0],
				diversityType == "EGC",
					(* Dharmawansa's exact PDF *)
					(*Which[
						n == 2,
							(2 Sqrt[\[Pi]] x^(2m - 1) Exp[-2m x / \[Gamma]0] / (2^(4m - 1))) (Gamma[2m] / (Gamma[m]^2 Gamma[2m + (1 / 2)])) (2m / \[Gamma]0)^(2m) Hypergeometric1F1[2m, 2m + (1 / 2), m x / \[Gamma]0],
						n == 3,
							(4 Sqrt[\[Pi]] Gamma[2m] Exp[-3m x / \[Gamma]0] / (Gamma[m]^3 2^(4m - 1))) Sum[((Gamma[2m+l] Gamma[4m+2l]) / (Gamma[2m+l+1/2] Gamma[6m+2l] Gamma[l+1] 2^l)) x^(3m+l-1) (3m / \[Gamma]0)^(3m+l) HypergeometricPFQ[{2m,4m+2l},{3m+l+1/2,3m+l},3m x / (2\[Gamma]0)],{l,0,\[Infinity]}],
						True,
							Undefined
					],*)
					(* Nakagami's approximate PDF *)					
					g[n, \[Gamma]0 (1 + (m (-1 + n) Gamma[1 / 2 + m]^2) / Gamma[1 + m]^2) / n],
				diversityType == "SC",
					(n / Gamma[m]) Sum[(-1)^(l) Binomial[n - 1, l] Sum[MultinomialCoefficient[l, k, m] (m / \[Gamma]0)^(m + k) x^(m + k - 1) Exp[-m (l + 1) x / \[Gamma]0],{k, 0, l (m - 1)}],{l, 0, n - 1}],
				diversityType == "SEC",
					Piecewise[{{(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0])^(n - 1) g[1, \[Gamma]0], x < \[Gamma]t}, {Sum[(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0])^j, {j, 0, n - 1}] g[1, \[Gamma]0], x >= \[Gamma]t}}],
				diversityType == "SLC",
					g[n, \[Gamma]0],
				diversityType == "SLS",
					Which[
						ListQ[\[Gamma]0],
							Product[NakagamiPDF[\[Gamma]0[[i]], m, x, n, DiversityType->"None"], {i, n}],
						!ListQ[\[Gamma]0],
							g[1, \[Gamma]0]^n,
						True,
							Undefined
					],
				True,
					Undefined
			],
		method == "Approximate" || StringTake[method, 11] == "Approximate",
			Which[
				diversityType == "None",
					g[1, \[Gamma]0],
				diversityType == "MRC",
					g[n, \[Gamma]0],
				diversityType == "EGC",
					(* Nakagami's approximate PDF *)
					PDF[NormalDistribution[(\[Gamma]0 + (m (-1 + n) \[Gamma]0 Gamma[1 / 2 + m]^2) / Gamma[1 + m]^2), (\[Gamma]0 + (m (-1 + n) \[Gamma]0 Gamma[1 / 2 + m]^2)/Gamma[1 + m]^2) Sqrt[m n]], x],
				diversityType == "SC",
					Undefined,
				diversityType == "SEC",
					Undefined,
				diversityType == "SLC",
					g[n, \[Gamma]0],
				diversityType == "SLS",
					Which[
						ListQ[\[Gamma]0],
							Product[NakagamiPDF[\[Gamma]0[[i]], m, x, n, DiversityType->"None"], {i, n}],
						!ListQ[\[Gamma]0],
							g[1]^n,
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
(*Detection probability*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


Options[NakagamiProbabilityOfDetection] = {Method->OptionValue[ProbabilityOfDetection,Method], DiversityType->OptionValue[ProbabilityOfDetection,DiversityType],Timed->OptionValue[ProbabilityOfDetection,Timed],MaxTime->OptionValue[ProbabilityOfDetection,MaxTime],MaxIterations->OptionValue[ProbabilityOfDetection,MaxIterations]};
NakagamiProbabilityOfDetection::usage = "NakagamiProbabilityOfDetection[M, \[Gamma], \[Lambda], m] calculates the probability of detection for a single energy detector operating on a Nakagami-m fading channel.
NakagamiProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n] calculates the probability of detection for energy detection with diversity reception in a Nakagami-m fading channel.\n\n"<>MethodHelp[NakagamiProbabilityOfDetection, {"\"ExactNumerical\"", "\"ExactAnnamalai\"", "\"ExactDigham\"", "\"ExactHerath\"", "\"ExactSun\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateIntegerMN\"", "\"ApproximateAsymptotic\"", "{\"ApproximateSwitchedMN\", mn} (where mn is the switching point)", "\"ApproximateSmallPf\"", "\"ApproximateLargeSNR\""}]<>"\n\n"<>DiversityTypeHelp[NakagamiProbabilityOfDetection]<>"\n\n"<>TimingHelp[NakagamiProbabilityOfDetection];
NakagamiProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,m_,OptionsPattern[]]:=Module[{n = 1, RelevantOptions},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[NakagamiProbabilityOfDetection]]]
]
NakagamiProbabilityOfDetection[M_,\[Gamma]_,\[Lambda]_,m_,n_,OptionsPattern[]]:=Module[{RelevantOptions, method, mn, limit, f, totaltime = 0, iterations = 0, time, result},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	limit = NakagamiLimit[M, \[Gamma], \[Lambda], m, n, Method->method, RelevantOptions[NakagamiLimit]];

	f := Which[
		(* Catch extreme values - they can cause errors *)
		\[Lambda] == -\[Infinity] || M == \[Infinity] || \[Gamma] == \[Infinity] || n == \[Infinity],
			1,
		method == "ExactAnnamalai",
			AnnamalaiNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,Limit->limit,RelevantOptions[AnnamalaiNakagamiProbabilityOfDetection]],
		method == "ExactDigham",
			DighamNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[DighamNakagamiProbabilityOfDetection]],
		method == "ExactHerath",
			HerathNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,Limit->limit,RelevantOptions[HerathNakagamiProbabilityOfDetection]],
		method == "ExactSun",
			SunNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,Limit->limit,RelevantOptions[SunNakagamiProbabilityOfDetection]],
		method == "ExactNumerical",
			NumericalNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NumericalNakagamiProbabilityOfDetection]],
		method == "ApproximateIntegerMN",
			IntegerMNNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[IntegerMNNakagamiProbabilityOfDetection]],
		(* Keep two names for legacy reasons *)
		method == "ApproximateLargeSNR" || method == "ApproximateSmallPf",
			LargeSNRNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[LargeSNRNakagamiProbabilityOfDetection]],
		method == "ApproximateAsymptotic",
			AsymptoticNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[AsymptoticNakagamiProbabilityOfDetection]],
		(* Keep two names for legacy reasons *)
		method == "Lopez-Benitez" || method == "ApproximateLargeMN",
			LopezBenitezNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[LopezBenitezNakagamiProbabilityOfDetection]],
		method == "ApproximateSwitchedMN",
			If[m n >= mn,
				AsymptoticNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[AsymptoticNakagamiProbabilityOfDetection]],
				IntegerMNNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[IntegerMNNakagamiProbabilityOfDetection]]
			],
		(* Keep two names for legacy reasons *)
		method == "ApproximateNumerical" || method == "NGaussian",
			NGaussianNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NGaussianNakagamiProbabilityOfDetection]],
		method == "ApproximateNumericalLowSNR",
			NGaussianLowSNRNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NGaussianLowSNRNakagamiProbabilityOfDetection]],
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


Options[NakagamiLimit] = {Method->"ExactAnnamalai", DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType], Tolerance->10^-6};
NakagamiLimit::usage = "NakagamiLimit[M, \[Gamma], \[Lambda], m] calculates the truncation point for use in the specified method for a single energy detector operating on a Nakagami-m channel.
NakagamiLimit[M, \[Gamma], \[Lambda], m, n] calculates the truncation point for use in the specified method for energy detection with diversity reception in a Nakagami-m channel.\n\n" <> MethodHelp[NakagamiLimit, {"\"ExactNumerical\"", "\"ExactAnnamalai\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateIntegerKN\"", "\"ApproximateAsymptotic\""}] <> "\n\n" <> DiversityTypeHelp[NakagamiLimit];
NakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1}, NakagamiLimit[M, \[Gamma], \[Lambda], m, n, DiversityType->"None", Tolerance->OptionValue[Tolerance]]]
NakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, \[Gamma]t, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType], method, mn},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NakagamiLimit][[All,1]]],Options[target][[All,1]]];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	Which[
		(* Catch extreme values - they can cause errors *)
		\[Lambda] == -\[Infinity] || M == \[Infinity] || \[Gamma] == \[Infinity] || n == \[Infinity],
			Undefined,
		method == "ExactAnnamalai",
			AnnamalaiNakagamiLimit[M, \[Gamma], \[Lambda], m, n, RelevantOptions[AnnamalaiNakagamiLimit]],
		method == "ExactHerath",
			HerathNakagamiLimit[M, \[Gamma], \[Lambda], m, n, RelevantOptions[HerathNakagamiLimit]],
		method == "ExactSun",
			SunNakagamiLimit[M, \[Gamma], \[Lambda], m, n, RelevantOptions[SunNakagamiLimit]],
		method == "ApproximateIntegerMN",
			IntegerMNNakagamiLimit[M, \[Gamma], \[Lambda], m, n, RelevantOptions[IntegerMNNakagamiLimit]],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Annamalai' s method*)


Options[AnnamalaiNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType], Limit->Null};
AnnamalaiNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[AnnamalaiNakagamiProbabilityOfDetection, "Annamalai"];
AnnamalaiNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AnnamalaiNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	AnnamalaiNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[AnnamalaiNakagamiProbabilityOfDetection]]]
]
AnnamalaiNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = NakagamiLimit[M, \[Gamma], \[Lambda], m, n, Method->"ExactAnnamalai", DiversityType->diversityType]];

	Which[
		diversityType == "None",
			1 - ((m / (m + (M / 2) \[Gamma]0))^m) (1 - GammaRegularized[M / 2, \[Lambda] / 2]) - Total[Table[(Gamma[m + k] / (Gamma[m] Gamma[k+1])) ((m / (m + (M / 2) \[Gamma]0))^m) ((((M / 2) \[Gamma]0) / (m + (M / 2) \[Gamma]0))^k) (1 - GammaRegularized[M / 2 + k, \[Lambda] / 2]),{k, 1, limit}]],
		diversityType == "MRC",
			1 - (2m / (2m + M \[Gamma]0))^(m n) (1 - GammaRegularized[M / 2, \[Lambda] / 2]) - Total[Table[Gamma[m n+k]/(Gamma[m n]Gamma[k+1]) (m/(m+M/2 \[Gamma]0))^(m n) (((M/2 \[Gamma]0)/(m+M/2 \[Gamma]0))^k) (1-GammaRegularized[M / 2 + k, \[Lambda] / 2]),{k,1,limit}]],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			1 - (m / (m + (M / 2) \[Gamma]0))^(m n) (1 - GammaRegularized[M n/2, \[Lambda] / 2]) - Total[Table[Gamma[m n+k]/(Gamma[m n]Gamma[k+1]) (m/(m+M/2 \[Gamma]0))^(m n) (((M/2 \[Gamma]0)/(m+M/2 \[Gamma]0))^k) (1-GammaRegularized[M n/2+k, \[Lambda] / 2]),{k,1,limit}]],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - AnnamalaiNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,Limit->limit[[i]],DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - AnnamalaiNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,Limit->limit,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


Options[AnnamalaiNakagamiLimit] = {DiversityType->OptionValue[AnnamalaiNakagamiProbabilityOfDetection,DiversityType], Tolerance->OptionValue[NakagamiLimit,Tolerance]};
AnnamalaiNakagamiLimit::usage = GenerateTruncationHelp[AnnamalaiNakagamiLimit, "Annamalai"];
AnnamalaiNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1},AnnamalaiNakagamiLimit[M,\[Gamma],\[Lambda],m,n,DiversityType->"None",Tolerance->OptionValue[Tolerance]]]
AnnamalaiNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType], f},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	Quiet[
		Which[
			diversityType == "None",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol] // Round;
				f[j_] := 1 - GammaRegularized[(M / 2) + j + 1, \[Lambda] / 2];
				Which[
					f[j0] > tol,
						While[f[j0] > tol, j0++];
						j0,
					f[j0] < tol,
						While[f[j0] < tol, j0--];
						j0,
					f[j0] == tol,
						j0,
					True,
						Undefined
				],
			diversityType == "MRC",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol] // Round;
				f[j_] := 1 - GammaRegularized[(M / 2) + j + 1, \[Lambda] / 2];
				Which[
					f[j0] > tol,
						While[f[j0] > tol, j0++];
						j0,
					f[j0] < tol,
						While[f[j0] < tol, j0--];
						j0,
					f[j0] == tol,
						j0,
					True,
						Undefined
				],
			diversityType == "EGC",
				Undefined,
			diversityType == "SC",
				Undefined,
			diversityType == "SEC",
				Undefined,
			diversityType == "SLC",
				j0 = (\[Lambda] / 2) - (M n / 2) - Sqrt[M n / 2] InverseQ[1 - tol] // Round;
				f[j_] := 1 - GammaRegularized[(M / 2) n + j + 1, \[Lambda] / 2];
				Which[
					f[j0] > tol,
						While[f[j0] > tol, j0++];
						j0,
					f[j0] < tol,
						While[f[j0] < tol, j0--];
						j0,
					f[j0] == tol,
						j0,
					True,
						Undefined
				],
			diversityType == "SLS",
				Table[AnnamalaiNakagamiLimit[M,\[Gamma][[i]],\[Lambda],m,DiversityType->"None",Tolerance->tol], {i, Length[\[Gamma]]}],
			True,
				Undefined
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*Digham' s method*)


Options[DighamNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
DighamNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[DighamNakagamiProbabilityOfDetection, "Digham"];
DighamNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[DighamNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	DighamNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[DighamNakagamiProbabilityOfDetection]]
]
DighamNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{limit, \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			Module[{A1, \[Beta]},
				A1 = Exp[-\[Lambda] \[Beta] / (2 m)] (\[Beta]^(m - 1) LaguerreL[m - 1, -\[Lambda] (1 - \[Beta]) / 2] + (1 - \[Beta]) Total[Table[\[Beta]^i LaguerreL[i,-\[Lambda] (1 - \[Beta]) / 2], {i, 0, m - 2}]]);
				\[Beta] = (2m) / (2m + M \[Gamma]0);
				A1 + \[Beta]^(m) Exp[-\[Lambda] / 2] Total[Table[((\[Lambda] / 2)^i / i!) Hypergeometric1F1[m, i + 1, \[Lambda] (1 - \[Beta]) / 2], {i, 1, (M / 2) - 1}]]
			],
		diversityType == "MRC",
			Undefined,
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			Module[{A1, \[Beta]},
				A1 = Exp[-\[Lambda] \[Beta] / (2 m n)] (\[Beta]^(m n - 1) LaguerreL[m n - 1, -\[Lambda] (1 - \[Beta]) / 2] + (1 - \[Beta]) Total[Table[\[Beta]^i LaguerreL[i,-\[Lambda] (1 - \[Beta]) / 2], {i, 0, m n - 2}]]);
				\[Beta] = (2m) / (2m + M \[Gamma]0);
				A1 + \[Beta]^(m n) Exp[-\[Lambda] / 2] Total[Table[((\[Lambda] / 2)^i / i!) Hypergeometric1F1[m n, i + 1, \[Lambda] (1 - \[Beta]) / 2], {i, 1, (M n / 2) - 1}]]
			],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - DighamNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - DighamNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Herath' s method*)


Options[HerathNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType], Limit->Null};
HerathNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[HerathNakagamiProbabilityOfDetection, "Herath"];
HerathNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[HerathNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	HerathNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[HerathNakagamiProbabilityOfDetection]]
]
HerathNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Gamma]0, \[Gamma]t, \[Psi], diversityType = OptionValue[DiversityType]},
	limit = HerathNakagamiLimit[M, \[Gamma], \[Lambda], m, n, DiversityType->diversityType];

	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = NakagamiLimit[M, \[Gamma], \[Lambda], m, n, Method->"ExactHerath", DiversityType->diversityType]];

	Which[
		diversityType == "None",
			1 - Exp[-\[Lambda] / 2] (m / ((M / 2) \[Gamma]0 + m))^(m) Total[Table[((\[Lambda] / 2)^j / j!) Hypergeometric1F1[m, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))],{j, (M / 2), limit}]],
		diversityType == "MRC",
			1 - Exp[-\[Lambda] / 2] (m / ((M / 2) \[Gamma]0 + m))^(n m) Total[Table[((\[Lambda] / 2)^j / j!) Hypergeometric1F1[m n, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))],{j, (M / 2), limit}]],
		!ListQ[diversityType] && diversityType == "EGC",
			Which[
				n == 2,
					1 - (2 Sqrt[\[Pi]] Gamma[2m] (m)^(2m) / (Gamma[m]^(2) Gamma[2m + (1 / 2)] 2^(4m - 1))) Exp[-\[Lambda] / 2] (2 / ((M / 2) \[Gamma]0))^(2m) Total[Table[(\[Lambda] / 2)^(j / 2) NIntegrate[x^(2m - (j / 2) - 1) Exp[-(2m / ((M / 2) \[Gamma]0) + 1) x] BesselI[j, Sqrt[2 \[Lambda] x]] Hypergeometric1F1[2m, 2m + (1 / 2), m x / ((M / 2) \[Gamma]0)],{x, 0, \[Infinity]}],{j, M / 2, limit}]],
				n == 3,
					(* Untested *)
					(*1 - Sqrt[\[Pi]] Exp[-\[Lambda] / 2] Total[Table[Total[Table[Total[Table[(\[Lambda] / 2)^j (3m / ((M / 2) \[Gamma]0 + 3m))^(3m + p + k) ((Gamma[2m + p] Gamma[2m + k] Gamma[3m + p] Gamma[3m + p + 1 / 2] Gamma[4m + 2p + k]) / (Gamma[m]^3 Gamma[n + 1] Gamma[p + 1] Gamma[2m + p + 1 / 2] Gamma[3m + p + k + 1 / 2] Gamma[6m + 2p])) (1 / (2^(4m + p + k - 3) k!)) Hypergeometric1F1[3m + p + k, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + 3m))], {k, 0, limit[[1]]}]], {p, 0, limit[[2]]}]], {j, M / 2, limit[[3]]}]]*)
					Undefined,
				True,
					Undefined
			],
		diversityType == "SC",
			(* Untested *)
			(*1 - n Exp[-\[Lambda] / 2] (m / \[Gamma]0)^(m) Total[Table[Total[Table[Binomial[n - 1, k] (-1)^(k) ((\[Lambda] / 2)^(j) / (j!)) Total[Table[MultinomialCoefficient[k, i, m] Pochhammer[m, i] (\[Gamma]0 / (\[Gamma]0 + m (k + 1)))^(i + m) Hypergeometric1F1[i + m, j + 1, \[Lambda] \[Gamma]0 / (2 (\[Gamma]0 + m (k + 1)))],{i, 0, k*(m - 1)}]],{k, 0, n - 1}]],{j, M / 2, limit}]],*)
			Undefined,
		diversityType == "SEC" && n == 2,
			(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0]) (1 - Exp[-\[Lambda] / 2] (m / (m + (M / 2) \[Gamma]0))^(m) Total[Table[((\[Lambda] / 2)^(j) / j!) Hypergeometric1F1[m, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))],{j, M / 2, limit[[1]]}]]) + (m / ((M / 2) \[Gamma]0 + m))^(m) (Exp[-\[Lambda] / 2] / Gamma[m]) Total[Table[Total[Table[((M / 2) \[Gamma]0 / ((M / 2) \[Gamma]0 + m))^(j) Gamma[j + m, (1 + m / ((M / 2) \[Gamma]0)) (M / 2) \[Gamma]t] (\[Lambda] / 2)^(k) / (j! k!),{k, 0, j + (M / 2) - 1}]],{j, 0, limit[[2]]}]],
		diversityType == "SLC",
			1 - Exp[-\[Lambda] / 2] (m / ((M / 2) \[Gamma]0 + m))^(n m) Total[Table[((\[Lambda] / 2)^j / j!) Hypergeometric1F1[m n, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))],{j, (M n / 2), limit}]],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - HerathNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,Limit->limit[[i]],DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - HerathNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,Limit->limit,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


Options[HerathNakagamiLimit] = {DiversityType->OptionValue[HerathNakagamiProbabilityOfDetection,DiversityType], Tolerance->OptionValue[NakagamiLimit,Tolerance]};
HerathNakagamiLimit::usage = GenerateTruncationHelp[HerathNakagamiLimit, "Herath"];
HerathNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1},HerathNakagamiLimit[M,\[Gamma],\[Lambda],m,n,Tolerance->OptionValue[Tolerance]]]
HerathNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType], \[Psi]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType]//N;

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			j0 = (\[Lambda] / 2) - 1 - InverseQ[1 - tol];
			j/.FindRoot[(m / ((M / 2) \[Gamma]0 + m))^(m) Hypergeometric1F1[m, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]) == tol,{j, j0, 1, \[Infinity]}],
		diversityType == "MRC",
			j0 = (\[Lambda] / 2) - 1 - InverseQ[1 - tol];
			j/.FindRoot[(m / ((M / 2) \[Gamma]0 + m))^(m) Hypergeometric1F1[m, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]) == tol,{j, j0, 1, \[Infinity]}],
		diversityType == "EGC",
			j0 = Max[M / 2, (\[Lambda] / 2) - 1 - InverseQ[1 - tol]];
			\[Psi][\[Alpha]_,\[Beta]_,\[Gamma]0_,\[Gamma]1_,x_,y_]:=NSum[Pochhammer[\[Alpha],m0+n0] Pochhammer[\[Beta],m0] / (Pochhammer[\[Gamma]0, m0] Pochhammer[\[Gamma]1, n0]) (x^m0 / m0!) (y^n0 / n0!),{m0,0,\[Infinity]},{n0,0,\[Infinity]}];
			Which[
				n == 2,
					j/.FindRoot[(Gamma[2m]^2 Gamma[1 / 2] / (Gamma[m]^(2) Gamma[2m + (1 / 2)] 2^(2m - 2))) (m / ((M / 2) \[Gamma]0 + 2m))^(2m) \[Psi][2m, 2m, 2m + 1 / 2, j + 1, m / ((M / 2) \[Gamma]0 + 2m), \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + 2m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]) == tol, {j, j0, M / 2, \[Infinity]}],
				n == 3,
					(* Untested *)
					(*{100, 100, M / 2 + 100}*)
					Undefined,
				True,
					Undefined
			],
		diversityType == "SC",
			(* Untested *)
			(*j0 = (\[Lambda] / 2) - 1 - InverseQ[1 - tol];
			j/.FindRoot[n (m / \[Gamma]0)^(m) (1 - GammaRegularized[j + 1, \[Lambda] / 2]) Total[Table[Binomial[n - 1, k] Total[Table[MultinomialCoefficient[k, i, m] Pochhammer[m, i] (\[Gamma]0 / (\[Gamma]0 + m (k + 1)))^(i + m) Hypergeometric1F1[i + m, j + 1, \[Lambda] \[Gamma]0 / (2 (\[Gamma]0 + m (k + 1)))],{i, 0, k*(m - 1)}]],{k, 0, n - 1}]] == tol,{j, j0, 1, \[Infinity]}]*)
			Undefined,
		diversityType == "SEC",
			j0 = {(\[Lambda] / 2) - 1 - InverseQ[1 - tol], Null};
			{j/.FindRoot[(1 - GammaRegularized[m, m (M / 2) \[Gamma]t / ((M / 2) \[Gamma]0)]) Hypergeometric1F1[m, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]) (m / (m + (M / 2) \[Gamma]0))^(m) == tol,{j, j0[[1]], 1, \[Infinity]}], InverseCDF[NegativeBinomialDistribution[m, m / ((M / 2) \[Gamma]0 + m)], 1 - tol]},
		diversityType == "SLC",
			j0 = (\[Lambda] / 2) - 1 - InverseQ[1 - tol];
			j/.FindRoot[(m / ((M / 2) \[Gamma]0 + m))^(m n) Hypergeometric1F1[m n, j + 1, \[Lambda] (M / 2) \[Gamma]0 / (2 ((M / 2) \[Gamma]0 + m))] (1 - GammaRegularized[j + 1, \[Lambda] / 2]) == tol,{j, j0, 1, \[Infinity]}],
		diversityType == "SLS",
			Table[HerathNakagamiLimit[M,\[Gamma][[i]],\[Lambda],m,DiversityType->"None", Tolerance->tol], {i, Length[\[Gamma]]}],
		True,
			Undefined
	]//Ceiling
]


(* ::Subsubsection::Closed:: *)
(*Sun' s method*)


Options[SunNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType], Limit->Null};
SunNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[SunNakagamiProbabilityOfDetection, "Sun"];
SunNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[SunNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	SunNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[SunNakagamiProbabilityOfDetection]]
]
SunNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	If[limit==Null, limit = NakagamiLimit[M, \[Gamma], \[Lambda], m, n, Method->"ExactSun", DiversityType->diversityType]];

	Which[
		diversityType == "None",
			GammaRegularized[M / 2, \[Lambda] / 2] + Exp[-\[Lambda] / 2] Total[Table[((\[Lambda] / 2)^(j) / j!) (1 - (m / (m + (M / 2) \[Gamma]0))^(m) Total[Table[(m + k - 1)! / (Gamma[m] k!) ((M / 2) \[Gamma]0 / (m + (M / 2) \[Gamma]0))^(k),{k, 0, j - M / 2}]]),{j, M / 2, M / 2 + limit}]],
		diversityType == "MRC",
			1 - (m / (m + (M / 2) \[Gamma]0))^(m n) Total[Table[(1 - GammaRegularized[(M / 2) + j, \[Lambda] / 2]) Pochhammer[m n, j] / j! ((M / 2) \[Gamma]0 / (m + (M / 2) \[Gamma]0))^(j),{j, 0, limit}]],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			1 - n Total[Table[(-1)^(l) Binomial[n - 1, l] Total[Table[MultinomialCoefficient[l, k, m] (m / ((M / 2) \[Gamma]0))^(m + k) Total[Table[(1 - GammaRegularized[i + M / 2, \[Lambda] / 2]) Pochhammer[m, i + k] / i! ((M / 2) \[Gamma]0 / ((M / 2) \[Gamma]0 + m (l + 1)))^(i + m + k),{i, 0, limit}]],{k, 0, l (m - 1)}]],{l, 0, n - 1}]],
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			1 - (m / (m + (M / 2) \[Gamma]0))^(m n) Total[Table[(1 - GammaRegularized[(M / 2) n + j, \[Lambda] / 2]) Pochhammer[m n, j] / j! ((M / 2) \[Gamma]0 / (m + (M / 2) \[Gamma]0))^(j),{j, 0, limit}]],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - SunNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,Limit->limit[[i]],DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - SunNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,Limit->limit,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


Options[SunNakagamiLimit] = {DiversityType->OptionValue[SunNakagamiProbabilityOfDetection,DiversityType], Tolerance->OptionValue[NakagamiLimit,Tolerance]};
SunNakagamiLimit::usage = GenerateTruncationHelp[SunNakagamiLimit, "Sun"];
SunNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{n = 1},SunNakagamiLimit[M,\[Gamma],\[Lambda],m,n,Tolerance->OptionValue[Tolerance]]]
SunNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, j, j0, tol = OptionValue[Tolerance], diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	Quiet[
		Which[
			diversityType == "None",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol];
				j/.FindRoot[(1 - GammaRegularized[M / 2 + j - 1, \[Lambda] / 2]) (1 - CDF[NegativeBinomialDistribution[m, (m / (m + (M / 2) \[Gamma]))^(m)], j + 1]) == tol,{j, j0, 1, \[Infinity]}],
			diversityType == "MRC",
				AnnamalaiNakagamiLimit[M, \[Gamma], \[Lambda], m, n, DiversityType->diversityType, Tolerance->tol],
			diversityType == "EGC",
				Undefined,
			diversityType == "SC",
				j0 = (\[Lambda] / 2) - (M / 2) - Sqrt[M / 2] InverseQ[1 - tol];
				j/.FindRoot[n Total[Table[(-1)^(l) Binomial[n - 1, l] Total[Table[MultinomialCoefficient[l, k, m] (Gamma[m + k] / Gamma[m]) (1 - GammaRegularized[M / 2 + j, \[Lambda] / 2]) (1 / (l + 1))^(m + k), {k, 0, l (m - 1)}]], {l, 0, n - 1}]] == tol, {j, j0, 1, \[Infinity]}],
			diversityType == "SEC",
				Undefined,
			diversityType == "SLC",
				AnnamalaiNakagamiLimit[M, \[Gamma], \[Lambda], m, n, DiversityType->diversityType, Tolerance->tol],
			diversityType == "SLS",
				Table[SunNakagamiLimit[M,\[Gamma][[i]],\[Lambda],m,DiversityType->"None",Tolerance->tol], {i, Length[\[Gamma]]}],
			True,
				Undefined
		]//Ceiling
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method*)


Options[NumericalNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
NumericalNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[NumericalNakagamiProbabilityOfDetection, "Numerical"];
NumericalNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NumericalNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NumericalNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,RelevantOptions[NumericalNakagamiProbabilityOfDetection]]
]
NumericalNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType], \[ScriptCapitalD], x},
	\[ScriptCapitalD] = NumericalNakagamiDistribution[M, \[Gamma], m, n, diversityType];

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
					1 - Product[1 - NumericalNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NumericalNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


NumericalNakagamiDistribution[M_?NumericQ, \[Gamma]_, m_?NumericQ, n_?IntegerQ, diversityType0_] := NumericalNakagamiDistribution[M, \[Gamma], m, n, diversityType0] = Module[{\[Gamma]t, \[Gamma]0, g, numberOfPoints = 100000, x, diversityType},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType0];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[a_,b_] := SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M a, M x], x \[Distributed] GammaDistribution[m b, \[Gamma]0 / m]], numberOfPoints]];

	Which[
		diversityType == "None",
			g[1, 1],
		diversityType == "MRC",
			g[1, n],
		diversityType == "EGC",
			x = Table[Unique[], {n}];
			SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M, M Sum[Sqrt[x[[i]]], {i, n}]^2/n], Table[x[[i]] \[Distributed] GammaDistribution[m, \[Gamma]0 / m], {i, n}]], numberOfPoints]],
		diversityType == "SC",
			SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M, M x], x \[Distributed] ProbabilityDistribution[Piecewise[{{(n / Gamma[m]) (m / \[Gamma]0)^m x^(m - 1) Exp[-m x / \[Gamma]0] (1 - GammaRegularized[m, m x / \[Gamma]0])^(n - 1), x > 0}}], {x, 0, \[Infinity]}]], numberOfPoints]],
		diversityType == "SEC",
			SmoothKernelDistribution[RandomVariate[ParameterMixtureDistribution[NoncentralChiSquareDistribution[M, M x], x \[Distributed] ProbabilityDistribution[Piecewise[{{(1 - GammaRegularized[m, m \[Gamma]t/\[Gamma]0])^(n - 1) PDF[GammaDistribution[m, \[Gamma]0/m], x], 0 <= x < \[Gamma]t}, {Sum[(1 - GammaRegularized[m, m \[Gamma]t/\[Gamma]0])^j, {j, 0, n - 1}] PDF[GammaDistribution[m, \[Gamma]0/m], x], x >= \[Gamma]t}}], {x, 0, \[Infinity]}]], numberOfPoints]],
		diversityType == "SLC",
			g[n, n],
		diversityType == "SLS",
			(* We don't need the SLS PDF here, the master function will use the no diversity PDF instead *)
			Null,
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Integer mn method*)


Options[IntegerMNNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType], Limit->Null};
IntegerMNNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[IntegerMNNakagamiProbabilityOfDetection, "IntegerMN"];
IntegerMNNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[IntegerMNNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	IntegerMNNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[IntegerMNNakagamiProbabilityOfDetection]]]
]
IntegerMNNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{limit = OptionValue[Limit], \[Gamma]0, \[Gamma]t, g, diversityType = OptionValue[DiversityType], x = Round[m n], tol = 10^-6},
	(* This method can only be used when m * n is an integer *)
	If[Abs[m n - x] > tol, Return[Undefined]];
	
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	(* Precision is set to 100 here - required for FaddeevaDerivative to be stable for large orders *)
	g[k_, a_, b_, c_] := (1 / 2) Exp[-a^2] Total[Table[((I c / (2 b))^l / l!) FaddeevaDerivative[l, N[-I (a + (c / (2 b))), 100]], {l, 0, k - 1}]] // Re;
	
	Which[
		diversityType == "None",
			AWGNProbabilityOfFalseAlarm[M, \[Lambda], Method->"Approximate", DiversityType->diversityType] + g[m, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / \[Gamma]0],
		diversityType == "MRC",
			AWGNProbabilityOfFalseAlarm[M, \[Lambda], Method->"Approximate", DiversityType->diversityType] + g[m n, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / \[Gamma]0],
		diversityType == "EGC",
			(* Dharmawansa's exact method *)
			(*If[limit==Null, limit = NakagamiLimit[M, \[Gamma], \[Lambda], m, n, Method->"ApproximateIntegerMN", DiversityType->diversityType]];
			Which[
				n == 2,
					With[{A = (\[Lambda] - M) / (2 Sqrt[M]), B = - Sqrt[M] \[Gamma]0 / (2 m n)},
						AWGNProbabilityOfFalseAlarm[M, \[Lambda]] + (2 Sqrt[\[Pi]] / 2^(4m - 1)) (Gamma[2m] / (Gamma[m]^2 Gamma[2m + 1 / 2])) Total[Table[(Pochhammer[2m, s] / Pochhammer[2m + 1 / 2, s]) (Gamma[2m + s] / s!) (1 / 2)^s g[A, B, 2m + s], {s, 0, limit}]]
					],
				n == 3,
					(* Untested *)
					(4 Sqrt[\[Pi]] / (2^(4m - 1))) (Gamma[2m] / (Gamma[m]^3)) Total[Table[(Gamma[2m+l] Gamma[4m+2l] / (Gamma[2m+l+1/2] Gamma[6m+2l] Gamma[l+1] (2^l))) Total[Table[(Pochhammer[2m,k] Pochhammer[4m+2l,k] / (Pochhammer[3m+l+1/2,k] Pochhammer[3m+l,k]))((1/2)^k/k!)(-1)^(k+2+3 m) (*D[(1+Erf[(M-\[Lambda])/(2 Sqrt[M])]+E^((3 m t (3 m t+M \[Gamma]0-\[Gamma]0 \[Lambda]))/(M \[Gamma]0^2)) Erfc[(6 m t+M \[Gamma]0-\[Gamma]0 \[Lambda])/(2 Sqrt[M] \[Gamma]0)])/(2 t),{t,3m+l+k-1}]*),{k,0,limit[[2]]}]],{l,0,limit[[1]]}]],
				True,
					Undefined
			],*)
			(* Nakagami's approximate method *)
			With[{\[Beta] = (Gamma[m + 1]^2 + m (n - 1) Gamma[m + 1 / 2]^2) / (n Gamma[m + 1]^2)},
				AWGNProbabilityOfFalseAlarm[M, \[Lambda], Method->"Approximate", DiversityType->diversityType] + g[m n, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / (\[Beta] \[Gamma]0)]
			],
		diversityType == "SC",
			With[{A = (\[Lambda] - M) / (2 Sqrt[M]), B = - Sqrt[M] \[Gamma]0 / (2 m)},
				AWGNProbabilityOfFalseAlarm[M, \[Lambda], Method->"Approximate", DiversityType->diversityType] + (n / Gamma[m]) Total[Table[(-1)^l Binomial[n - 1, l] Total[Table[MultinomialCoefficient[l, k, m] (1 / (l + 1)^(m + k)) g[A, B / (l + 1), m + k], {k, 0, l (m - 1)}]], {l, 0, n - 1}]]
			],
		diversityType == "SEC",
			With[{A = (\[Lambda] - M) / (2 Sqrt[M]), B = - Sqrt[M] \[Gamma]0 / (2 m), C = (\[Lambda] - M (1 + \[Gamma]t)) / (2 Sqrt[M])},
				(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0])^(n - 1) IntegerMNNakagamiProbabilityOfDetection[M, \[Gamma]0, \[Lambda], m, DiversityType->"None"] + Sum[(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0])^j, {j, 0, n - 2}] (GammaRegularized[m, m \[Gamma]t / \[Gamma]0] ((1 / 2) (1 - Erf[C])) + ((1 / 2) Exp[-C^2] Total[Table[((I/(2 B))^l / l!) GammaRegularized[m - l, m \[Gamma]t / \[Gamma]0] FaddeevaDerivative[l, N[-I (C + (1 / (2 B))),30]], {l, 0, m - 1}]] // Re))
			],
		diversityType == "SLC",
			AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType] + g[m n, (\[Lambda] - M n) / (2 Sqrt[M n]), - Sqrt[M] / (2 Sqrt[n]), m / \[Gamma]0],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - IntegerMNNakagamiProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], m, DiversityType->"None"], {i, n}],
				!ListQ[\[Gamma]0],
					1 - (1 - IntegerMNNakagamiProbabilityOfDetection[M, \[Gamma]0, \[Lambda], m, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]




(* This function applies to Dharmawansa's method only *)
Options[IntegerMNNakagamiLimit]={DiversityType->OptionValue[IntegerMNNakagamiProbabilityOfDetection,DiversityType], Tolerance->OptionValue[NakagamiLimit,DiversityType]};
IntegerMNNakagamiLimit::usage = GenerateTruncationHelp[IntegerMNNakagamiLimit, "IntegerMN"];
IntegerMNNakagamiLimit[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{j, j0, tol = OptionValue[Tolerance]},
	Which[
		n == 2,
			j0 = 0;
			Quiet[j/.FindRoot[(1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda]]) (1 / Gamma[m]^2) 2^(1 - 4m - j) Sqrt[\[Pi]] Gamma[2m + j + 1]^2 HypergeometricPFQRegularized[{1, 2m + j + 1, 2m + j + 1},{2 + j, 2m + 3 / 2 + j}, 1 / 2] == tol, {j, j0, 0, \[Infinity]}]],
		n == 3,
			Module[{l0, k0, f},
				k0 = 0;
				l0 = 0;
				While[
					(4 Sqrt[\[Pi]] Gamma[2m] / (Gamma[m]^3 2^(4m))) Sum[((Gamma[2m + l] Gamma[4m + 2l]) / (Gamma[2m + l + 1 / 2] Gamma[6m + 2l])) ((1 / 2)^l / l!) Gamma[3m + l] Sum[(Pochhammer[2m, k] Pochhammer[4m + 2l, k] / Pochhammer[3m + l + 1 / 2, k]) ((1 / 2)^k / k!), {k, k0 + 1, \[Infinity]}], {l, 0, \[Infinity]}] > tol,
					k0++
				];
				f[l_?NumericQ]:=Total[Table[(Pochhammer[2m, k] Pochhammer[4m + 2l, k] / Pochhammer[3m + l + 1 / 2, k]) ((1 / 2)^k / k!), {k, 0, k0}]];
				While[
					(4 Sqrt[\[Pi]] Gamma[2m] / (Gamma[m]^3 2^(4m))) NSum[((Gamma[2m + l] Gamma[4m + 2l]) / (Gamma[2m + l + 1 / 2] Gamma[6m + 2l])) ((1 / 2)^l / l!) Gamma[3m + l] f[l], {l, l0 + 1, \[Infinity]}] > tol,
					l0++
				];
				{l0, k0}
			],
		True,
			Undefined
	]//Ceiling
]


(* ::Subsubsection::Closed:: *)
(*Large SNR method*)


Options[LargeSNRNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
LargeSNRNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[LargeSNRNakagamiProbabilityOfDetection, "LargeSNR"];
LargeSNRNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[LargeSNRNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	LargeSNRNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[LargeSNRNakagamiProbabilityOfDetection]]]
]
LargeSNRNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{limit, g, \[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType], x = Round[m n], tol = 10^-6},
	(* This method can only be used when m * n is an integer *)
	If[Abs[m n - x] > tol, Return[Undefined]];
	
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[k_, a_, b_, c_] := AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType] + (1 - AWGNProbabilityOfFalseAlarm[M, \[Lambda], n, Method->"Approximate", DiversityType->diversityType]) GammaRegularized[k, - a c / b];

	Which[
		diversityType == "None",
			g[m, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / \[Gamma]0],
		diversityType == "MRC",
			g[m n, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / \[Gamma]0],
		diversityType == "EGC",
			With[{\[Beta] = (Gamma[m + 1]^2 + m (n - 1) Gamma[m + 1 / 2]^2) / (n Gamma[m + 1]^2)},
				g[m n, (\[Lambda] - M) / (2 Sqrt[M]), - Sqrt[M] / 2, m / (\[Beta] \[Gamma]0)]
			],
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			g[m n, (\[Lambda] - M n) / (2 Sqrt[M n]), - Sqrt[M] / (2 Sqrt[n]), m / \[Gamma]0],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - LargeSNRNakagamiProbabilityOfDetection[M, \[Gamma]0[[i]], \[Lambda], m, DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - LargeSNRNakagamiProbabilityOfDetection[M, \[Gamma]0, \[Lambda], m, DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]




(* ::Subsubsection::Closed:: *)
(*Asymptotic method*)


Options[AsymptoticNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
AsymptoticNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[AsymptoticNakagamiProbabilityOfDetection, "Asymptotic"];
AsymptoticNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[AsymptoticNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	AsymptoticNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[AsymptoticNakagamiProbabilityOfDetection]]]
]
AsymptoticNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{g, \[Gamma]0, \[Gamma]t, time, result, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	g[A_,B_] := Q[A Sqrt[2 / (B^2 + 1)]];

	Which[
		diversityType == "None",
			g[(\[Lambda] - M (1 + \[Gamma]0)) / (2 Sqrt[M]), - Sqrt[M / (2 m)] \[Gamma]0],
		diversityType == "MRC",
			g[(\[Lambda] - M (1 + n \[Gamma]0)) / (2 Sqrt[M]), - Sqrt[M n / (2 m)] \[Gamma]0],
		diversityType == "EGC",
			With[{c = Gamma[m + 1]^2 / (Gamma[m + 1]^2 + m (n - 1) Gamma[m + 1 / 2]^2)},
				g[(\[Lambda] - M (1 + \[Gamma]0 / c)) / (2 Sqrt[M]), - Sqrt[M / (2 m n)] \[Gamma]0 / c]
			],
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Undefined,
		diversityType == "SLC",
			g[(\[Lambda] - M n (1 + \[Gamma]0)) / (2 Sqrt[M n]), - Sqrt[M / (2 m)] \[Gamma]0],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - AsymptoticNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - AsymptoticNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Lopez-Benitez asymptotic method*)


Options[LopezBenitezNakagamiProbabilityOfDetection] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
LopezBenitezNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[LopezBenitezNakagamiProbabilityOfDetection, "LopezBenitez"];
LopezBenitezNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[LopezBenitezNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	LopezBenitezNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[LopezBenitezNakagamiProbabilityOfDetection]]]
]
LopezBenitezNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?NumericQ,OptionsPattern[]]:=Module[{\[Gamma]0, \[Gamma]t, diversityType = OptionValue[DiversityType]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			(1/2 (1+Erf[(m (M (1+\[Gamma]0)-\[Lambda]))/(Sqrt[2] M Sqrt[m] \[Gamma]0)])+1/2 E^((4 c M (-m+a M \[Gamma]0^2)+b (-b M^2 \[Gamma]0^2+2 Sqrt[2] m Sqrt[M] (M (1+\[Gamma]0)-\[Lambda]))-2 a m (-M (1+\[Gamma]0)+\[Lambda])^2)/(4 M (-m+a M \[Gamma]0^2))) Sqrt[m/(m-a M \[Gamma]0^2)] (Erf[(b Sqrt[M^3] \[Gamma]0^2+Sqrt[2] m (-M (1+\[Gamma]0)+\[Lambda]))/(2 M \[Gamma]0 Sqrt[(m-a M \[Gamma]0^2)])]-Erf[(-2 m+\[Gamma]0 (-2 a M+Sqrt[2] b Sqrt[M]+2 a \[Lambda]))/(2 Sqrt[2] Sqrt[(m-a M \[Gamma]0^2)])])+1/2 E^(-((4 c M (m-a M \[Gamma]0^2)+b (b M^2 \[Gamma]0^2+2 Sqrt[2] m Sqrt[M] (M (1+\[Gamma]0)-\[Lambda]))+2 a m (-M (1+\[Gamma]0)+\[Lambda])^2)/(4 M (-m+a M \[Gamma]0^2)))) Sqrt[m/(m-a M \[Gamma]0^2)] (-2+Erfc[(b Sqrt[M^3] \[Gamma]0^2+Sqrt[2] m (M (1+\[Gamma]0)-\[Lambda]))/(2 M \[Gamma]0 Sqrt[(m-a M \[Gamma]0^2)])]))/.LopezBenitezParameters[(-M (1+\[Gamma]0)+\[Lambda])/(2 Sqrt[M])],
		diversityType == "MRC",
			(-((E^(-((b^2 M n \[Gamma]0^2+4 c (m-a M n \[Gamma]0^2)+(2 Sqrt[2] b m (M+M n \[Gamma]0-\[Lambda]))/Sqrt[M]+(2 a m (M+M n \[Gamma]0-\[Lambda])^2)/M)/(4 (-m+a M n \[Gamma]0^2)))) m (1+Erf[(b M^(3/2) n \[Gamma]0^2+Sqrt[2] m (M+M n \[Gamma]0-\[Lambda]))/(2 M \[Gamma]0 Sqrt[n (m-a M n \[Gamma]0^2)])]))/(2 Sqrt[m (m-a M n \[Gamma]0^2)]))+(E^(-((b^2 M n \[Gamma]0^2+4 c (m-a M n \[Gamma]0^2)-(2 Sqrt[2] b m (M+M n \[Gamma]0-\[Lambda]))/Sqrt[M]+(2 a m (M+M n \[Gamma]0-\[Lambda])^2)/M)/(4 (-m+a M n \[Gamma]0^2)))) m (Erf[(Sqrt[2] b M^(3/2) n \[Gamma]0^2-2 m (M+M n \[Gamma]0-\[Lambda]))/(2 M n Sqrt[-2 a M+(2 m)/(n \[Gamma]0^2)] \[Gamma]0^2)]-Erf[(-2 m+\[Gamma]0 (Sqrt[2] b Sqrt[M]-2 a M+2 a \[Lambda]))/(2 Sqrt[-2 a M+(2 m)/(n \[Gamma]0^2)] \[Gamma]0)]))/(2 Sqrt[m (m-a M n \[Gamma]0^2)])+1/2 Erfc[-((m (M+M n \[Gamma]0-\[Lambda]))/(Sqrt[2] M Sqrt[m n] \[Gamma]0))])/.LopezBenitezParameters[(-M (1+\[Gamma]0)+\[Lambda])/(2 Sqrt[M])],
		diversityType == "EGC",
			Undefined,
		diversityType == "SC",
			Undefined,
		diversityType == "SEC",
			Which[
				-((M-\[Lambda])/M) < \[Gamma]t,
					(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0]) LopezBenitezNakagamiProbabilityOfDetection[M, \[Gamma]0, \[Lambda], m] + (-(1/(2 Sqrt[2])) E^(-((b^2 M^2 \[Gamma]0^2+4 c M (m-a M \[Gamma]0^2)+2 Sqrt[2] b m Sqrt[M] (M+M \[Gamma]0-\[Lambda])+2 a m (M+M \[Gamma]0-\[Lambda])^2)/(4 M (-m+a M \[Gamma]0^2)))) Sqrt[m] ((Sqrt[2] (1+Erf[(Sqrt[2] m+\[Gamma]0 (b Sqrt[M]+Sqrt[2] a (M-\[Lambda])))/(2 Sqrt[m-a M \[Gamma]0^2])]))/Sqrt[m-a M \[Gamma]0^2]+((2 m (\[Gamma]0-\[Gamma]t)+\[Gamma]0^2 (Sqrt[2] b Sqrt[M]+2 a (M+M \[Gamma]t-\[Lambda]))) Erf[Abs[2 m (\[Gamma]0-\[Gamma]t)+\[Gamma]0^2 (Sqrt[2] b Sqrt[M]+2 a (M+M \[Gamma]t-\[Lambda]))]/(2 Sqrt[2] \[Gamma]0 Sqrt[m-a M \[Gamma]0^2])])/Sqrt[(m-a M \[Gamma]0^2) (2 m^2 (\[Gamma]0-\[Gamma]t)^2-2 m \[Gamma]0^2 (-\[Gamma]0+\[Gamma]t) (Sqrt[2] b Sqrt[M]+2 a (M+M \[Gamma]t-\[Lambda]))+\[Gamma]0^4 (b^2 M+2 Sqrt[2] a b Sqrt[M] (M+M \[Gamma]t-\[Lambda])+2 a^2 (M+M \[Gamma]t-\[Lambda])^2))]-(Sqrt[2] (2 m+\[Gamma]0 (Sqrt[2] b Sqrt[M]+2 a (M-\[Lambda]))) Erf[Abs[2 m+Sqrt[2] b Sqrt[M] \[Gamma]0+2 a M \[Gamma]0-2 a \[Gamma]0 \[Lambda]]/(2 Sqrt[2] Sqrt[m-a M \[Gamma]0^2])])/(Sqrt[m-a M \[Gamma]0^2] Abs[2 m+Sqrt[2] b Sqrt[M] \[Gamma]0+2 a M \[Gamma]0-2 a \[Gamma]0 \[Lambda]]))+1/2 Erfc[(Sqrt[m] (-\[Gamma]0+\[Gamma]t))/(Sqrt[2] \[Gamma]0)])/.LopezBenitezParameters[(-M (1+\[Gamma]0)+\[Lambda])/(2 Sqrt[M])],
				-(M-\[Lambda])/M >= \[Gamma]t,
					(1 - GammaRegularized[m, m \[Gamma]t / \[Gamma]0]) LopezBenitezNakagamiProbabilityOfDetection[M, \[Gamma]0, \[Lambda], m] + (1/2 (-Erf[(Sqrt[m] (-\[Gamma]0+\[Gamma]t))/(Sqrt[2] \[Gamma]0)]-Erf[(Sqrt[m] (M+M \[Gamma]0-\[Lambda]))/(Sqrt[2] M \[Gamma]0)])-1/(2 Sqrt[2-(2 a M \[Gamma]0^2)/m]) E^(-((b^2 M^2 \[Gamma]0^2+4 c M (m-a M \[Gamma]0^2)+2 Sqrt[2] b m Sqrt[M] (M+M \[Gamma]0-\[Lambda])+2 a m (M+M \[Gamma]0-\[Lambda])^2)/(4 M (-m+a M \[Gamma]0^2)))) (Sqrt[2] (1+Erf[(Sqrt[2] m+\[Gamma]0 (b Sqrt[M]+Sqrt[2] a (M-\[Lambda])))/(2 Sqrt[m-a M \[Gamma]0^2])])+((Sqrt[2] b M^(3/2) \[Gamma]0^2+2 m (M+M \[Gamma]0-\[Lambda])) Erf[Abs[Sqrt[2] b M^(3/2) \[Gamma]0^2+2 m (M+M \[Gamma]0-\[Lambda])]/(2 Sqrt[2] M \[Gamma]0 Sqrt[m-a M \[Gamma]0^2])])/Sqrt[b^2 M^3 \[Gamma]0^4+2 Sqrt[2] b m M^(3/2) \[Gamma]0^2 (M+M \[Gamma]0-\[Lambda])+2 m^2 (M+M \[Gamma]0-\[Lambda])^2]-(Sqrt[2] (2 m+\[Gamma]0 (Sqrt[2] b Sqrt[M]+2 a (M-\[Lambda]))) Erf[Abs[2 m+Sqrt[2] b Sqrt[M] \[Gamma]0+2 a M \[Gamma]0-2 a \[Gamma]0 \[Lambda]]/(2 Sqrt[2] Sqrt[m-a M \[Gamma]0^2])])/Abs[2 m+Sqrt[2] b Sqrt[M] \[Gamma]0+2 a M \[Gamma]0-2 a \[Gamma]0 \[Lambda]])-(E^((-b^2 M^2 \[Gamma]0^2+4 c M (-m+a M \[Gamma]0^2)+2 Sqrt[2] b m Sqrt[M] (M+M \[Gamma]0-\[Lambda])-2 a m (M+M \[Gamma]0-\[Lambda])^2)/(4 M (-m+a M \[Gamma]0^2))) Sqrt[m/(m-a M \[Gamma]0^2)] (((2 m (\[Gamma]0-\[Gamma]t)+\[Gamma]0^2 (-Sqrt[2] b Sqrt[M]+2 a (M+M \[Gamma]t-\[Lambda]))) Abs[Sqrt[2] b M^(3/2) \[Gamma]0^2-2 m (M+M \[Gamma]0-\[Lambda])] Erf[Abs[2 m (\[Gamma]0-\[Gamma]t)+\[Gamma]0^2 (-Sqrt[2] b Sqrt[M]+2 a (M+M \[Gamma]t-\[Lambda]))]/(2 Sqrt[2] \[Gamma]0 Sqrt[m-a M \[Gamma]0^2])])/(m-a M \[Gamma]0^2)+2 (Sqrt[2] b M^(3/2) \[Gamma]0^2-2 m (M+M \[Gamma]0-\[Lambda])) Abs[\[Gamma]t+(\[Gamma]0 (-2 m+\[Gamma]0 (Sqrt[2] b Sqrt[M]+2 a (-M+\[Lambda]))))/(2 (m-a M \[Gamma]0^2))] Erf[Abs[Sqrt[2] b Sqrt[M] \[Gamma]0^2+m (-2-2 \[Gamma]0+(2 \[Lambda])/M)]/(2 Sqrt[2] \[Gamma]0 Sqrt[m-a M \[Gamma]0^2])]))/(4 Abs[Sqrt[2] b M^(3/2) \[Gamma]0^2-2 m (M+M \[Gamma]0-\[Lambda])] Abs[\[Gamma]t+(\[Gamma]0 (-2 m+\[Gamma]0 (Sqrt[2] b Sqrt[M]+2 a (-M+\[Lambda]))))/(2 (m-a M \[Gamma]0^2))])+1/2 Erfc[-((Sqrt[m] (M+M \[Gamma]0-\[Lambda]))/(Sqrt[2] M \[Gamma]0))])/.LopezBenitezParameters[(-M (1+\[Gamma]0)+\[Lambda])/(2 Sqrt[M])],
				True,
					Undefined
			],
		diversityType == "SLC",
			(1/2 (1+Erf[(m (M n (1+\[Gamma]0)-\[Lambda]))/(Sqrt[2] M Sqrt[m n] \[Gamma]0)])+1/2 E^((4 c M n (-m+a M \[Gamma]0^2)+b (-b M^2 n \[Gamma]0^2+2 Sqrt[2] m Sqrt[M n] (M n (1+\[Gamma]0)-\[Lambda]))-2 a m (-M n (1+\[Gamma]0)+\[Lambda])^2)/(4 M n (-m+a M \[Gamma]0^2))) Sqrt[m/(m-a M \[Gamma]0^2)] (Erf[(b Sqrt[M^3 n] \[Gamma]0^2+Sqrt[2] m (-M n (1+\[Gamma]0)+\[Lambda]))/(2 M \[Gamma]0 Sqrt[n (m-a M \[Gamma]0^2)])]-Erf[(-2 m n+\[Gamma]0 (-2 a M n+Sqrt[2] b Sqrt[M n]+2 a \[Lambda]))/(2 Sqrt[2] Sqrt[n (m-a M \[Gamma]0^2)])])+1/2 E^(-((4 c M n (m-a M \[Gamma]0^2)+b (b M^2 n \[Gamma]0^2+2 Sqrt[2] m Sqrt[M n] (M n (1+\[Gamma]0)-\[Lambda]))+2 a m (-M n (1+\[Gamma]0)+\[Lambda])^2)/(4 M n (-m+a M \[Gamma]0^2)))) Sqrt[m/(m-a M \[Gamma]0^2)] (-2+Erfc[(b Sqrt[M^3 n] \[Gamma]0^2+Sqrt[2] m (M n (1+\[Gamma]0)-\[Lambda]))/(2 M \[Gamma]0 Sqrt[n (m-a M \[Gamma]0^2)])]))/.LopezBenitezParameters[(-M (n+\[Gamma]0)+\[Lambda])/(2 Sqrt[M n])],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - LopezBenitezNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - LopezBenitezNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate)*)


Options[NGaussianNakagamiProbabilityOfDetection] = {DiversityType -> OptionValue[NakagamiProbabilityOfDetection, DiversityType]};
NGaussianNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[NGaussianNakagamiProbabilityOfDetection, "NGaussian"];
NGaussianNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NGaussianNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[NGaussianNakagamiProbabilityOfDetection]]]
]
NGaussianNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, diversityType = OptionValue[DiversityType], method = "ApproximateNumerical"},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "MRC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "EGC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SEC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->{diversityType, \[Gamma]t}],{x,0,\[Infinity]}],
		diversityType == "SLC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - NGaussianNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NGaussianNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
				True,
					Undefined
			],
		True,
			Undefined
	]
]


(* ::Subsubsection::Closed:: *)
(*Numerical method (approximate, low SNR)*)


Options[NGaussianLowSNRNakagamiProbabilityOfDetection] = {DiversityType -> OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
NGaussianLowSNRNakagamiProbabilityOfDetection::usage = GenerateMethodHelp[NGaussianLowSNRNakagamiProbabilityOfDetection, "NGaussian"];
NGaussianLowSNRNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NGaussianLowSNRNakagamiProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];
	NGaussianLowSNRNakagamiProbabilityOfDetection[M,\[Gamma],\[Lambda],m,n,#/.(DiversityType/.#)->"None"&[RelevantOptions[NGaussianLowSNRNakagamiProbabilityOfDetection]]]
]
NGaussianLowSNRNakagamiProbabilityOfDetection[M_?NumericQ,\[Gamma]_,\[Lambda]_,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, diversityType = OptionValue[DiversityType], method = "ApproximateNumericalLowSNR"},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	Which[
		diversityType == "None",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "MRC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "EGC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SEC",
			NIntegrate[AWGNProbabilityOfDetection[M,x,\[Lambda],Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->{diversityType, \[Gamma]t}],{x,0,\[Infinity]}],
		diversityType == "SLC",
			NIntegrate[AWGNProbabilityOfDetection[M,x/n,\[Lambda],n,Method->method,DiversityType->diversityType]NakagamiPDF[\[Gamma]0,m,x,n,DiversityType->diversityType],{x,0,\[Infinity]}],
		diversityType == "SLS",
			Which[
				ListQ[\[Gamma]0],
					1 - Product[1 - NGaussianNakagamiProbabilityOfDetection[M,\[Gamma]0[[i]],\[Lambda],m,DiversityType->"None"],{i,n}],
				!ListQ[\[Gamma]0],
					1 - (1 - NGaussianNakagamiProbabilityOfDetection[M,\[Gamma]0,\[Lambda],m,DiversityType->"None"])^n,
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


Options[LowSNRAssumptionErrorNakagami] = {DiversityType->OptionValue[NakagamiProbabilityOfDetection,DiversityType]};
LowSNRAssumptionErrorNakagami::usage="LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m, n] calculates the upper bound for the low SNR approximation error.\n\n"<>DiversityTypeHelp[LowSNRAssumptionErrorNakagami];
LowSNRAssumptionErrorNakagami[M_,\[Gamma]_,\[Lambda]_,m_] := LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m] = Module[{n = 1},
	LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m, n, DiversityType->"None"]
]
LowSNRAssumptionErrorNakagami[M_,\[Gamma]_,\[Lambda]_,m_,n_,OptionsPattern[]] := Module[{diversityType = OptionValue[DiversityType], \[Gamma]t, g, \[Epsilon], \[Epsilon]max = 10},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];

	(*g[\[Epsilon]_?NumericQ] := Which[
		diversityType == "None" || diversityType == "SLC",
			Abs[AWGNProbabilityOfDetection[M, \[Epsilon], \[Lambda], n, DiversityType->diversityType, Method->"ApproximateNumerical"] - AWGNProbabilityOfDetection[M, \[Epsilon], \[Lambda], n, DiversityType->diversityType, Method->"ApproximateNumericalLowSNR"]],
		diversityType == "MRC" || diversityType == "EGC" || diversityType == "SC" || diversityType == "SEC",
			LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m],
		True,
			Undefined
	];

	If[diversityType == "SLS",
		(* This bound is VERY loose *)
		With[{Pm = (1 - ProbabilityOfDetection[M, \[Gamma], \[Lambda], ChannelType -> {"Nakagami", m}, Method -> "ApproximateNumerical", DiversityType -> "None"])},
			Min[Abs[{Pm^n - (Pm - LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m])^n, Pm^n - (Pm + LowSNRAssumptionErrorNakagami[M, \[Gamma], \[Lambda], m])^n}]]
		],
		NMaximize[{g[\[Epsilon]], 0 <= \[Epsilon] <= \[Epsilon]max}, {\[Epsilon], 0, \[Epsilon]max}][[1]]
	]*)

	(*Abs[Q[(\[Lambda] - M (1 + n \[Gamma]))/Sqrt[2 M (1 + 2 n \[Gamma])]] - Q[(\[Lambda] - M (1 +n \[Gamma]))/Sqrt[2 M]]]*)
	(*Abs[-(E^((-((\[Lambda]-M(1+n \[Gamma]))/Sqrt[2M(1+2 n \[Gamma])])^2/2))/Sqrt[2\[Pi]])((\[Lambda]-M(1+n \[Gamma]))/Sqrt[2M(1+2n \[Gamma])])(1-Sqrt[1+2n \[Gamma]])]//N*)
	(*Abs[E^(-((\[Lambda]-M(1+n \[Gamma]))/Sqrt[2M])^2/2)/Sqrt[2\[Pi]] \[Gamma]]//N*)
	(2^(-(1/2)+m n) E^(-((M-\[Lambda])^2/(4 M))) (2+2 Sqrt[2] Sqrt[M]+M)^(-((m n)/2)) (m/\[Gamma])^(m n) (-(((2 m Sqrt[M]+(Sqrt[2]+Sqrt[M]) \[Gamma] (M-\[Lambda])) Gamma[1+(m n)/2] Hypergeometric1F1[1+(m n)/2,3/2,(2 m Sqrt[M]+(Sqrt[2]+Sqrt[M]) \[Gamma] (M-\[Lambda]))^2/(4 M (2+2 Sqrt[2] Sqrt[M]+M) \[Gamma]^2)])/(Sqrt[M] (2+2 Sqrt[2] Sqrt[M]+M) \[Gamma]))+(Gamma[1/2 (1+m n)] Hypergeometric1F1[1/2 (1+m n),1/2,(2 m Sqrt[M]+(Sqrt[2]+Sqrt[M]) \[Gamma] (M-\[Lambda]))^2/(4 M (2+2 Sqrt[2] Sqrt[M]+M) \[Gamma]^2)])/Sqrt[2+2 Sqrt[2] Sqrt[M]+M]))/(Sqrt[\[Pi]] Gamma[m n])//N
]


(* ::Subsubsection::Closed:: *)
(*Approximaton error for Nakagami's PDF for EGC diversity*)


NakagamiEGCPDFApproximationError::usage="LowSNRAssumptionErrorNakagami[M, \[Lambda], n] calculates the upper bound for the low SNR approximation error.\n\n"<>DiversityTypeHelp[LowSNRAssumptionErrorNakagami];
NakagamiEGCPDFApproximationError[Pf_,\[Gamma]_,m_] := NakagamiEGCPDFApproximationError[Pf, \[Gamma], m] = Module[{n = 1},
	NakagamiEGCPDFApproximationError[Pf, \[Gamma], m, n]
]
NakagamiEGCPDFApproximationError[Pf_,\[Gamma]_,m_,n_] := Module[{diversityType = "EGC", \[Gamma]0, g, \[Epsilon], \[Epsilon]max = 10},
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	g[\[Epsilon]_?NumericQ] := (1 - Pf) Abs[Integrate[PDF[NakagamiEGCSumDistribution[\[Gamma], m , n], z], {z, 0, \[Epsilon]}] - Integrate[NakagamiPDF[\[Gamma], m, z, n, DiversityType->diversityType], {z, 0, \[Epsilon]}]];

	NMaximize[{g[\[Epsilon]], 0 <= \[Epsilon] <= \[Epsilon]max}, {\[Epsilon], 0, \[Epsilon]max}][[1]]
]


NakagamiEGCSumDistribution[\[Gamma]_,m_,n_] := NakagamiEGCSumDistribution[\[Gamma], m, n] = Module[{diversityType = "EGC", \[Gamma]0, x, pdf, numberOfPoints = 100000},
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	x = Table[Unique[], {n}];
	SmoothKernelDistribution[RandomVariate[TransformedDistribution[Sum[Sqrt[x[[i]]], {i, n}]^2 / n, Table[x[[i]] \[Distributed] GammaDistribution[m, \[Gamma]0 / m], {i, n}]], numberOfPoints]]
]


(* ::Subsubsection::Closed:: *)
(*Asymptotic approximation error*)


AsymptoticErrorNakagami::usage="AsymptoticErrorNakagami[Pf, m n] gives the upper bound for the error of the asymptotic method for the specified parameters.";
AsymptoticErrorNakagami[Pf_,mn_] := AsymptoticErrorNakagami[Pf, mn] = Module[{f, z},
	f[z_] := (Pf / 2) Erfc[Sqrt[mn / 2]] + (1 - Pf) (GammaRegularized[mn, z] - (1 / 2) Erfc[(z - mn)/Sqrt[2 mn]]);

	NMaximize[{Abs[f[z]], z >= 0}, {z, 0, mn}][[1]]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[NakagamiSampleComplexity] = {DiversityType->OptionValue[SampleComplexity,DiversityType], Method->OptionValue[SampleComplexity,Method]};
NakagamiSampleComplexity::usage="NakagamiSampleComplexity[\[Gamma], Pf, Pd, m, n] calculates the number of samples required for the specified decision probabilities and signal to noise ratio in a Nakagami-m channel.\n\n"<>MethodHelp[NakagamiSampleComplexity, {"\"ExactNumerical\"", "\"ApproximateNumerical\"", "\"ApproximateNumericalLowSNR\"", "\"ApproximateAsymptotic\"", "\"ApproximateAsymptoticUpperBound\"", "\"ApproximateAsymptoticLowerBound\""}]<>"\n\n"<>DiversityTypeHelp[NakagamiSampleComplexity];
NakagamiSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,m_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NakagamiSampleComplexity][[All,1]]],Options[target][[All,1]]];
	NakagamiSampleComplexity[\[Gamma], Pf, Pd, m, n, RelevantOptions[NakagamiSampleComplexity]]
]
NakagamiSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,m_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{\[Gamma]t, \[Gamma]0, RelevantOptions, diversityType = OptionValue[DiversityType], method, mn, g, M, \[Epsilon]},
	(* Handle both lists and scalar values for diversityType *)
	{diversityType, \[Gamma]t} = ProcessDiversityType[diversityType];
	
	(* Convert lists of SNR values to averages or maxima, depending on the specified diversity type *)
	\[Gamma]0 = ProcessSNR[\[Gamma], diversityType];

	{method, mn} = ProcessMethod[OptionValue[Method]];

	(* Check for invalid combinations of inputs *)
	If[diversityType == "None" && n > 1, Return[Undefined]];
	If[\[Gamma]0 == Undefined, Return[Undefined]];

	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[NakagamiSampleComplexity][[All,1]]],Options[target][[All,1]]];

	Quiet[Which[
		method == "ExactNumerical",
			f[x_?NumericQ] := NakagamiProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], m, n, RelevantOptions[NakagamiProbabilityOfDetection]/."ExactNumerical"->"ExactAnnamalai"];
			M/.FindRoot[f[M] == Pd, {M, NakagamiSampleComplexity[\[Gamma]0, Pf, Pd, m, n, Method->"ApproximateNumerical", DiversityType->diversityType], 1, \[Infinity]}],
		method == "ApproximateNumerical",
			f[x_?NumericQ] := NakagamiProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], m, n, RelevantOptions[NakagamiProbabilityOfDetection]];
			M/.FindRoot[f[M] == Pd, {M, AWGNSampleComplexity[\[Gamma]0, Pf, Pd, n, RelevantOptions[AWGNSampleComplexity]], 1, \[Infinity]}],
		method == "ApproximateNumericalLowSNR",
			f[x_?NumericQ] := NakagamiProbabilityOfDetection[x, \[Gamma]0, \[Lambda][x, Pf, n, RelevantOptions[\[Lambda]]], m, n, RelevantOptions[NakagamiProbabilityOfDetection]];
			M/.FindRoot[f[M] == Pd, {M, AWGNSampleComplexity[\[Gamma]0, Pf, Pd, n, RelevantOptions[AWGNSampleComplexity]], 1, \[Infinity]}],
		(* Keep two names for legacy reasons *)
		method == "ApproximateLargeSNR" || method == "ApproximateSmallPf",
			Which[
				diversityType == "None",
					2 (m InverseQ[Pf] / (\[Gamma]0 InverseGammaRegularized[m, (Pd - Pf) / (1 - Pf)]))^2,
				diversityType == "MRC",
					2 (m InverseQ[Pf] / (\[Gamma]0 InverseGammaRegularized[m n, (Pd - Pf) / (1 - Pf)]))^2,
				diversityType == "EGC",
					With[{\[Beta] = (Gamma[m + 1]^2 + m (n - 1) Gamma[m + 1 / 2]^2) / (n Gamma[m + 1]^2)},
						2 (m InverseQ[Pf] / (\[Beta] \[Gamma]0 InverseGammaRegularized[m n, (Pd - Pf) / (1 - Pf)]))^2
					],
				diversityType == "SLC",
					2 n (m InverseQ[Pf] / (\[Gamma]0 InverseGammaRegularized[m n, (Pd - Pf) / (1 - Pf)]))^2,
				True,
					Undefined
			],
		method == "ApproximateAsymptotic",
			Which[
				diversityType == "None",
					(2 m (-InverseQ[Pd]^4+m InverseQ[Pf]^2+InverseQ[Pd]^2 (m+InverseQ[Pf]^2)+2 Sqrt[m InverseQ[Pd]^2 InverseQ[Pf]^2 (m-InverseQ[Pd]^2+InverseQ[Pf]^2)]))/(\[Gamma]0^2 (m-InverseQ[Pd]^2)^2),
				diversityType == "MRC",
					(4 \[Sqrt](m^3 n^3 \[Gamma]0^4 InverseQ[Pd]^2 InverseQ[Pf]^2 (m n-InverseQ[Pd]^2+InverseQ[Pf]^2))+2 m n \[Gamma]0^2 (-InverseQ[Pd]^4+m n InverseQ[Pf]^2+InverseQ[Pd]^2 (m n+InverseQ[Pf]^2)))/(n^2 \[Gamma]0^4 (-m n+InverseQ[Pd]^2)^2),
				diversityType == "EGC",
					With[{c = Gamma[m + 1]^2 / (Gamma[m + 1]^2 + m (n - 1) Gamma[m + 1 / 2]^2)},
						(4 \[Sqrt](c^4 m^3 n^3 \[Gamma]0^4 InverseQ[Pd]^2 InverseQ[Pf]^2 (m n-InverseQ[Pd]^2+InverseQ[Pf]^2))+2 c^2 m n \[Gamma]0^2 (-InverseQ[Pd]^4+m n InverseQ[Pf]^2+InverseQ[Pd]^2 (m n+InverseQ[Pf]^2)))/(\[Gamma]0^4 (-m n+InverseQ[Pd]^2)^2)
					],
				diversityType == "SLC",
					(4 \[Sqrt](m^3 n InverseQ[Pd]^2 InverseQ[Pf]^2 (m n-InverseQ[Pd]^2+InverseQ[Pf]^2))+2 m (-InverseQ[Pd]^4+m n InverseQ[Pf]^2+InverseQ[Pd]^2 (m n+InverseQ[Pf]^2)))/(\[Gamma]0^2 (-m n+InverseQ[Pd]^2)^2),
				True,
					Undefined
			],
		method == "ApproximateAsymptoticUpperBound",
			\[Epsilon] = AsymptoticErrorNakagami[Pf, m n];
			If[Pd + \[Epsilon] >= 1, \[Infinity], NakagamiSampleComplexity[\[Gamma], Pf, Pd + \[Epsilon], m, n, Method->"ApproximateAsymptotic", DiversityType->diversityType]],
		method == "ApproximateAsymptoticLowerBound",
			\[Epsilon] = AsymptoticErrorNakagami[Pf, m n];
			If[Pd - \[Epsilon] <= 0, 0, NakagamiSampleComplexity[\[Gamma], Pf, Pd - \[Epsilon], m, n, Method->"ApproximateAsymptotic", DiversityType->diversityType]],
		True,
			Undefined
	]]
]


(* ::Subsection::Closed:: *)
(*Miscellaneous*)


MultinomialCoefficient::usage="Computes the multinomial coefficient for the specified inputs.";
MultinomialCoefficient[l_?NumericQ,k_?NumericQ,m_?NumericQ]:=Which[
	k == 0,
		1,
	k == 1,
		l,
	k == (m - 1) l,
		1 / (Gamma[m]^l),
	2 <= k <= (m - 1) l - 1,
		(1 / k) Total[Table[MultinomialCoefficient[l, k - j, m] (j (l + 1) - k) / j!, {j, 1, Min[k, m - 1]}]],
	True,
		Undefined
]


End[];


EndPackage[];

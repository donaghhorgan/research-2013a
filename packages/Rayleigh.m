(* ::Package:: *)

(* ::Title:: *)
(*Rayleigh channel functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for cooperative energy detection in Rayleigh channels.*)
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
(*13/12/2012*)
(*1.12*)


(* ::Subsection:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.12: Updated help functions and cleaned up code.*)
(*Version 1.11: Moved database logging functions to the Network package.*)
(*Version 1.1: Introduced RelevantOptions function and changed function definitions, so that child options are inherited from parents.*)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Rayleigh`"];


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


RayleighPDF;


(* ::Subsection::Closed:: *)
(*Detection probability*)


RayleighProbabilityOfDetection;


(* ::Subsection::Closed:: *)
(*Sample complexity*)


RayleighSampleComplexity;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<Nakagami`;


(* ::Subsection::Closed:: *)
(*PDF of the signal to noise ratio*)


Options[RayleighPDF] = Options[NakagamiPDF];
RayleighPDF::usage = StringReplace[NakagamiPDF::usage, {"Nakagami-m"->"Rayleigh", "Nakagami"->"Rayleigh", ", m"->""}];
RayleighPDF[\[Gamma]_,x_,OptionsPattern[]]:=Module[{n = 1}, RayleighPDF[\[Gamma], x, n, Method->OptionValue[Method]]]
RayleighPDF[\[Gamma]_,x_,n_,OptionsPattern[]]:=Module[{m = 1}, NakagamiPDF[\[Gamma], m, x, n, Method->OptionValue[Method]]]


(* ::Subsection::Closed:: *)
(*Detection probability*)


Options[RayleighProbabilityOfDetection] = Options[NakagamiProbabilityOfDetection];
RayleighProbabilityOfDetection::usage = StringReplace[NakagamiProbabilityOfDetection::usage, {"Nakagami-m"->"Rayleigh", "Nakagami"->"Rayleigh", ", m"->""}];
RayleighProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RayleighProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];

	RayleighProbabilityOfDetection[M, \[Gamma], \[Lambda], n, RelevantOptions[RayleighProbabilityOfDetection]]
]
RayleighProbabilityOfDetection[M_?NumericQ,\[Gamma]_?NumericQ,\[Lambda]_,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, m = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RayleighProbabilityOfDetection][[All,1]]],Options[target][[All,1]]];

	NakagamiProbabilityOfDetection[M, \[Gamma], \[Lambda], m, n, RelevantOptions[NakagamiProbabilityOfDetection]]
]


(* ::Subsection::Closed:: *)
(*Sample complexity*)


Options[RayleighSampleComplexity] = Options[NakagamiSampleComplexity];
RayleighSampleComplexity::usage = StringReplace[NakagamiSampleComplexity::usage, {"Nakagami-m"->"Rayleigh", "Nakagami"->"Rayleigh", ", m"->""}];
RayleighSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,OptionsPattern[]]:=Module[{RelevantOptions, n = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RayleighSampleComplexity][[All,1]]],Options[target][[All,1]]];

	RayleighSampleComplexity[\[Gamma], Pf, Pd, n, RelevantOptions[RayleighSampleComplexity]]
]
RayleighSampleComplexity[\[Gamma]_?NumericQ,Pf_?NumericQ,Pd_?NumericQ,n_?IntegerQ,OptionsPattern[]]:=Module[{RelevantOptions, m = 1},
	RelevantOptions[target_]:=FilterRules[Table[#[[i]]->OptionValue[#[[i]]],{i,Length[#]}]&[Options[RayleighSampleComplexity][[All,1]]],Options[target][[All,1]]];

	NakagamiSampleComplexity[\[Gamma], Pf, Pd, m, n, RelevantOptions[NakagamiSampleComplexity]]
]


End[];


EndPackage[];

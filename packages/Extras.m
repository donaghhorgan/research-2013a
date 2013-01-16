(* ::Package:: *)

(* ::Title:: *)
(*Extra functions*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions to support other packages.*)
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
(*1.05*)


(* ::Subsection:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.05: Added ProcessChannelType function.*)
(*Version 1.04: Added ProcessMethod function.*)
(*Version 1.03: Fixed a bug in ProcessSNR.*)
(*Version 1.02: Moved FaddeevaDerivative function from the Nakagami package.*)
(*Version 1.01: Added SEC support, and removed SSC support.*)
(*Version 1.0: First working version.*)


(* ::Section:: *)
(*Public*)


BeginPackage["Extras`"]; 


ProcessSNR;


ProcessMethod;


ProcessChannelType;


ProcessDiversityType;


FaddeevaDerivative;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


ProcessSNR::usage="ProcessSNR[\[Gamma], diversityType] processes lists of SNR values for the specified diversity type.";
ProcessSNR[\[Gamma]_,diversityType_:"SLC"]:=Which[
	diversityType == "None" || diversityType == "MRC" || diversityType == "EGC" || diversityType == "SEC" || diversityType == "SLC" || diversityType == "SC",
		If[ListQ[\[Gamma]], Undefined, \[Gamma]],
	diversityType == "SLS",
		\[Gamma],
	True,
		Undefined
]


ProcessMethod::usage="ProcessMethod[method] processes the specified method and returns a list.";
ProcessMethod[method_]:=If[ListQ[method], {method[[1]], method[[2]]}, {method, Null}, Undefined]


ProcessChannelType::usage="ProcessChannelType[channelType] processes the specified channel type and returns a list."
ProcessChannelType[channelType_]:=If[ListQ[channelType], {channelType[[1]], channelType[[2]]}, {channelType, Null}, Undefined]


ProcessDiversityType::usage="ProcessDiversityType[diversityType] processes the specified diversity type and returns a list.";
ProcessDiversityType[diversityType_]:=If[ListQ[diversityType], {diversityType[[1]], diversityType[[2]]}, {diversityType, Null}, Undefined]


FaddeevaDerivative::usage="Computes the \!\(\*SuperscriptBox[\(k\), \(th\)]\) derivative of the Faddeeva function w(z).";
FaddeevaDerivative[0, z_] := Exp[-z^2] Erfc[-I z];
FaddeevaDerivative[1, z_] := -2 z FaddeevaDerivative[0, z] + (2 I)/Sqrt[\[Pi]];
FaddeevaDerivative[k_?IntegerQ, z_] := FaddeevaDerivative[k, z] = -2 z FaddeevaDerivative[k - 1, z] - 2 (k - 1) FaddeevaDerivative[k - 2, z] // Simplify;


End[];


EndPackage[];

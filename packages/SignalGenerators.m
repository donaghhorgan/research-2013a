(* ::Package:: *)

(* ::Title:: *)
(*Signal generators*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica function definitions for different signal generator types.*)
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


(* ::Subsection::Closed:: *)
(*Version information*)


(* ::Text:: *)
(*28/08/2012*)
(*1.0*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.0: First working version, minor bug fixes to follow.*)


(* ::Section:: *)
(*Public*)


BeginPackage["SignalGenerators`"]; 


(* ::Subsection::Closed:: *)
(*Binary signal generator*)


BinarySignal;


(* ::Subsection::Closed:: *)
(*Noise signal generator*)


NoiseSignal;


(* ::Subsection::Closed:: *)
(*Modulators*)


Modulate;


(* ::Subsection::Closed:: *)
(*Demodulators*)


Demodulate;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Binary signal generator*)


BinarySignal::usage = "BinarySignal[n] generates a random binary signal of length n.
BinarySignal[n, type] generates a binary signal of the specified type of length n. Type can be either \"Random\" or \"Alternating\". By default, the type is \"Random\".";
BinarySignal[n_,type_:"Random"]:=Which[
	type=="Random",
	RandomInteger[1,n],
	type=="Alternating",
	Table[Mod[i,2],{i,n}],
	True,
	Undefined
]


(* ::Subsection::Closed:: *)
(*Noise signal generator*)


NoiseSignal::usage = "NoiseSignal[n, P] generates an AWGN signal of length n and power P.
NoiseSignal[n, P, type] generates a noise signal of the specified type of length n and power P. Type can be either \"AWGN\" or \"White\". By default, the type is \"AWGN\".";
NoiseSignal[n_,P_,type_:"AWGN"]:=Which[
	type=="AWGN",
	Sqrt[P] RandomReal[NormalDistribution[],n],
	type=="White",
	Sqrt[3P] RandomReal[{-1,1},n],
	True,
	Undefined
]


(* ::Subsection:: *)
(*Modulators*)


Options[Modulate]={ModulationScheme->"BPSK",SamplingFrequency->"Default"};
Modulate::usage="Modulate[signal] modulates the given binary signal using BPSK.

The following options may be given for ModulationScheme:

ModulationScheme->\"BPSK\".";
Modulate[signal_,power_,symbolPeriod_,frequency_,OptionsPattern[]]:=Module[{modulationScheme=OptionValue[ModulationScheme],Ts=symbolPeriod,fc=frequency,n=Length[signal],P=power,fs=If[OptionValue[SamplingFrequency]=="Default",2frequency,OptionValue[SamplingFrequency]]},
	Which[
		modulationScheme=="BPSK",
			Flatten[Table[Sqrt[2 P]Cos[2\[Pi] fc t + \[Pi](1 - i)],{i,signal},{t,Range[0, Ts - 1 / fs, 1 / fs]}]],		
		True,
			Undefined
	]
]


(* ::Subsection:: *)
(*Demodulators*)


Options[Demodulate]=Options[Modulate];
Demodulate::usage="Demodulate[signal] demodulates the given BPSK signal to binary.

The following options may be given for ModulationScheme:

ModulationScheme->\"BPSK\".";
Demodulate[signal_,power_,symbolPeriod_,frequency_,OptionsPattern[]]:=Module[{modulationScheme=OptionValue[ModulationScheme],Ts=symbolPeriod,fc=frequency,n=Length[signal],temp,P=power,fs=If[OptionValue[SamplingFrequency]=="Default",2frequency,OptionValue[SamplingFrequency]]},
	Which[
		modulationScheme=="BPSK",
		temp = signal / P * Table[Cos[2\[Pi] fc t],{t,Range[0, (n - 1) / fs, 1 / fs]}];
		temp = Table[Mean[temp[[(i - 1)Ts fs + 1 ;; i Ts fs]]],{i,n/(Ts fs)}];
		Table[If[i<0,0,1],{i,temp}],
		True,
		Undefined
	]
]


End[];


EndPackage[];

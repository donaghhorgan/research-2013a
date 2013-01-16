(* ::Package:: *)

(* ::Title:: *)
(*Q-function*)


(* ::Subsection::Closed:: *)
(*Copyright notice*)


(* ::Text:: *)
(*Mathematica implementation of the inverse Marcum Q-function.*)
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
(*07/01/2013*)
(*1.0*)


(* ::Subsection::Closed:: *)
(*Changelog*)


(* ::Text:: *)
(*Version 1.0: Basic implemenation of the inverse of the Marcum Q-function.*)


(* ::Section:: *)
(*Public*)


BeginPackage["InverseMarcumQ`"];


(* ::Subsection::Closed:: *)
(*Inverse Marcum Q-function*)


InverseMarcumQ;


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


<<QFunction`


(* ::Subsection::Closed:: *)
(*Inverse Marcum Q-function*)


InverseMarcumQ::usage="InverseMarcumQ[m, a, z] calculates the value of the inverse of the Marcum Q-function, \!\(\*SubscriptBox[\(Q\), \(m\)]\)(a, b). If \!\(\*SubscriptBox[\(Q\), \(m\)]\)(a, b) = z, then this function will return b = \!\(\*SuperscriptBox[SubscriptBox[\(Q\), \(m\)], \(-1\)]\)(a, z).";
InverseMarcumQ[m_?NumericQ, a_?NumericQ, z_?NumericQ]:=Module[{lowerBound, upperBound, b},
	Which[
		m == 1,
			(* From Digital Communications, 2nd Edition, by Simon and Alouini, Equation 4.56, we have the following bounds on b = InverseMarcumQ[1, a, z] *)
			lowerBound = -a - Sqrt[-2 Log[z]];
			upperBound = a + Sqrt[-2 Log[z]];
			b/.FindRoot[MarcumQ[m, a, b] == z, {b, (lowerBound + upperBound) / 2 + $MachineEpsilon,\[NonBreakingSpace]lowerBound, upperBound}],
		True,
			Undefined
	]
]


End[];


EndPackage[];

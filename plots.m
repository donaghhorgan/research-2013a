(* ::Package:: *)

(* ::Title:: *)
(*Plotting functions*)


(* ::Text:: *)
(*Some functions for consistent custom plots across notebooks.*)


<<PlotLegends`


(* ::Section:: *)
(*Standard options*)


plotOptions={Joined->True,PlotRange->Full,ImageSize->400,GridLines->Automatic,Frame->True};
histogramOptions={PlotRange->Full,ImageSize->400,GridLines->Automatic,ChartStyle->Gray,Frame->True};
plotLegendOptions={LegendShadow->None};
plotMarkers={"",Graphics[{Black,Rectangle[]},ImageSize->7],Graphics[{Black,Disk[]},ImageSize->7],"",Graphics[{Black,Line[{{0,0},{0,1},{1,1},{1,0},{0,0}}]},ImageSize->7],Graphics[{Black,Circle[]},ImageSize->7],"",Graphics[{Black,Polygon[{{0,0},{0.5,1},{1,0}}]},ImageSize->7],Graphics[{Black,Polygon[{{0,0},{0.5,-1},{1,0}}]},ImageSize->7],"",Graphics[{Black,Line[{{0,0},{0.5,1},{1,0},{0,0}}]},ImageSize->7],Graphics[{Black,Line[{{0,0},{0.5,-1},{1,0},{0,0}}]},ImageSize->7]};


(* ::Section:: *)
(*XY plot*)


CustomPlot[x_,y_,label_]:=Module[{a0,b0,c0,d0,e0},
	{a0,b0,c0,d0} = Dimensions[y];	
	ListLogPlot[Table[{x[[a]],y[[a,b,c,d]]},{d,d0},{c,c0},{b,b0},{a,a0}],plotOptions,FrameLabel->label]
]


(* ::Section:: *)
(*Error histogram*)


CustomHistogram[x_,label_]:=Histogram[Flatten[x],Automatic,"PDF",histogramOptions,FrameLabel->label]


(* ::Section:: *)
(*Statistics table*)


StatsTable[x_]:=TableForm[{{Mean[x//Flatten],StandardDeviation[x//Flatten]}},TableHeadings->{{},{"\!\(\*StyleBox[\"\[Mu]\", \"TR\"]\)","\[Sigma]"}}]


(* ::Section:: *)
(*Error table*)


MaxErrorsTable[x_,n_]:=Module[{a0,b0,c0,d0,y,p},
	{a0,b0,c0,d0} = Dimensions[x];
	y = Sort[Flatten[x],Greater][[1;;n]];
	p = Flatten[Table[Position[x,y[[i]]],{i,Length[y]}],1];
	TableForm[Flatten[{{{"\[Gamma] (dB)","n","m","M","pf","\[Epsilon]","\!\(\*SubscriptBox[\(\[Epsilon]\), \(r\)]\)"}},Table[{snrdbRange[[p[[i,1]]]],nRange[[p[[i,2]]]],If[mRange//ListQ,mRange[[p[[i,2]]]],1],samplesRange[[p[[i,3]]]],pfRange[[p[[i,4]]]],Extract[error,p[[i]]],Extract[relError,p[[i]]]}//N,{i,Length[p]}]},1]]
]


(* ::Section:: *)
(*ROC plot*)


ROCPlot[x_,legend_]:=ListPlot[x,PlotRangePadding->None,ImageSize->450,Joined->{True,False,False,True,False,False,True,False,False},PlotRange->{{0,1},{0,1}},PlotStyle->{Black,Black,Black,{Black,Dotted},Black,Black,{Black,Dashed},Black,Black,{Black,DotDashed},Black,Black},PlotLegend->legend,LegendShadow->None,LegendPosition->{1,-0.5},LegendSize->{1.3,1},LegendTextSpace->10,GridLines->Automatic,PlotMarkers->plotMarkers,Frame->True,FrameLabel->{"\!\(\*SubscriptBox[\(P\), \(f\)]\)","\!\(\*SubscriptBox[\(P\), \(d\)]\)"},AspectRatio->1]


(* ::Section:: *)
(*Optimum mn plot*)


MNPlot[x1_,x2_,legend_]:=Module[{f},
	f[y_,n_]:=Select[y,Round[#[[1]]]==n&]//Mean;
	ListLogLogPlot[{Table[f[{mnRange,Table[Mean[Flatten[Transpose[#][[i]]]]//Abs,{i,Length[Transpose[#]]}]&[x1]}//Transpose,mn],{mn,mnRange//Round//Union}],Table[f[{mnRange,Table[StandardDeviation[Flatten[Transpose[#][[i]]]]//Abs,{i,Length[Transpose[#]]}]&[x1]}//Transpose,mn],{mn,mnRange//Round//Union}],Table[f[{mnRange,Mean[Mean[x2]//Transpose]}//Transpose,mn],{mn,mnRange//Round//Union}],Table[f[{mnRange,StandardDeviation[StandardDeviation[x2]//Transpose]}//Transpose,mn],{mn,mnRange//Round//Union}]},Joined->True,PlotLegend->legend,plotLegendOptions,LegendPosition->{1,-0.35/GoldenRatio},LegendSize->{0.3,0.5},plotOptions,FrameLabel->{"mn","Error"},PlotStyle->{Black,Black,{Black,Dashed},{Black,Dashed}},PlotMarkers->plotMarkers]
]

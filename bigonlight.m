(* ::Package:: *)

(*************************************** Package info ******************************************)


(* :Title: BiGONLight *)

(* :Author: Michele Grasso (grasso@cft.edu.pl)*)

(* :Description: 
	BiGONLight (Bi-local geodesic operators framework for numerical light ptopagation) is
	a Mathematica package which encode the Bi-local Geodesic Operators formalism (BGO)
	to study light propagation in the geometric optics regime in General Relativit.
	More information about the BGO framework can be found in Phys.Rev.D 99,064038 (2019)
    (arXiv:1811.10284) or Phys.Rev.D 101,063506 (2020) (arXiv:1912.04988). 
    More information about the package can be found in Grasso, Villa (2021) and in 
    Grasso, Villa, Korzynski, Matarrese (2021). *)

(* :Discussion:
   - Perform the 3+1 splitting of a given 4D metric tensor Subscript[g, \[Mu]\[Nu]]
   - Define Christoffel symbols, covariant derivatives, Lie derivatives, Riemann tensor and 
        optical tidal matrix in terms of 3+1 quantities.
   - Define Gauss, Codazzi and Ricci relatins to compute the optical tidal matrix
   - Compute initial conditions for null tangent vectors
   - Compute the geodesic equation, parallel transport equations and geodesic deviation
            equation for the BGO.
   - Define bases and charts. Compute components. Package xCoba.
   - Define the concept of tensor equation and tensor rule.
*)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 10.0 and later *)


(* This program is free software; you can redistribute it and/or modify it under the terms of the
    GNU General Public License.
	Permission is hereby granted, free of charge, to any person obtaining a copy of this software
	and associated documentation files (the \[OpenCurlyDoubleQuote]Software\[CloseCurlyDoubleQuote]), to deal in the Software without
	restriction, including without limitation the rights to use, copy, modify, merge and to permit
	persons to whom the Software is furnished to do so, subject to the ollowing conditions:
	The above copyright notice and this permission notice shall be included in all copies or
    substantial portions of the Software.

	The Software is provided \[OpenCurlyDoubleQuote]as is\[CloseCurlyDoubleQuote], without warranty of any kind, expressed or implied, including
	but not limited to the warranties of merchantability, fitness for a particular purpose and 
	noninfringement. In no event shall the authors or copyright holders be liable for any claim,
    damages or other liability, whether in an action of contract, tort or otherwise, arising from,
    out of or in connection with the Software or the use or other dealings in the Software. 
*)


BeginPackage["bigonlight`"];



TransformedMetric::usage="TransformedMetric[\!\(\*SubscriptBox[\(g\), \(\[Mu]\[Nu]\)]\),\!\(\*SuperscriptBox[\(oldX\), \(\[Mu]\)]\), \!\(\*SuperscriptBox[\(newX\), \(\[Mu]\)]\), transformation]. Given the coordinate transformations, returns the full covariant metric in the new coordinates"
TransformedMetric::dimension="Trubles with matrix and coordinates dimensions.";
ADM::usage = "ADM[g\[Mu]\[Nu],\!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\)]. Returns the 3-metric \!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), the lapse \[Alpha] and the shift \!\(\*SuperscriptBox[\(\[Beta]\), \(i\)]\) given the full covariant metric g\[Mu]\[Nu]."
SNFtoADM::usage="SNFtoADM[\!\(\*SuperscriptBox[\(\[ScriptCapitalV]\), \(\[Mu]\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(e\), \(i\)], \((a)\)]\),\!\(\*SuperscriptBox[SubscriptBox[\(e\), \(0\)], \((a)\)]\),\!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\),\!\(\*SubscriptBox[\(n\), \(\[Mu]\)]\)]. Returns the 3+1 form of a vector \[ScriptCapitalV] given the coframe\!\(\*SuperscriptBox[\(\\\ \), \((4)\)]\)\!\(\*SubscriptBox[SuperscriptBox[\(e\), \(\[Mu]\)], \(A\)]\) and\!\(\*SuperscriptBox[\(\\\ \), \((4)\)]\)\!\(\*SuperscriptBox[\(\[ScriptCapitalV]\), \(A\)]\) the 4D vector in the frame "
Vsplit::usage="Vsplit[\!\(\*SuperscriptBox[\(\[ScriptCapitalV]\), \(\[Mu]\)]\),\!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\),\!\(\*SubscriptBox[\(n\), \(\[Mu]\)]\)]. Returns the 3+1 form of a four-vector as \[ScriptCapitalV]=vecN( \!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\) + \!\(\*SuperscriptBox[\(vecT\), \(\[Mu]\)]\)) given the normal vector \!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\) and \!\(\*SubscriptBox[\(n\), \(\[Mu]\)]\)]"
Christoffel::usage="Christoffel[\!\(\*SubscriptBox[\(g\), \(\[Mu]\[Nu]\)]\), \!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\)]. Returns the Christoffel simbols \!\(\*SubscriptBox[SuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\)], \(\[Mu]\[Nu]\)]\) for the given metric in the given coordinates."
ConD::usage="Returns the covariant derivative of fully contravariant tensor T with Chrisoffel symbols \[CapitalGamma] in coordinates 'chart'."
CovD::usage = "Returns the covariant derivative of fully covariant tensor \[Omega] with Chrisoffel symbols \[CapitalGamma] in coordinates 'chart'."
LieD::usage = "LieD[\[Omega],X,chart,g,\[CapitalGamma]] returns the covariant derivative of fully covariant tensor \[Omega] with Chrisoffel symbols \[CapitalGamma] in coordinates 'chart' with metric g along vector field X.
Limited to tensors of rank no higher than 2."
Riemann::usage = "Riemann[\!\(\*SubscriptBox[\(g\), \(\[Mu]\[Nu]\)]\),\!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\),\[CapitalGamma]] returns the Riemann tensor for the metric g in coordinates with Christoffel symbols \[CapitalGamma]."
GeodesicEquations::usage="GeodesicEquations[\!\(\*SuperscriptBox[\(X\), \(i\)]\),\!\(\*SuperscriptBox[\(V\), \(i\)]\),\!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), \!\(\*SubscriptBox[\(K\), \(ij\)]\), \[Alpha], \[Beta], time]. Returns the geodesic equation in 3+1 form. It contains the RHS for dx/dt and dV/dt"
EnergyEquations::usage="EnergyEquations[\!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), \!\(\*SubscriptBox[\(K\), \(ij\)]\), \[Alpha], \[Beta], {E,\[Lambda]},\!\(\*SuperscriptBox[\(X\), \(i\)]\),\!\(\*SuperscriptBox[\(V\), \(i\)]\), time]. Returns the equations for d\[Lambda]/dt and dE/dt"
InitialConditions::usage="InitialConditions[\!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), {(\!\(\*SuperscriptBox[\(V\), \(2\)]\)\!\(\*SubscriptBox[\()\), \(in\)]\),(\!\(\*SuperscriptBox[\(V\), \(3\)]\)\!\(\*SubscriptBox[\()\), \(in\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(X\), \(i\)], \(in\)]\)}, \!\(\*SuperscriptBox[\(V\), \(i\)]\), param]. Compute the initial value for the first spatial component \!\(\*SuperscriptBox[\(V\), \(1\)]\) using the spatial null condition \!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\)\!\(\*SuperscriptBox[\(V\), \(i\)]\)\!\(\*SuperscriptBox[\(V\), \(j\)]\)==1" 
ParallelTransport::usage="ParallelTransport[\!\(\*SuperscriptBox[\(X\), \(i\)]\),\!\(\*SuperscriptBox[\(V\), \(i\)]\), \!\(\*SuperscriptBox[\(e\), \(i\)]\), \!\(\*SuperscriptBox[\(e\), \(0\)]\), \!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), \!\(\*SubscriptBox[\(K\), \(ij\)]\), \[Alpha], \!\(\*SuperscriptBox[\(\[Beta]\), \(i\)]\), time]. Returns the parallel transport equations in 3+1 form for a vector \!\(\*SuperscriptBox[\(e\), \(\[Mu]\)]\)=\!\(\*SuperscriptBox[\(e\), \(0\)]\) \!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\) + {0, \!\(\*SuperscriptBox[\(e\), \(i\)]\)}"
OpticalTidalMatrix::usage="OpticalTidalMatrix[\!\(\*SuperscriptBox[\(X\), \(i\)]\),\!\(\*SuperscriptBox[\(V\), \(i\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(i\)], \((a)\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(0\)], \((a)\)]\), Q, \!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), \!\(\*SubscriptBox[\(K\), \(ij\)]\), \[Alpha], \!\(\*SuperscriptBox[\(\[Beta]\), \(i\)]\), time]. Returns the optical tidal matrix computed using only 3+1 quantities. The basic idea is to take the contraction of the full Riemann with the tangent vector and the frame vectors \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(\[Mu]\)], \((a)\)]\) splitted in 3+1. This will give the different parts of the Riemann 
decomposition along \!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\) and \!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\). Using the Riemann's symmetries and the expression for the non-vanishing projection (known as the Gauss, Codazzi e Ricci relations), we obtain the optical tidal matrix in the frame \!\(\*SubscriptBox[SuperscriptBox[\(R\), \((a)\)], \(kk \((b)\)\)]\)"
SolveGeodesic::usage="SolveGeodesic[geod_eq, initial_conditions, \!\(\*SuperscriptBox[\(X\), \(i\)]\),\!\(\*SuperscriptBox[\(V\), \(i\)]\),param, time, tin,tfin, method, workingprecision_,precision_,interpolation_, number_maxsteps]. Solve the ODE for the geodesic equations and returns the position and the tangent vectors"
SolveEnergy::usage="SolveEnergy[energy_eq, initial_conditions, {E,\[Lambda]}, geodesic_solution, param, time, tin,tfin, method, workingprecision_,precision_,interpolation_, number_maxsteps]. Solve the ODE for the energy"
SolvePTransport::usage="SolvePTransport[parallelTransport_eq, geodesic_solution, initial_conditions, \!\(\*SuperscriptBox[\(e\), \(i\)]\), \!\(\*SuperscriptBox[\(e\), \(0\)]\), param, time, tin,tfin, method, workingprecision_,precision_,interpolation_, number_maxsteps]. Solve the parallel transport equations in 3+1 form, and returns the value of C and the components"
PTransportedFrame::usage="PTransportedFrame[\!\(\*SubscriptBox[\(\[Gamma]\), \(ij\)]\), \!\(\*SubscriptBox[\(K\), \(ij\)]\), \[Alpha], \!\(\*SuperscriptBox[\(\[Beta]\), \(i\)]\), geod_sol, frame_initialCond, \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(i\)], \((a)\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(0\)], \((a)\)]\), param, time, tin,tfin, method, workingprecision_,precision_,interpolation_, number_maxsteps]. Compute and solve the parallel transport equation for the semi-null frame \!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(\[Mu]\)], \((a)\)]\)={\!\(\*SuperscriptBox[\(u\), \(\[Mu]\)]\),\!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(\[Mu]\)], \((1)\)]\),\!\(\*SubscriptBox[SuperscriptBox[\(\[Phi]\), \(\[Mu]\)], \((2)\)]\),\!\(\*SuperscriptBox[\(k\), \(\[Mu]\)]\)} in 3+1 form."
SolveBGO::usage="SolveBGO[]. Solve the geodesic deviation equations in the parallel transported frame for a group of the Bilocal Geodesic Operators."
GDE::usage="GDE[geodesic_sol, \!\(\*SuperscriptBox[\(\[Xi]\), \((a)\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(R\), \((a)\)], \(kk \((b)\)\)]\), (\!\(\*SuperscriptBox[\(\[Xi]\), \((a)\)]\)\!\(\*SubscriptBox[\()\), \(in\)]\), \[Alpha], E, param, time, ti,tf, method, workingprecision_,precision_,interpolation_, number_maxsteps]. Compute and solve the geodesic deviation equation in the parallel transported frame along the geodesic. It can accept a vector \!\(\*SuperscriptBox[\(\[Xi]\), \((a)\)]\) or a matrix \!\(\*SubsuperscriptBox[\(M\), \(\(\\\ \\\ \)\((b)\)\), \((a)\)]\).
WORKS WELL SO FAR, BUT IT NEEDS MORE TESTING. HOWEVER, IT IS BETTER TO USE BGOEQUATIONS[]!"
BGOequations::usage="BGOequations[\!\(\*SubscriptBox[\(W\), \(x\)]\),\!\(\*SubscriptBox[\(W\), \(y\)]\), OPT[t], \[Alpha][t], E[t], param, t]. Compute the geodesic deviation equations in the parallel transported frame for a group of the Bilocal Geodesic Operators {\!\(\*SubscriptBox[\(\[ScriptCapitalW]\), \(XL\)]\),\!\(\*SubscriptBox[\(\[ScriptCapitalW]\), \(LL\)]\)} and {\!\(\*SubscriptBox[\(\[ScriptCapitalW]\), \(XX\)]\),\!\(\*SubscriptBox[\(\[ScriptCapitalW]\), \(LX\)]\)} in the form {\!\(\*SubscriptBox[\(W\), \(x\)]\)'[t]= \!\(\*FractionBox[\(\[Alpha]\), \(E\)]\) \!\(\*SubscriptBox[\(W\), \(y\)]\)[t], \!\(\*SubscriptBox[\(W\), \(y\)]\)'[t]= \!\(\*FractionBox[\(\[Alpha]\), \(E\)]\) OPT[t].\!\(\*SubscriptBox[\(W\), \(x\)]\)[t] "
CodazziRelation::usage="CodazziRelation[\!\(\*SubscriptBox[\(K\), \(\[Mu]\[Nu]\)]\), \[CapitalGamma], \!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\)]. Returns the Codazzi-Minardi relation as a tensor \!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\[Mu]\[Alpha]\[Beta]\)]\) . NOTE: the 3+1 quantities has to be given as 4 dimensional where time components are zero"
GaussRelation::usage="GaussRelation[\!\(\*SubscriptBox[\(K\), \(\[Mu]\[Nu]\)]\), \!\(\*SubscriptBox[SuperscriptBox[\(R\), \(\[Mu]\)], \(\[Alpha]\[Beta]\[Nu]\)]\), \!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\[Nu]\)]\), \!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\)]. Returns the Gauss-Codazzi relation as a tensor \!\(\*SubscriptBox[\(\[ScriptCapitalG]\), \(\[Mu]\[Alpha]\[Beta]\[Nu]\)]\) . NOTE: the 3+1 quantities has to be given as 4 dimensional where time components are zero"
RicciRelation::usage="RicciRelation[\!\(\*SubscriptBox[\(K\), \(\[Mu]\[Nu]\)]\),\!\(\*SubscriptBox[SuperscriptBox[\(K\), \(\[Mu]\)], \(\[Nu]\)]\),\[CapitalGamma],\[Alpha], \!\(\*SuperscriptBox[\(n\), \(\[Mu]\)]\),\!\(\*SuperscriptBox[\(X\), \(\[Mu]\)]\)]. Returns the Ricci relation as a tensor \!\(\*SubscriptBox[\(\[ScriptCapitalR]\), \(\[Beta]\[Alpha]\)]\) . NOTE: the 3+1 quantities has to be given as 4 dimensional where time components are zero"
CheckRiemann::usage="CheckRiemann[coordinates_, metric_, extrinsic_, lapse_, shift_, t_] uses the Gauss Subscript[\[ScriptCapitalG], \[Rho]\[Mu]\[Nu]\[Sigma]], Codazzi Subscript[\[ScriptCapitalC], \[Rho]\[Mu]\[Nu]] and Ricci Subscript[\[ScriptCapitalR], \[Mu]\[Nu]] relations to calculate the 4D Riemann Subscript[R, \[Rho]\[Mu]\[Nu]\[Sigma]] = Subscript[\[ScriptCapitalG], \[Rho]\[Mu]\[Nu]\[Sigma]]-(Subscript[\[ScriptCapitalC], \[Rho]\[Mu]\[Nu]]Subscript[n, \[Sigma]]-Subscript[\[ScriptCapitalC], \[Rho]\[Mu]\[Sigma]]Subscript[n, \[Nu]])-(Subscript[\[ScriptCapitalC], \[Nu]\[Sigma]\[Rho]]Subscript[n, \[Mu]]-Subscript[\[ScriptCapitalC], \[Nu]\[Sigma]\[Mu]]Subscript[n, \[Rho]]) + (Subscript[\[ScriptCapitalR], \[Rho]\[Nu]]Subscript[n, \[Sigma]]Subscript[n, \[Mu]]-Subscript[\[ScriptCapitalR], \[Rho]\[Sigma]]Subscript[n, \[Nu]]Subscript[n, \[Mu]])-(Subscript[\[ScriptCapitalR], \[Mu]\[Nu]]Subscript[n, \[Sigma]]Subscript[n, \[Rho]]-Subscript[\[ScriptCapitalR], \[Mu]\[Sigma]]Subscript[n, \[Nu]]Subscript[n, \[Rho]]). It also returns Subscript[\[ScriptCapitalG], \[Rho]\[Mu]\[Nu]\[Sigma]], Subscript[\[ScriptCapitalC], \[Rho]\[Mu]\[Nu]] and Subscript[\[ScriptCapitalR], \[Mu]\[Nu]] accessible through the output components (i.e. CheckRiemann[[1]]==Subscript[R, \[Rho]\[Mu]\[Nu]\[Sigma]], CheckRiemann[[2]]==Subscript[\[ScriptCapitalG], \[Rho]\[Mu]\[Nu]\[Sigma]], CheckRiemann[[3]]==Subscript[\[ScriptCapitalC], \[Rho]\[Mu]\[Nu]], CheckRiemann[[4]]==Subscript[\[ScriptCapitalR], \[Mu]\[Nu]])"


Begin["`Private`"];


TransformedMetric[metric_, oldcoordinates_, newcoordinates_, transformation_]:=Module[{met=metric, oldX=oldcoordinates, newX=newcoordinates, transf=transformation,dim=Length@oldcoordinates, olddX, jacob},
olddX= oldX//.transf;
jacob=Table[D[olddX[[indi]],newX[[indj]]],{indi,1,dim},{indj,1,dim}];
If[!MatchQ[Dimensions[met], Dimensions[jacob]], Message[TransformedMetric::dimension]];
Transpose[jacob].met.jacob//Simplify
]



ADM[metric_, coordinates_]:=Module[{dim=Length@coordinates, met=metric, coord=coordinates, Alpha, \[Alpha]\[Alpha], \[Beta]\[Beta], \[Beta]\[Beta]1, \[Beta]\[Beta]2, \[Beta]\[Beta]3, \[Gamma]\[Gamma],\[CapitalGamma]\[CapitalGamma], metuu,Kext,Nup,Ndown,Dndd},
\[CapitalGamma]\[CapitalGamma]=Christoffel[met,coord];
\[Gamma]\[Gamma]=Table[met[[indi,indj]],{indi,2,dim},{indj,2,dim}];
metuu=Inverse[met];
\[Alpha]\[Alpha]= -( 1/(metuu[[1,1]]));
Alpha=Sqrt[(\[Alpha]\[Alpha])]//Simplify;
\[Beta]\[Beta]={\[Alpha]\[Alpha] metuu[[1,2]],\[Alpha]\[Alpha] metuu[[1,3]],\[Alpha]\[Alpha] metuu[[1,4]]}//Simplify;
Nup=Join[{1/Alpha},-(\[Beta]\[Beta]/Alpha)];
Ndown={-Alpha,0,0,0};
Dndd=Table[D[Ndown[[indi]],coord[[indj]]]-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(3\)]\(\[CapitalGamma]\[CapitalGamma][\([indk, indi, indj]\)] Ndown[\([indk]\)]\)\),{indi,1,dim},{indj,1,dim}];
Kext=Table[-(Dndd[[indi,indj]]+Ndown[[indi]](\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(dim\)]\(Nup[\([indk]\)] Dndd[\([indk, indj]\)]\)\))),{indi,1,dim},{indj,1,dim}];
{Alpha,\[Beta]\[Beta],\[Gamma]\[Gamma],Take[Kext,-(dim-1),-(dim-1)],Nup,Ndown,Kext}
]


SNFtoADM[vectorinframe_,coframe_,coframe0_,normup_,normdown_]:=Module[{vectA=vectorinframe,EE=coframe,EE0=coframe0,nu=normup,nd=normdown,EE4,vect4,vectN,vectT},
EE4=Table[Flatten[{0,Table[EE[[ind\[Mu],indA]],{ind\[Mu],1,3}]}],{indA,1,4}];
vect4=Table[(Sum[vectA[[indA]](EE0[[indA]]nu[[ind\[Mu]]]+EE4[[ind\[Mu],indA]]),{indA,1,4}]),{ind\[Mu],1,4}];
vectN=Sum[-vect4[[indi]]nd[[indi]],{indi,1,4}];
vectT=Table[vect4[[ind\[Mu]]]-vectN,{ind\[Mu],1,4}];
{vectN,Take[vectT,-3],vect4}
]


Vsplit[vector_,normup_,normdown_]:=Module[{vect=vector,nu=normup,nd=normdown,vectN,vectT},
vectN=Sum[-vect[[indi]]nd[[indi]],{indi,1,4}];
vectT=Table[vect[[indi]]/vectN-nu[[indi]],{indi,1,4}];
{vectN,Take[vectT,-3]}
]


Christoffel[metric_,coordinates_]:=Block[{met=metric, coord=coordinates, n=Length@coordinates}, 
Table[(1/2 \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Sigma] = 1\), \(n\)]\(\(Inverse[met]\)[\([ind\[Lambda], ind\[Sigma]]\)] \((D[met[\([ind\[Sigma], ind\[Nu]]\)], coord[\([ind\[Mu]]\)]] + D[met[\([ind\[Sigma], ind\[Mu]]\)], coord[\([ind\[Nu]]\)]] - D[met[\([ind\[Mu], ind\[Nu]]\)], coord[\([ind\[Sigma]]\)]])\)\)\))//Simplify,{ind\[Lambda],1,n},{ind\[Mu],1,n},{ind\[Nu],1,n}]
]


CovD[\[Tau]_,coordinates_,christoffel_]:=
Block[{x=coordinates,n=Length@coordinates,\[CapitalGamma]=christoffel,rk=Assuming[\[Tau]\[Element]Reals,TensorRank[\[Tau]]],ind,a},
ind=Table[Subscript[a, i],{i,1,rk}];
Transpose[Table[D[\[Tau][[Sequence@@ind]],x[[der]]]-Plus@@Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Sigma] = 1\), \(n\)]\((\[CapitalGamma][\([\[Sigma], der, ind[\([j]\)]]\)] \[Tau][\([Sequence @@ ReplacePart[ind, j -> \[Sigma]]]\)])\)\),{j,1,rk}],Evaluate@(Sequence@@Table[{ind[[i]],1,n},{i,1,rk}]),{der,1,n}]]
]


ConD[T_,coordinates_,christoffel_]:=Block[{x=coordinates,n=Length@coordinates,\[CapitalGamma]=christoffel,rk=Assuming[T\[Element]Reals,TensorRank[T]],ind,a},
ind=Table[Subscript[a, i],{i,1,rk}];
Transpose[Table[D[T[[Sequence@@ind]],x[[der]]]+Plus@@Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Sigma] = 1\), \(n\)]\((\[CapitalGamma][\([\[Sigma], der, ind[\([j]\)]]\)] T[\([Sequence @@ ReplacePart[ind, j -> \[Sigma]]]\)])\)\),{j,1,rk}],Evaluate@(Sequence@@Table[{ind[[i]],1,n},{i,1,rk}]),{der,1,n}]]
]



LieD[\[Omega]_,X_,coordinates_,metric_,christoffel___]:=
Module[{g=metric,x=coordinates,\[CapitalGamma]={christoffel},n=Length@coordinates,rk=Assuming[\[Omega]\[Element]Reals,TensorRank[\[Omega]]],ind,a},
If[\[CapitalGamma]=={},\[CapitalGamma]=Christoffel[g,x],\[CapitalGamma]=First@\[CapitalGamma]];
ind=Table[Subscript[a, i],{i,1,rk}];
If[
rk==0,\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Sigma] = 1\), \(n\)]\(X[\([\[Sigma]]\)] Transpose[\(CovD[\[Omega], x, \[CapitalGamma]]\)[\([\[Sigma]]\)]]\)\),
Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Sigma] = 1\), \(n\)]\((X[\([\[Sigma]]\)] \(Transpose[CovD[\[Omega], x, \[CapitalGamma]]]\)[\([Sequence @@ Join[ind, {\[Sigma]}]]\)] + Plus @@ Table[
\*SubsuperscriptBox[\(\[Sum]\), \(\[Tau] = 1\), \(n\)]\((\(Inverse[g]\)[\([\[Sigma], \[Tau]]\)] \[Omega][\([Sequence @@ ReplacePart[ind, j -> \[Tau]]]\)] \(Transpose[CovD[g . X, x, \[CapitalGamma]]]\)[\([\[Sigma], ind[\([j]\)]]\)])\), {j, 1, rk}])\)\),Evaluate@(Sequence@@Table[{ind[[i]],1,n},{i,1,rk}])]
]
]


Riemann[metric_,coordinates_,christoffel___]:=Module[{g=metric,x=coordinates,n=Length@coordinates,\[CapitalGamma]={christoffel}},
If[\[CapitalGamma]=={},\[CapitalGamma]=Christoffel[g,x],\[CapitalGamma]=First@\[CapitalGamma]];
Table[(D[\[CapitalGamma][[ind\[Lambda],ind\[Mu],ind\[Sigma]]],x[[ind\[Nu]]]]-D[\[CapitalGamma][[ind\[Lambda],ind\[Mu],ind\[Nu]]],x[[ind\[Sigma]]]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Eta] = 1\), \(n\)]\(\[CapitalGamma][\([ind\[Eta], ind\[Mu], ind\[Sigma]]\)] \[CapitalGamma][\([ind\[Lambda], ind\[Eta], ind\[Nu]]\)]\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Eta] = 1\), \(n\)]\(\[CapitalGamma][\([ind\[Eta], ind\[Mu], ind\[Nu]]\)] \[CapitalGamma][\([ind\[Lambda], ind\[Eta], ind\[Sigma]]\)]\)\)),{ind\[Lambda],1,n},{ind\[Mu],1,n},{ind\[Nu],1,n},{ind\[Sigma],1,n}]
]


GeodesicEquations[coordinates_, velocity_, metric_, extrinsic_, lapse_, shift_, t_]:=
Module[{coord=coordinates, vel=velocity, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, met=metric, ext=extrinsic, n=Length@coordinates, coordt, \[CapitalGamma], \[Gamma]uu, Kud},
coordt=Flatten[Table[{coord[[indi]]->coord[[indi]][t], vel[[indi]]->vel[[indi]][t] },{indi,1,n}]];
\[CapitalGamma]=Christoffel[met,coord];
\[Gamma]uu=Inverse[met];
Kud=Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\(\(Inverse[met]\)[\([indi, indk]\)] ext[\([indk, indj]\)]\)\),{indi,1,n},{indj,1,n}];
Flatten[Table[{D[coord[[indi]][t],t]== Evaluate[(\[Alpha]\[Alpha] vel[[indi]]-\[Beta]\[Beta][[indi]])/.coordt], D[vel[[indi]][t],t]== Evaluate[(\[Alpha]\[Alpha]((\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((vel[\([indj]\)] D[Log[\[Alpha]\[Alpha]], coord[\([indj]\)]])\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((ext[\([indj, indk]\)] vel[\([indj]\)] vel[\([indk]\)])\)\)\))vel[[indi]]+ \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((2\ Kud[\([indi, indj]\)] vel[\([indj]\)])\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((\[CapitalGamma][\([indi, indj, indk]\)] vel[\([indj]\)] vel[\([indk]\)])\)\)\))-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((\[Gamma]uu[\([indi, indj]\)] D[\[Alpha]\[Alpha], coord[\([indj]\)]])\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((vel[\([indj]\)] D[\[Beta]\[Beta][\([indi]\)], coord[\([indj]\)]])\)\))/.coordt]},{indi,1,n}]]
]




EnergyEquations[metric_,extrinsic_,lapse_,shift_,energyvariable_,coordinates_,velocity_,t_]:=
Module[{met=metric, ext=extrinsic, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, coord=coordinates, vel=velocity,n=Length@coordinates,  energy=energyvariable,RHS,coordt},
coordt=Flatten[Table[{coord[[indi]]->coord[[indi]][t], vel[[indi]]->vel[[indi]][t]},{indi,1,n}]];
RHS[t]=energy[[1]][t](Evaluate[(\[Alpha]\[Alpha] \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((ext[\([indj, indk]\)] vel[\([indj]\)] vel[\([indk]\)])\)\)\) - \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((vel[\([indj]\)] D[\[Alpha]\[Alpha], coord[\([indj]\)]])\)\))/.coordt]);
{D[energy[[1]][t],t]== Evaluate[RHS[t]], D[energy[[2]][t],t]==(Evaluate[\[Alpha]\[Alpha]/.coordt]/energy[[1]][t])}
]


InitialConditions[metric_, initialconditions_, velocity_,parameters_]:=Module[{met=metric, initial=initialconditions,constants=parameters, vel=velocity, n=Length@velocity, totcond},
totcond=Flatten[Join[initial,constants]];
Flatten[Solve[Evaluate[(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indi = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((met[\([indi, indj]\)] vel[\([indi]\)] vel[\([indj]\)])\)\)\))/.totcond]==1, vel[[1]]]]
]



ParallelTransport[coordinates_, velocity_,vector_,vectortime_, metric_, extrinsic_, lapse_, shift_, t_]:=
Module[{coord=coordinates, vel=velocity, vect=vector, vect0=vectortime, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, met=metric, ext=extrinsic, n=Length@coordinates, vectt, \[CapitalGamma], \[Gamma]uu, Kud},
vectt=Join[Flatten[Table[{coord[[indi]]->coord[[indi]][t], vel[[indi]]->vel[[indi]][t]},{indi,1,n}]],Table[vect[[indi]]->vect[[indi]][t],{indi,1,n}]];
\[CapitalGamma]=Christoffel[met,coord];
\[Gamma]uu=Inverse[met];
Kud=Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\(\(Inverse[met]\)[\([indi, indk]\)] ext[\([indk, indj]\)]\)\),{indi,1,n},{indj,1,n}];
Flatten[Join[{D[vect0[t],t]== Evaluate[(\[Alpha]\[Alpha] \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((\ ext[\([indj, indk]\)]\ vel[\([indj]\)]\ vect[\([indk]\)])\)\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((vect[\([indj]\)]\ D[\[Alpha]\[Alpha], coord[\([indj]\)]])\)\))/.vectt]},Table[ D[vect[[indi]][t],t]== Evaluate[(vect0[t](\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((\[Alpha]\[Alpha]\ Kud[\([indi, indk]\)] vel[\([indk]\)] - \[Gamma]uu[\([indi, indk]\)] D[\[Alpha]\[Alpha], coord[\([indk]\)]])\)\))+ \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((\[Alpha]\[Alpha]\ Kud[\([indi, indk]\)] vect[\([indk]\)])\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indj = 1\), \(n\)]\((\[Alpha]\[Alpha]\ \[CapitalGamma][\([indi, indj, indk]\)] vel[\([indj]\)]\ vect[\([indk]\)])\)\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\((vect[\([indk]\)] D[\[Beta]\[Beta][\([indi]\)], \ coord[\([indk]\)]])\)\))/.vectt],{indi,1,n}]]]
]



OpticalTidalMatrix[coordinates_, velocity_,vector_,vectortime_,UtimesK_, metric_, extrinsic_, lapse_, shift_, t_]:=
Module[{coord=coordinates, vel=velocity, f3=vector, \[Phi]=vectortime,QQ=UtimesK, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, met=metric, ext=extrinsic, n=(Length@coordinates+1),vectt, \[Gamma]uu, Ric,Cod,Gaus,\[Beta]\[Beta]3,\[CapitalGamma]\[CapitalGamma]3,\[Gamma]\[Gamma]3,I\[Gamma]\[Gamma]3,KK3,udKK3,KKud,RR3,x4,vel4,f,uno,due,S00,S01,S02,S11,S12,S22},
x4=Join[{t},coord];
vel4=Join[{0},vel];
f=Table[Join[{0},f3[[indi,All]]],{indi,1,4}];
vectt=Join[Flatten[Table[{coord[[indi]]->coord[[indi]][t], vel[[indi]]->vel[[indi]][t]},{indi,1,n-1}]],Flatten[Table[f3[[indi,indj]]->f3[[indi,indj]][t],{indi,1,Dimensions[f3][[1]]},{indj,1,Dimensions[f3][[2]]}]],Table[\[Phi][[indi]]->\[Phi][[indi]][t],{indi,1,Dimensions[\[Phi]][[1]]}]];
\[Beta]\[Beta]3=Join[{0},\[Beta]\[Beta]];
\[CapitalGamma]\[CapitalGamma]3=Table[If[inda==1||indb==1||indc==1,0,Christoffel[met,coord][[inda-1,indb-1,indc-1]]],{inda,1,4},{indb,1,4},{indc,1,4}];
\[Gamma]\[Gamma]3=Table[If[inda==1||indb==1,0,met[[inda-1,indb-1]]],{inda,1,4},{indb,1,4}];
\[Gamma]uu=Table[If[inda==1||indb==1,0,Inverse[met][[inda-1,indb-1]]],{inda,1,4},{indb,1,4}];
KK3=Table[If[inda==1||indb==1,0,ext[[inda-1,indb-1]]],{inda,1,4},{indb,1,4}];
KKud=Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(3\)]\(\(Inverse[met]\)[\([indi, indk]\)] ext[\([indk, indj]\)]\)\),{indi,1,3},{indj,1,3}];
udKK3=Table[If[inda==1||indb==1,0,KKud[[inda-1,indb-1]]],{inda,1,4},{indb,1,4}];
RR3=Table[If[inda==1||indb==1||indc==1||indd==1,0,Riemann[met,coord,Christoffel[met,coord]][[inda-1,indb-1,indc-1,indd-1]]],{inda,1,4},{indb,1,4},{indc,1,4},{indd,1,4}];
Ric=RicciRelation[KK3,udKK3,\[CapitalGamma]\[CapitalGamma]3,\[Alpha]\[Alpha],Join[{1/\[Alpha]\[Alpha]},-(\[Beta]\[Beta]/\[Alpha]\[Alpha])],x4];
Cod=CodazziRelation[KK3,\[CapitalGamma]\[CapitalGamma]3,x4];
Gaus=GaussRelation[KK3,RR3,\[Gamma]\[Gamma]3,x4];
uno=Table[Simplify[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indalpha = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indbeta = 1\), \(n\)]\((Ric[\([indbeta, indalpha]\)]\ \((\[Phi][\([indaa]\)] f[\([indbb, indbeta]\)] vel4[\([indalpha]\)] + \[Phi][\([indbb]\)] f[\([indaa, indbeta]\)] vel4[\([indalpha]\)] - \[Phi][\([indaa]\)] \[Phi][\([indbb]\)] vel4[\([indalpha]\)] vel4[\([indbeta]\)] - f[\([indbb, indbeta]\)] f[\([indaa, indalpha]\)])\))\)\)\)],{indaa,1,n-1},{indbb,1,n-1}];
due=Table[Simplify[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indnnu = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indalpha = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(indbeta = 1\), \(n\)]\((Cod[\([indnnu, indbeta, indalpha]\)] \((\[Phi][\([indaa]\)] f[\([indbb, indnnu]\)] vel4[\([indalpha]\)] vel4[\([indbeta]\)] + \[Phi][\([indbb]\)] f[\([indaa, indnnu]\)] vel4[\([indalpha]\)] vel4[\([indbeta]\)] - vel4[\([indbeta]\)] f[\([indbb, indnnu]\)] f[\([indaa, indalpha]\)] - vel4[\([indbeta]\)] f[\([indaa, indnnu]\)] f[\([indbb, indalpha]\)])\))\)\)\)\)],{indaa,1,n-1},{indbb,1,n-1}];
(*Table[(Simplify[(uno[[aa,bb]]+due[[aa,bb]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(\[Delta] = 1\), \(n\)]\((f[\([aa, \[Lambda]]\)]Gaus[\([\[Lambda], \[Rho], \[Epsilon], \[Delta]]\)]vel4[\([\[Rho]]\)]vel4[\([\[Epsilon]]\)]f[\([bb, \[Delta]]\)])\)\)\)\)\))])/.vectt,{aa,1,n},{bb,1,n}]*)
S00=Evaluate[((uno[[1,1]]+due[[1,1]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([1, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([1, ind\[Delta]]\)])\)\)\)\)\)))];
S01=Evaluate[((uno[[1,2]]+due[[1,2]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([1, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([2, ind\[Delta]]\)])\)\)\)\)\)))];
S02=Evaluate[((uno[[1,3]]+due[[1,3]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([1, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([3, ind\[Delta]]\)])\)\)\)\)\)))];
S11=Evaluate[((uno[[2,2]]+due[[2,2]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([2, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([2, ind\[Delta]]\)])\)\)\)\)\)))];
S12=Evaluate[((uno[[2,3]]+due[[2,3]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([2, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([3, ind\[Delta]]\)])\)\)\)\)\)))];
S22=Evaluate[((uno[[3,3]]+due[[3,3]]+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Lambda] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Rho] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Epsilon] = 1\), \(n\)]\(
\*SubsuperscriptBox[\(\[Sum]\), \(ind\[Delta] = 1\), \(n\)]\((f[\([3, ind\[Lambda]]\)] Gaus[\([ind\[Lambda], ind\[Rho], ind\[Epsilon], ind\[Delta]]\)] vel4[\([ind\[Rho]]\)] vel4[\([ind\[Epsilon]]\)] f[\([3, ind\[Delta]]\)])\)\)\)\)\)))];
{{0,0,0,0},{S01,S11,S12,0},{S02,S12,S22,0},{S00/QQ,S01/QQ,S02/QQ,0}}/.vectt
]


SolveGeodesic[geodesic_, initialconditions_, coordinates_, velocity_, parameters_,time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,numbersteps_]:=
Module[{coord=coordinates, vel=velocity, n=Length@coordinates, trajectory=geodesic, initial=initialconditions, constants=parameters,t=time,ti=tinitial,tf=tfinal,M=method,met,sol, wp=workingprecision, p=precision, int=interpolation, nstep=numbersteps},
If[M=={},Print["Please, specify a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic"]];
Which[M=="RK",met={"ExplicitRungeKutta","DifferenceOrder"->4} ,M=="SS",met={ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},M=="A",met={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],met={"Automatic"}}];
sol=NDSolve[Join[trajectory,initial]/.parameters,Join[coord,vel],{t,ti,tf},WorkingPrecision->wp,Method->met, PrecisionGoal->p, InterpolationOrder->int ,MaxSteps->10^nstep];
sol]


SolveEnergy[energyeq_,initialconditions_,energyvariable_,geodesic_, parameters_,time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,numbersteps_]:=
Module[{equations=energyeq, energy0=initialconditions,  constants=parameters,trajectory=geodesic, energy=energyvariable,t=time,ti=tinitial,tf=tfinal,M=method,met,sol,coordt, wp=workingprecision, p=precision, int=interpolation, nstep=numbersteps},
If[M=={},Print["Please, specify a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic"]];
Which[M=="RK",met={"ExplicitRungeKutta","DifferenceOrder"->4} ,M=="SS",met={ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},M=="A",met={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],met={"Automatic"}}];
sol=NDSolve[Join[equations,energy0]/.Join[trajectory,constants], energy,{t,ti, tf},WorkingPrecision->wp,Method->met, PrecisionGoal->p, InterpolationOrder->int,MaxSteps->10^nstep];
sol
]


SolvePTransport[paralleltransportequations_, geodesic_, initialconditions_,vector_,vectortime_, parameters_,time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,numbersteps_]:=
Module[{vect=vector,vect0=vectortime, partransequ=paralleltransportequations, trajectory=geodesic, initial=initialconditions, constants=parameters, t=time,ti=tinitial,tf=tfinal,M=method,wp=workingprecision, p=precision, int=interpolation,met,sol, nstep=numbersteps},
If[M=={},Print["Please, specify a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic"]];
Which[M=="RK",met={"ExplicitRungeKutta","DifferenceOrder"->4} ,M=="SS",met={ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},M=="A",met={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],met={"Automatic"}}];
sol=NDSolve[Join[partransequ,initial]/.Flatten[Join[constants,trajectory]],AppendTo[vect,vect0],{t,ti,tf},WorkingPrecision->wp,Method->met, PrecisionGoal->p, InterpolationOrder->int,MaxSteps->10^nstep];
sol
]


(*PTransportedFrame[paralleltransportequations_, geodesic_, initialconditions_,vector_,vectortime_, parameters_,time_,tinitial_,tfinal_,precision_,numbersteps_]:=
Module[{vect=vector,vect0=vectortime, partransequ=paralleltransportequations, trajectory=geodesic, initial=initialconditions, constants=parameters, t=time,ti=tinitial,tf=tfinal,sol,newvec,n=Length@vectortime,p=precision,nstep=numbersteps},
newvec=Table[Flatten[{e0[t]->vect0[[i]][t],e0'[t]->vect0[[i]]'[t],e1[t]->vect[[i,1]][t],e2[t]->vect[[i,2]][t],e3[t]->vect[[i,3]][t],e1'[t]->vect[[i,1]]'[t],e2'[t]->vect[[i,2]]'[t],e3'[t]->vect[[i,3]]'[t]}],{i,1,Length[vect0]}];
Table[Flatten[NDSolve[Join[Evaluate[partransequ/.newvec[[i]]],initial[[i]]]/.Flatten[Join[constants,trajectory]],Flatten[Join[vect[[i]],{vect0[[i]]}]],{t,ti,tf},WorkingPrecision->p,Method->{ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},MaxSteps->10^nstep ]],{i,1,n}]
]*)
PTransportedFrame[met_,extrinsic_,lapse_,shift_, geodesic_, initialconditions_,vector_,vectortime_, parameters_,time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,numbersteps_]:=
Module[{g=met,ext=extrinsic, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, vect=vector,vect0=vectortime, trajectory=geodesic, initial=initialconditions, constants=parameters, t=time,ti=tinitial,tf=tfinal,M=method,meth,wp=workingprecision, p=precision, int=interpolation,equations,solution,ii,coord,velocity,n=Length@vectortime,nstep=numbersteps},
solution={0,0,0,0};
coord=Table[geodesic[[in,1]],{in,1,3}];
velocity=Table[geodesic[[in,1]],{in,4,6}];
Which[M=="RK",meth={"ExplicitRungeKutta","DifferenceOrder"->4},M=="SS",meth={ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},M=="A",meth={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],meth={"Automatic"}}];
For[ii=1,ii<=4,ii++,equations=ParallelTransport[coord,velocity,vect[[ii]],vect0[[ii]],g,ext,\[Alpha]\[Alpha],\[Beta]\[Beta],t];
solution[[ii]]=NDSolve[Join[equations,initial[[ii]]]/.Flatten[Join[constants,trajectory]],Flatten[Join[vect[[ii]],{vect0[[ii]]}]],{t,ti,tf},WorkingPrecision->wp,Method->meth, PrecisionGoal->p, InterpolationOrder->int,MaxSteps->10^nstep]
];
solution
]



SolveBGO[equations_, initialconditions_, BGO1_, BGO2_, frame_, parameters_,time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,numbersteps_,timeconst_]:=
Module[{W1=BGO1, W2=BGO2, fr=frame, BGOeq=equations, initial=initialconditions, constants=parameters,t=time,ti=tinitial,tf=tfinal,M=method,met,sol,constraint=timeconst, wp=workingprecision, p=precision, int=interpolation, nstep=numbersteps},
If[M=={},Print["Please, specify a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic"]];
Which[M=="RK",met={"ExplicitRungeKutta","DifferenceOrder"->4} ,M=="SS",met={"StiffnessSwitching"},M=="A",met={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],met={"Automatic"}}];
With[{opts=SystemOptions[]},Internal`WithLocalSettings[SetSystemOptions["NDSolveOptions"->"DefaultSolveTimeConstraint"->constraint], sol=Flatten[NDSolve[Flatten[Join[BGOeq,initial]]/.Flatten[Join[constants,fr]],Join[Flatten[W1],Flatten[W2]],{t,ti,tf},WorkingPrecision->wp,Method->met, PrecisionGoal->p, InterpolationOrder->int ,MaxSteps->10^nstep]];,SetSystemOptions[opts]]];
sol]


GDE[geodesic_,deviation_,opticaltidalmatrix_, initialconditions_, lapse_, energy_,parameters_, time_,tinitial_,tfinal_,method_,workingprecision_,precision_,interpolation_,steps_]:=
Module[{ OPT=opticaltidalmatrix, \[Alpha]\[Alpha]=lapse,dev=deviation, n=Length@deviation,trajectory=geodesic, ener=energy, initial=initialconditions, constants=parameters, t=time,ti=tinitial,tf=tfinal,M=method,meth,sol,equations,vectdt,wp=workingprecision, p=precision, int=interpolation,step=steps},
vectdt=Table[trajectory[[All,1]][[i]]->trajectory[[All,1]][[i]][t],{i,1,Length[trajectory[[All,1]]]}];
If[TensorRank[deviation]==1,
equations=Table[(dev[[A]]''[t]==-dev[[A]]'[t]((ener'[t]/(ener[t]))-((\[Alpha]\[Alpha]'[t])/(\[Alpha]\[Alpha][t])))+((\[Alpha]\[Alpha][t]/ener[t])^2)(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(B = 1\), \(n\)]\((OPT[\([A, B]\)] \(dev[\([B]\)]\)[t])\)\))),{A,1,n}],equations=Table[(dev[[A,B]]''[t]==-dev[[A,B]]'[t]((ener'[t]/(ener[t]))-((\[Alpha]\[Alpha]'[t])/(\[Alpha]\[Alpha][t])))+((\[Alpha]\[Alpha][t])^2)(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(C = 1\), \(n\)]\((OPT[\([A, C]\)] \(dev[\([C, B]\)]\)[t])\)\))),{A,1,n},{B,1,Length[dev[[1]]]}]];
If[M=={},Print["Please, specify a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic"]];
Which[M=="RK",meth={"ExplicitRungeKutta","DifferenceOrder"->4},M=="SS",meth={ "StiffnessSwitching", Method->{"ExplicitRungeKutta",Automatic}},M=="A",meth={"Automatic"},M!="RK"&&M!="SS"&&M!="A", {Print["Invalid method specified. Please, choose a numerical method: RK= ExplicitRungeKutta, SS=StiffnessSwitching, A=Automatic. The defoult is Automatic"],meth={"Automatic"}}];
sol=NDSolve[Join[equations,initial]/.Flatten[Join[constants,trajectory]],Join[Flatten[dev]],{t,ti,tf},WorkingPrecision->wp,Method->meth, PrecisionGoal->p, InterpolationOrder->int,MaxSteps->10^step ];
sol
]




BGOequations[deviation_,Ddeviation_,opticaltidalmatrix_,lapse_,energy_,parameters_, time_]:=
Module[{ OPT=opticaltidalmatrix, W1=deviation,W2=Ddeviation,\[Alpha]\[Alpha]=lapse,ener=energy, n=Length@deviation, t=time,equations,vectdt},

Flatten[Table[{W1[[inda,indb]]'[t]==\[Alpha]\[Alpha]/ener W2[[inda,indb]][t],W2[[inda,indb]]'[t]==(\[Alpha]\[Alpha] /ener)(\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indc = 1\), \(n\)]\((OPT[\([inda, indc]\)] \(W1[\([indc, indb]\)]\)[t])\)\))},{inda,1,n},{indb,1,Length[W1[[1]]]}]]
]


CodazziRelation[extrinsic_,chritstoffel_,coordinates_]:=Module[{ext=extrinsic,\[CapitalGamma]=chritstoffel,x=coordinates,n=Length@coordinates, coddd},
coddd=Table[(D[ext[[ind1,ind3]],x[[ind2]]]-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\(\[CapitalGamma][\([sum1, ind1, ind2]\)] ext[\([sum1, ind3]\)]\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\(\[CapitalGamma][\([sum1, ind3, ind2]\)] ext[\([ind1, sum1]\)]\)\) -(D[ext[[ind2,ind3]],x[[ind1]]]-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\(\[CapitalGamma][\([sum1, ind2, ind1]\)] ext[\([sum1, ind3]\)]\)\)-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\(\[CapitalGamma][\([sum1, ind3, ind1]\)] ext[\([ind2, sum1]\)]\)\) ) )//Simplify,{ind1,1,n},{ind2,1,n},{ind3,1,n}];
Table[If[ind1==1||ind2==1||ind3==1,0,coddd[[ind1,ind2,ind3]]],{ind1,1,4},{ind2,1,4},{ind3,1,4}]]


GaussRelation[extrinsic_,riemann_,metric_,coordinates_]:=Module[{ext=extrinsic,R=riemann,coord=coordinates,met=metric,n=Length@coordinates,Rd},
Rd=Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indi = 1\), \(n\)]\(met[\([inda, indi]\)] R[\([indi, indb, indc, indd]\)]\)\),{inda,1,n},{indb,1,n},{indc,1,n},{indd,1,n}];
Table[(Rd[[inda,indb,indc,indd]]+ext[[inda,indc]]ext[[indb,indd]]-ext[[inda,indd]]ext[[indc,indb]] ),{inda,1,n},{indb,1,n},{indc,1,n},{indd,1,n}]
]



RicciRelation[extrinsic_,extrinsicud_,chritstoffel_,lapse_,nu_,coordinates_]:=Block[{coord=coordinates,nup=nu,n=Length@coordinates,\[CapitalGamma]=chritstoffel,\[Alpha]\[Alpha]=lapse,ext=extrinsic, Kud=extrinsicud,D\[Alpha]\[Alpha], DD\[Alpha]\[Alpha]},
D\[Alpha]\[Alpha]=Table[D[\[Alpha]\[Alpha],coord[[ind1]]],{ind1,1,n}];
DD\[Alpha]\[Alpha]=Table[Simplify[D[D\[Alpha]\[Alpha][[ind2]],coord[[ind1]]]-\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\(\[CapitalGamma][\([sum1, ind2, ind1]\)] D\[Alpha]\[Alpha][\([sum1]\)]\)\)],{ind1,1,n},{ind2,1,n}];
Table[Simplify[((\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\((nup[\([sum1]\)] D[ext[\([ind1, ind2]\)], coord[\([sum1]\)]])\)\)+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\((ext[\([ind1, sum1]\)] D[nup[\([sum1]\)], coord[\([ind2]\)]])\)\)+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\((ext[\([ind2, sum1]\)] D[nup[\([sum1]\)], coord[\([ind1]\)]])\)\))+\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(sum1 = 1\), \(n\)]\((ext[\([ind1, sum1]\)] Kud[\([sum1, ind2]\)])\)\)+DD\[Alpha]\[Alpha][[ind1,ind2]]/\[Alpha]\[Alpha] )],{ind1,1,n},{ind2,1,n}]
]




CheckRiemann[coordinates_, metric_, extrinsic_, lapse_, shift_, t_]:=
Module[{coord=coordinates, \[Alpha]\[Alpha]=lapse, \[Beta]\[Beta]=shift, met=metric, ext=extrinsic, n=Length@coordinates, \[Gamma]uu, Ric,Cod,Gaus,\[Beta]\[Beta]3,\[CapitalGamma]\[CapitalGamma]3,\[Gamma]\[Gamma]3,I\[Gamma]\[Gamma]3,KK3,udKK3,KKud,RR3,x4,dim4,normd},
x4=Join[{t},coord];
dim4=Length@x4;
normd={-\[Alpha]\[Alpha],0,0,0};
\[Beta]\[Beta]3=Join[{0},\[Beta]\[Beta]];
\[CapitalGamma]\[CapitalGamma]3=Table[If[inda==1||indb==1||indc==1,0,Christoffel[met,coord][[inda-1,indb-1,indc-1]]],{inda,1,dim4},{indb,1,dim4},{indc,1,dim4}];
\[Gamma]\[Gamma]3=Table[If[inda==1||indb==1,0,met[[inda-1,indb-1]]],{inda,1,dim4},{indb,1,dim4}];
\[Gamma]uu=Table[If[inda==1||indb==1,0,Inverse[met][[inda-1,indb-1]]],{inda,1,dim4},{indb,1,dim4}];
KK3=Table[If[inda==1||indb==1,0,ext[[inda-1,indb-1]]],{inda,1,dim4},{indb,1,dim4}];
KKud=Table[\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(indk = 1\), \(n\)]\(\(Inverse[met]\)[\([indi, indk]\)] ext[\([indk, indj]\)]\)\),{indi,1,n},{indj,1,n}];
udKK3=Table[If[inda==1||indb==1,0,KKud[[inda-1,indb-1]]],{inda,1,dim4},{indb,1,dim4}];
RR3=Table[If[inda==1||indb==1||indc==1||indd==1,0,Riemann[met,coord,Christoffel[met,coord]][[inda-1,indb-1,indc-1,indd-1]]],{inda,1,dim4},{indb,1,dim4},{indc,1,dim4},{indd,1,dim4}];
Ric=RicciRelation[KK3,udKK3,\[CapitalGamma]\[CapitalGamma]3,\[Alpha]\[Alpha],Join[{1/\[Alpha]\[Alpha]},-(\[Beta]\[Beta]/\[Alpha]\[Alpha])],x4];
Cod=CodazziRelation[KK3,\[CapitalGamma]\[CapitalGamma]3,x4];
Gaus=GaussRelation[KK3,RR3,\[Gamma]\[Gamma]3,x4];
{Table[(Gaus[[ia,ib,ic,id]]-(Cod[[ia,ib,ic]]normd[[id]]-Cod[[ia,ib,id]]normd[[ic]])-(Cod[[ic,id,ia]]normd[[ib]]-Cod[[ic,id,ib]]normd[[ia]])+(Ric[[ia,ic]]normd[[id]]normd[[ib]]-Ric[[ia,id]]normd[[ic]]normd[[ib]])-(Ric[[ib,ic]]normd[[id]]normd[[ia]]-Ric[[ib,id]]normd[[ic]]normd[[ia]]))//Simplify,{ia,1,dim4},{ib,1,dim4},{ic,1,dim4},{id,1,dim4}],Gaus, Cod, Ric}


]


End[];
EndPackage[];

(* ::Package:: *)

(* ::Section:: *)
(*Initialization*)


<< FeynCalc`


SetDirectory @ DirectoryName @ If[$FrontEnd === Null, $InputFileName, NotebookFileName[]];


(* ::Section:: *)
(*Vertices*)


iVh3Expr = Get @ FileNameJoin[Directory[], "data", "iVh3.m"];
VertexHHH[p1_, \[Mu]_, \[Nu]_, p2_, \[Rho]_, \[Sigma]_, p3_, \[Alpha]_, \[Beta]_] := iVh3Expr /. {
    ph1 -> p1, ph2 -> p2, ph3 -> p3,
    a -> \[Mu], b -> \[Nu], c -> \[Rho], d -> \[Sigma], e -> \[Alpha], f -> \[Beta]
}


iV\[Phi]\[Phi]hExpr = Get @ FileNameJoin[Directory[], "data", "iV\[Phi]\[Phi]h.m"];
VertexSSH[p1_, p2_, m_, \[Mu]_, \[Nu]_] := iV\[Phi]\[Phi]hExpr /. {
    p\[Phi]1 -> p1, p\[Phi]2 -> p2, a -> \[Mu], b -> \[Nu], m\[Phi] -> m
}


iV\[Phi]\[Phi]hhExpr = Get @ FileNameJoin[Directory[], "data", "iV\[Phi]\[Phi]hh.m"];
VertexSSHH[p1_, p2_, m_, \[Mu]_, \[Nu]_, \[Rho]_, \[Sigma]_] := iV\[Phi]\[Phi]hhExpr /. {
    p\[Phi]1 -> p1, p\[Phi]2 -> p2,
    a -> \[Mu], b -> \[Nu],
    c -> \[Rho], d -> \[Sigma],
    m\[Phi] -> m
}


iV\[Phi]\[CurlyPhi]\[CurlyPhi]Expr = Get @ FileNameJoin[Directory[], "data", "iV\[Phi]\[CurlyPhi]\[CurlyPhi].m"];
VertexS1S2S2[] := iV\[Phi]\[CurlyPhi]\[CurlyPhi]Expr


iV\[Phi]\[CurlyPhi]\[CurlyPhi]hExpr = Get @ FileNameJoin[Directory[], "data", "iV\[Phi]\[CurlyPhi]\[CurlyPhi]h.m"];
VertexS1S2S2H[\[Mu]_, \[Nu]_] := iV\[Phi]\[CurlyPhi]\[CurlyPhi]hExpr /. {a -> \[Mu], b -> \[Nu]}


(* ::Section:: *)
(*Propagators*)


PropagatorGraviton[p_, \[Mu]_, \[Nu]_, \[Rho]_, \[Sigma]_] := I (
    MT[\[Mu], \[Sigma]] MT[\[Nu], \[Rho]] + MT[\[Mu], \[Rho]] MT[\[Nu], \[Sigma]]-MT[\[Mu], \[Nu]] MT[\[Rho], \[Sigma]]
) / (2 SP[p])


PropagatorScalar[k_, m_] := I / (SP[k] - m^2)


(* ::Section:: *)
(*Polarization Summation*)


tildeMT[\[Mu]_, \[Nu]_, p_, r_] := MT[\[Mu], \[Nu]] -
    (FV[p, \[Mu]] FV[r, \[Nu]] + FV[p, \[Nu]] FV[r, \[Mu]]) / SP[p, r]
SumPol[\[Mu]_, \[Nu]_, \[Alpha]_, \[Beta]_, p_, r_] := (
	tildeMT[\[Mu], \[Alpha], p, r] tildeMT[\[Nu], \[Beta], p, r] +
	tildeMT[\[Mu], \[Beta], p, r] tildeMT[\[Nu], \[Alpha], p, r] -
	tildeMT[\[Mu], \[Nu], p, r] tildeMT[\[Alpha], \[Beta], p, r]
) / 2

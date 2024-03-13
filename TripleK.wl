(* ::Package:: *)

(* :Name: TripleK` *)


(* :Title: Tools for manipulation and evaluation of triple-K integrals and conformal correlation functions. *)


(* :Author: Adam Bzowski *)


(* :Summary:
     This package provides basic tools for manipulation and evaluation of triple-K integrals and conformal correlation functions.
     The package delivers tools for evaluation of triple-K integrals with integral and half-integral indices as well as divergences of arbitrary triple-K integrals. 
*)


(* :Context: TripleK` *)


(* :Package Version: 1.0.1 *)


(* :History:
     Version 1.0.1 by Adam Bzowski, June 2022.
*)


(* :Copyright: GNU General Public License v3.0, Adam Bzowski, 2020 *)


(* :Keywords:
     Triple-K, 
     Conformal field theory, 
     Feynman diagrams, 
     Loop integrals, 
     Dimensional regularization
*)


(* :Mathematica Version: 11.2 *)


BeginPackage["TripleK`", "Global`"];
Unprotect @@ Names["TripleK`*"];
ClearAll @@ Names["TripleK`*"];


(* ::Section:: *)
(*Interface*)


(* ::Text:: *)
(*Definitions.*)


p::usage = "\!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], 
StyleBox[\"j\",\nFontSlant->\"Italic\"], \(\[Mu]\)]\) for \!\(\*
StyleBox[\"j\",\nFontSlant->\"Italic\"]\)=1,2 represent external momenta.
\!\(\*
StyleBox[SubscriptBox[\"p\", \"j\"],\nFontSlant->\"Italic\"]\) for \!\(\*
StyleBox[\"j\",\nFontSlant->\"Italic\"]\)=1,2,3 represent magnitudes of vectors \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\), \(\[Mu]\)]\), \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\), \(\[Mu]\)]\), and \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(3\), \(\[Mu]\)]\)=-\!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\), \(\[Mu]\)]\)-\!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\), \(\[Mu]\)]\).";
\[Delta]::usage = "\!\(\*SubscriptBox[\(\[Delta]\), \(\[Mu]\[Nu]\)]\) represents the Euclidean metric.";
d::usage = "\!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\) represents the number of spacetime dimensions.";
\[Epsilon]::usage = "\[Epsilon] represents the regulator.";


LoopIntegral::usage = "LoopIntegral[\!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\), {\!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\),\!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\),\!\(\*SubscriptBox[\(\[Delta]\), \(3\)]\)}][\!\(\*
StyleBox[\"numerator\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)] represents the 3-point 1-loop integral over momentum \!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\) in \!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\) dimensions with denominator \!\(\*SuperscriptBox[
StyleBox[\"k\",\nFontSlant->\"Italic\"], \(2 \*SubscriptBox[\(\[Delta]\), \(3\)]\)]\)(\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)-\!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\)]\)\!\(\*SuperscriptBox[\()\), \(2 \*SubscriptBox[\(\[Delta]\), \(2\)]\)]\)(\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)+\!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\)]\)\!\(\*SuperscriptBox[\()\), \(2 \*SubscriptBox[\(\[Delta]\), \(1\)]\)]\) and the given \!\(\*
StyleBox[\"numerator\",\nFontSlant->\"Italic\"]\). External momenta are \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\), \(\[Mu]\)]\), \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\), \(\[Mu]\)]\), and \!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(3\), \(\[Mu]\)]\)=-\!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\), \(\[Mu]\)]\)-\!\(\*SubsuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\), \(\[Mu]\)]\).
LoopIntegral[\!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\), {\!\(\*SubscriptBox[\(\[Delta]\), \(1\)]\),\!\(\*SubscriptBox[\(\[Delta]\), \(2\)]\)}][\!\(\*
StyleBox[\"numerator\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)] represents the 2-point 1-loop integral over momentum \!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\) in \!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\) dimensions with denominator \!\(\*SuperscriptBox[
StyleBox[\"k\",\nFontSlant->\"Italic\"], \(2 \*SubscriptBox[\(\[Delta]\), \(1\)]\)]\)(\!\(\*
StyleBox[\"k\",\nFontSlant->\"Italic\"]\)-\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)\!\(\*SuperscriptBox[\()\), \(2 \*SubscriptBox[\(\[Delta]\), \(2\)]\)]\) and the given \!\(\*
StyleBox[\"numerator\",\nFontSlant->\"Italic\"]\). External momenta are \!\(\*SuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(\[Mu]\)]\) and -\!\(\*SuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(\[Mu]\)]\).";
i::usage = "\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)[\!\(\*
StyleBox[\"\[Alpha]\",\nFontSlant->\"Italic\"]\),{\!\(\*SubscriptBox[
StyleBox[\"\[Beta]\",\nFontSlant->\"Italic\"], \(1\)]\), \!\(\*SubscriptBox[
StyleBox[\"\[Beta]\",\nFontSlant->\"Italic\"], \(2\)]\)}][\!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)] represents the double-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integral with parameters \[Alpha], \!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\), \!\(\*SubscriptBox[\(\[Beta]\), \(2\)]\) and momentum \!\(\*SuperscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(\[Mu]\)]\).
\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)[\!\(\*
StyleBox[\"\[Alpha]\",\nFontSlant->\"Italic\"]\),{\!\(\*SubscriptBox[
StyleBox[\"\[Beta]\",\nFontSlant->\"Italic\"], \(1\)]\), \!\(\*SubscriptBox[
StyleBox[\"\[Beta]\",\nFontSlant->\"Italic\"], \(2\)]\), \!\(\*SubscriptBox[
StyleBox[\"\[Beta]\",\nFontSlant->\"Italic\"], \(3\)]\)}][\!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\)]\),\!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\)]\),\!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(3\)]\)] represents the triple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integral with parameters \[Alpha], \!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\), \!\(\*SubscriptBox[\(\[Beta]\), \(2\)]\), \!\(\*SubscriptBox[\(\[Beta]\), \(3\)]\) and momentum magnitudes \!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(1\)]\), \!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(2\)]\), \!\(\*SubscriptBox[
StyleBox[\"p\",\nFontSlant->\"Italic\"], \(3\)]\).";
NL::usage = "\!\(\*
StyleBox[\"NL\",\nFontSlant->\"Italic\"]\) represents the non-local part of the 3-point function.";
\[Lambda]::usage = "\[Lambda] represents the K\[ADoubleDot]llen lambda function.";


(* ::Text:: *)
(*Simplifications.*)


Swap::usage = "Swap[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)] swaps \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\) in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


KToIntegrand::usage = "KToIntegrand[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)] replaces each multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integral in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\) by its integrand with the integration variable \!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\).";


KSimplify::usage = "KSimplify[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] simplifies multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integrals and loop integerals in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


KExpand::usage = "KExpand[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] reduces derivatives of multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integrals to a combination of more elementary functions.";


KExpand::lev = "Level `1` is an invalid expansion level. Level 1 assumed.";


KFullExpand::usage = "KFullExpand[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] expands \!\(\*
StyleBox[\"NL\",\nFontSlant->\"Italic\"]\), its derivatives, and derivatives of multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integrals.";


PolyGammaToHarmonic::usage = "PolyGammaToHarmonic[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] replaces all PolyGamma functions of integral and half-integral arguments in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\) by harmonic numbers.";


TripleK::nonlin = "Warning: `1` parameter in `2` is non-linear in the regulator, \[Epsilon]. Expansion may be invalid.";


(* ::Text:: *)
(*Momentum manipulations.*)


Contract::usage = "Contract[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] contracts all repeated indices of recognized vectors in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontSlant->\"Italic\"]\)
Contract[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \[Mu]] contracts all repeated \[Mu]'s in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).
Contract[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), {\[Mu], ...}] contracts all repreated \[Mu], ... in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).
Contract[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \[Mu], \[Nu]] contracts indices \[Mu] and \[Nu] in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


Contract::argbt = "`1` called with `2` arguments; between `3` and `4` arguments are expected.";


Diff::usage = "Diff[\!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\), \[Mu]] calculates the derivative of \!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\) with respect to \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\)[\[Mu]].";


Diff::darg = "Cannot differentiate over the variable `1`, which appears as the integration variable in `2`.";


(* ::Text:: *)
(*Momentum integrals to multiple-K.*)


LoopToK::usage = "LoopToK[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] replaces all loop integrals by multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integrals.";


LoopToK::barg = "Option `1` `2` should be a Boolean variable.";
LoopToK::loopfail = "Warning: some loop integrals in `1` have not been converted to triple-K integrals. The resulting expression is likely to be incorrect.";


(* ::Text:: *)
(*Divergences.*)


IsDivergent::usage = "IsDivergent[i[\[Alpha],{\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}]] tests if the multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)integral i[\[Alpha],{\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}] is divergent.";


KDivergence::usage = "KDivergence[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] calculates divergent and scheme dependent terms in the \[Epsilon]-expansion of \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


KDivergence::iarg = "Option ExpansionOrder `1` should be an Integer.";
KDivergence::targ = "Option Type `1` should be All, a list of three signs, {\[PlusMinus]1, \[PlusMinus]1, \[PlusMinus]1}, or a list of such lists.";
KDivergence::ksing = "Warning: expression may be indeterminate since it exhibits the `1` singularity.";


(* ::Text:: *)
(*Reduction scheme.*)


KEvaluate::usage = "KEvaluate[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] evaluates all multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integrals in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


LoopEvaluate::usage = "LoopEvaluate[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\)] evaluates all momentum space loop integrals in \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\).";


IsSolvable::usage = "IsSolvable[\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)[\[Alpha],{\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}]] returns True if the multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integral \!\(\*SubscriptBox[
StyleBox[\"i\",\nFontSlant->\"Italic\"], \(\[Alpha], {\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}\)]\) is computable by KEvaluate.
IsSolvable[\!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)[\[Alpha],{\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}][\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(3\)]\)]] returns True if the multiple-\!\(\*
StyleBox[\"K\",\nFontSlant->\"Italic\"]\) integral \!\(\*SubscriptBox[
StyleBox[\"i\",\nFontSlant->\"Italic\"], \(\[Alpha], {\!\(\*SubscriptBox[\(\[Beta]\), \(1\)]\),...}\)]\)[\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(3\)]\)] is computable by KEvaluate.";


(* ::Text:: *)
(*Conformal operators.*)


ScalarKOp::usage = "KOp[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\), \[Beta]] applies single K\!\(\*
StyleBox[\"(\",\nFontSize->12,\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Beta]\",\nFontSize->12,\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontSize->12,\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSize->12]\) operator to \!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\) with respect to momentum magnitude \!\(\*
StyleBox[\"p\",\nFontSlant->\"Italic\"]\) and with parameter \[Beta].";
ScalarKKOp::usage = "KKOp[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"p\", \"i\"],\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"p\", \"j\"],\nFontSlant->\"Italic\"]\), \!\(\*SubscriptBox[\(\[Beta]\), 
StyleBox[\"i\",\nFontSlant->\"Italic\"]]\), \!\(\*SubscriptBox[\(\[Beta]\), 
StyleBox[\"j\",\nFontSlant->\"Italic\"]]\)] = KOp[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"p\", \"i\"],\nFontSlant->\"Italic\"]\), \!\(\*SubscriptBox[\(\[Beta]\), 
StyleBox[\"i\",\nFontSlant->\"Italic\"]]\)] - KOp[\!\(\*
StyleBox[\"expr\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"p\", \"j\"],\nFontSlant->\"Italic\"]\), \!\(\*SubscriptBox[\(\[Beta]\), 
StyleBox[\"j\",\nFontSlant->\"Italic\"]]\)] applies the conformal Ward identity operator in its scalar form.";


(* ::Text:: *)
(*Options.*)


Options[KExpand] = { Level -> 1 };
Options[KSimplify] = { Assumptions :> $Assumptions };
Options[PolyGammaToHarmonic] = { Assumptions :> $Assumptions };
Options[Contract] = { Dimension -> d, Vectors -> Automatic, Indices -> All };
Options[LoopToK] = { Recursive -> True };
Options[IsDivergent] = { Assumptions :> $Assumptions };
Options[KDivergence] = { Type -> All, ExpansionOrder -> 0, uParameter -> u, vParameters -> {v,v,v}, Assumptions :> $Assumptions };
Options[KEvaluate] = { uParameter -> u, vParameters -> {v,v,v}, ExpansionLevel -> 1 };
Options[LoopEvaluate] = { uParameter -> u, vParameters -> {v,v,v}, ExpansionLevel -> 1, Recursive -> True };


Begin["TripleK`Private`"];


(* ::Section:: *)
(*Definitions*)


Attributes[p] = { NHoldAll };
Attributes[\[Delta]] = { Orderless };
Attributes[NL] = { Orderless, NumericFunction };
Attributes[\[Lambda]] = { Orderless, NumericFunction };
Attributes[i] = { NumericFunction };
Attributes[CenterDot] = { Orderless, OneIdentity };


Attributes[NLfun] = { Orderless, NumericFunction };
Attributes[\[Lambda]fun] = { Orderless, NumericFunction };
Attributes[Xfun] = { NumericFunction };
Attributes[Yfun] = { NumericFunction };


NL /: N[NL[a_,b_,c_], opts___] := N[NLfun[a,b,c], opts];
\[Lambda] /: N[\[Lambda][a_,b_,c_], opts___] := N[\[Lambda]fun[a,b,c], opts];
i /: N[i[a_,{b1_,b2_}][p_], opts___] := N[DoubleKValue[a,{b1,b2}][p], opts];
i /: N[i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_], opts___] /; 0<Re[a]+1-Abs[Re[b1]]-Abs[Re[b2]]-Abs[Re[b3]] && 0<Re[p1+p2+p3] := 
		Module[{x}, NIntegrate[x^a p1^b1 p2^b2 p3^b3 BesselK[b1, p1 x] BesselK[b2, p2 x] BesselK[b3, p3 x], {x,0,\[Infinity]}]];


NL[a_?NumberQ,b_?NumberQ,c_?NumberQ] /; AnyTrue[{a,b,c}, InexactNumberQ] := NLfun[a,b,c];
\[Lambda][a_?NumberQ,b_?NumberQ,c_?NumberQ] /; AnyTrue[{a,b,c}, InexactNumberQ] := \[Lambda]fun[a,b,c];
i[a_?NumberQ,{b1_?NumberQ,b2_?NumberQ}][p_?NumberQ] /; AnyTrue[{a,b1,b2,p}, InexactNumberQ] := DoubleKValue[a,{b1,b2}][p];
i[a_?NumberQ,{b1_?NumberQ,b2_?NumberQ,b3_?NumberQ}][p1_?NumberQ,p2_?NumberQ,p3_?NumberQ] /; 
	0<Re[a]+1-Abs[Re[b1]]-Abs[Re[b2]]-Abs[Re[b3]] && 0<Re[p1+p2+p3] && AnyTrue[{a,b1,b2,b3,p1,p2,p3}, InexactNumberQ] := 
	Module[{x}, NIntegrate[x^a p1^b1 p2^b2 p3^b3 BesselK[b1, p1 x] BesselK[b2, p2 x] BesselK[b3, p3 x], {x,0,\[Infinity]}]];


$Assumptions = If[$Assumptions === True, 
	{ p[1]>0, p[2]>0, p[3]>0 }, 
	Union[$Assumptions, { p[1]>0, p[2]>0, p[3]>0 }]];


Format[p[j_][\[Mu]_]] := Subsuperscript[p,j,\[Mu]];
Format[p[j_]] := Subscript[p,j];
Format[\[Delta][i_,j_]] := Subscript[\[Delta],i,j];
Format[i[a_,b_]] := Subscript[i,a,b];
Format[i[a_,b_][p[1],p[2],p[3]]] := Subscript[i,a,b];
Format[i[a_,b_][Subscript[p,1],Subscript[p,2],Subscript[p,3]]] := Subscript[i,a,b];
Format[NL[p[1],p[2],p[3]]] = NL;
Format[\[Lambda][p[1],p[2],p[3]]] = \[Lambda];


Format[LoopIntegral[d_,{\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_,k_]] := With[{numref=num /. { k[\[Mu]_]->Superscript[k,\[Mu]] }, pow1=2\[Delta]1, pow2=2\[Delta]2, pow3=2\[Delta]3}, 
	HoldForm["\!\(\*StyleBox[\"\[Integral]\",FontSize->18]\)"*HoldForm[("\[DifferentialD]"^d*k)/(2*Pi)^d]*HoldForm[numref/(k^(pow3) (k-Subscript[p, 1])^(pow2) (k+Subscript[p, 2])^(pow1))]]];
Format[LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][num_,k_,-p_]] := With[{numref=num /. { k[\[Mu]_]->Superscript[k,\[Mu]] }, pow1=2\[Delta]1, pow2=2\[Delta]2}, 
	HoldForm["\!\(\*StyleBox[\"\[Integral]\",FontSize->18]\)"*HoldForm[("\[DifferentialD]"^d*k)/(2*Pi)^d]*HoldForm[numref/(k^(pow1) (k+p)^(pow2))]]];
Format[LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][num_,k_,p_]] := With[{numref=num /. { k[\[Mu]_]->Superscript[k,\[Mu]] }, pow1=2\[Delta]1, pow2=2\[Delta]2}, 
	HoldForm["\!\(\*StyleBox[\"\[Integral]\",FontSize->18]\)"*HoldForm[("\[DifferentialD]"^d*k)/(2*Pi)^d]*HoldForm[numref/(k^(pow1) (k-p)^(pow2))]]];


(* ::Section:: *)
(*Simplifications*)


Swap[exp_,toswap1stterm_,toswap2ndterm_] := Module[{toswap3rdterm},((exp/.toswap1stterm->toswap3rdterm)/.toswap2ndterm->toswap1stterm)/.toswap3rdterm->toswap2ndterm];
Swap[x___] := Null /; Message[Swap::argrx, "Swap", Length[{x}], 3];


KToIntegrand[exp_, x_] := exp /. { 
	i[a_, {b1_,b2_}][p_] :> x^a p^(b1+b2) BesselK[b1,p x] BesselK[b2,p x],
	i[a_, {b1_,b2_,b3_}][p1_,p2_,p3_] :> 
		x^a p1^b1 p2^b2 p3^b3 BesselK[b1, p1 x] BesselK[b2, p2 x] BesselK[b3, p3 x]
	} /. {
	i[a_, {b1_,b2_,b3_}] :> 
		x^a p[1]^b1 p[2]^b2 p[3]^b3 BesselK[b1, p[1] x] BesselK[b2, p[2] x] BesselK[b3, p[3] x]
	};


KToIntegrand[x___] := Null /; Message[KToIntegrand::argrx, "KToIntegrand", Length[{x}], 2];


suMomentumUnformat = { 			 
			Subsuperscript[p, j_, \[Mu]_] :> p[j][\[Mu]],
			Superscript[Subscript[p, j_], \[Mu]_] :> p[j][\[Mu]],
			Subscript[p, j_] :> p[j],
			Superscript[p[j_], \[Mu]_] :> p[j][\[Mu]] };


Unformat[exp_] := Module[{tmpNL,tmp\[Lambda],ii}, 
	exp /. { Subscript[i, a_, b_] :> i[a,b],
			 Subsuperscript[p, j_, \[Mu]_] :> p[j][\[Mu]],
			 Superscript[Subscript[p, j_], \[Mu]_] :> p[j][\[Mu]],
			 Subscript[p, j_] :> p[j],
			 Superscript[p[j_], \[Mu]_] :> p[j][\[Mu]],
			 Subscript[\[Delta], i_, j_] :> \[Delta][i,j],
			 NL[x__] :> tmpNL[x],
			 Derivative[i__][NL][x__] :> Derivative[i][tmpNL][x],
			 \[Lambda][x__] :> tmp\[Lambda][x],
			 Derivative[i__][\[Lambda]][x__] :> Derivative[i][tmp\[Lambda]][x] }
		/. { NL -> NL[p[1],p[2],p[3]], \[Lambda] -> \[Lambda][p[1],p[2],p[3]] }
		/. { tmpNL -> NL, tmp\[Lambda] -> \[Lambda] }
		/. { i[a_,b_][x___] :> ii[a,b][x],
			 Derivative[k_,m_,n_][i[a_,b_]][x___] :> Derivative[k,m,n][ii[a,b]][x] }
		/. { i[a_,b_] :> i[a,b][p[1],p[2],p[3]] }
		/. { ii -> i }
];


Unformat[x___] := Null /; Message[Unformat::argrx, "Unformat", Length[{x}], 1];


suKDiff = { 
	Derivative[k_,0,0][i[a_,{b1_,b2_,b3_}]][p1_,p2_,p3_] :>
				-D[p1 * i[a+1, {b1-1,b2,b3}][p1,p2,p3],{p1,k-1}],
	Derivative[k_,m_,0][i[a_,{b1_,b2_,b3_}]][p1_,p2_,p3_] :>
				-D[p2 * Derivative[k,0,0][i[a+1, {b1,b2-1,b3}]][p1,p2,p3],{p2,m-1}],
	Derivative[k_,m_,n_][i[a_,{b1_,b2_,b3_}]][p1_,p2_,p3_] :>
				-D[p3 * Derivative[k,m,0][i[a+1, {b1,b2,b3-1}]][p1,p2,p3],{p3,n-1}] };


suKDiffReg = { 
	Derivative[k_,0,0][i[a_,{b1_,b2_,b3_}, u_,v_]][p1_,p2_,p3_] :>
				-D[p1 * i[a+1, {b1-1,b2,b3}, u,v][p1,p2,p3],{p1,k-1}],
	Derivative[k_,m_,0][i[a_,{b1_,b2_,b3_}], u_,v_][p1_,p2_,p3_] :>
				-D[p2 * Derivative[k,0,0][i[a+1, {b1,b2-1,b3}, u,v]][p1,p2,p3],{p2,m-1}],
	Derivative[k_,m_,n_][i[a_,{b1_,b2_,b3_}], u_,v_][p1_,p2_,p3_] :>
				-D[p3 * Derivative[k,m,0][i[a+1, {b1,b2,b3-1}, u,v]][p1,p2,p3],{p3,n-1}] };


suKNegative[assume___] := { 
	i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_] :> 
		p1^(2b1) i[a,{-b1,b2,b3}][p1,p2,p3] /; Simplify[b1<0 /. \[Epsilon]->0,assume, TimeConstraint->0.1],
	i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_] :> 
		p2^(2b2) i[a,{b1,-b2,b3}][p1,p2,p3] /; Simplify[b2<0 /. \[Epsilon]->0,assume, TimeConstraint->0.1],
	i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_] :> 
		p3^(2b3) i[a,{b1,b2,-b3}][p1,p2,p3] /; Simplify[b3<0 /. \[Epsilon]->0,assume, TimeConstraint->0.1] };


suKNegativeReg = { 
	i[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}][p1_,p2_,p3_] :> 
		p1^(2b1) i[a,{-b1,b2,b3}, u,{-v1,v2,v3}][p1,p2,p3] /; b1<0,
	i[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}][p1_,p2_,p3_] :> 
		p2^(2b2) i[a,{b1,-b2,b3}, u,{v1,-v2,v3}][p1,p2,p3] /; b2<0,
	i[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}][p1_,p2_,p3_] :> 
		p3^(2b3) i[a,{b1,b2,-b3}, u,{v1,v2,-v3}][p1,p2,p3] /; b3<0 };


suKOrder[assume___] := { 
	i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_] :> i[a,{b1,b3,b2}][p1,p3,p2] /; Simplify[b2<b3 /. \[Epsilon]->0,assume, TimeConstraint->0.1], 
	i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_] :> i[a,{b2,b1,b3}][p2,p1,p3] /; Simplify[b1<b2 /. \[Epsilon]->0,assume, TimeConstraint->0.1] };


suKInverseOrder = { 
	i[a_,{b1_,b3_,b2_}][p[1],p[3],p[2]] :> i[a,{b1,b2,b3}][p[1],p[2],p[3]],
	i[a_,{b3_,b2_,b1_}][p[3],p[2],p[1]] :> i[a,{b1,b2,b3}][p[1],p[2],p[3]],
	i[a_,{b2_,b1_,b3_}][p[2],p[1],p[3]] :> i[a,{b1,b2,b3}][p[1],p[2],p[3]],
	i[a_,{b2_,b3_,b1_}][p[2],p[3],p[1]] :> i[a,{b1,b2,b3}][p[1],p[2],p[3]],
	i[a_,{b3_,b1_,b2_}][p[3],p[1],p[2]] :> i[a,{b1,b2,b3}][p[1],p[2],p[3]] };


suKOrderReg = { 
	i[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}][p1_,p2_,p3_] :> i[a,{b1,b3,b2}, u,{v1,v3,v2}][p1,p3,p2] /; b2<b3, 
	i[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}][p1_,p2_,p3_] :> i[a,{b2,b1,b3}, u,{v2,v1,v3}][p2,p1,p3] /; b1<b2 };


suKIdentities[assume___] := {
	c_. i[a_,{b1p2_,b2_,b3_}][p1_,p2_,p3_] + cmtwob1p1_. i[am1_,{b1p1_,b2_,b3_}][p1_,p2_,p3_] :> 
		c p1^2 i[a,{b1p1-1,b2,b3}][p1,p2,p3] /; Simplify[b1p2-b1p1==1 && a-am1==1 && cmtwob1p1==-2b1p1*c, assume, TimeConstraint->0.1],
	c_. i[a_,{b1_,b2p2_,b3_}][p1_,p2_,p3_] + cmtwob2p1_. i[am1_,{b1_,b2p1_,b3_}][p1_,p2_,p3_] :> 
		c p2^2 i[a,{b1,b2p1-1,b3}][p1,p2,p3] /; Simplify[b2p2-b2p1==1 && a-am1==1 && cmtwob2p1==-2b2p1*c, assume, TimeConstraint->0.1],
	c_. i[a_,{b1_,b2_,b3p2_}][p1_,p2_,p3_] + cmtwob3p1_. i[am1_,{b1_,b2_,b3p1_}][p1_,p2_,p3_] :>
		c p3^2 i[a,{b1,b2,b3p1-1}][p1,p2,p3] /; Simplify[b3p2-b3p1==1 && a-am1==1 && cmtwob3p1==-2b3p1*c, assume, TimeConstraint->0.1],
	cp_. i[a_,{b1m1_,b2_,b3_}][p1_,p2_,p3_] + c_. i[am1_,{b1_,b2_,b3_}][p1_,p2_,p3_] :>
		c/(2b1) i[a,{b1+1,b2,b3}][p1,p2,p3] /; Simplify[b1-b1m1==1 && a-am1==1 && 2b1 cp==c*p1^2, assume, TimeConstraint->0.1],
	cp_. i[a_,{b1_,b2m1_,b3_}][p1_,p2_,p3_] + c_. i[am1_,{b1_,b2_,b3_}][p1_,p2_,p3_] :>
		c/(2b2) i[a,{b1,b2+1,b3}][p1,p2,p3] /; Simplify[b2-b2m1==1 && a-am1==1 && 2b2 cp==c*p2^2, assume, TimeConstraint->0.1],
	cp_. i[a_,{b1_,b2_,b3m1_}][p1_,p2_,p3_] + c_. i[am1_,{b1_,b2_,b3_}][p1_,p2_,p3_] :>
		c/(2b3) i[a,{b1,b2,b3+1}][p1,p2,p3] /; Simplify[b3-b3m1==1 && a-am1==1 && 2b3 cp==c*p3^2, assume, TimeConstraint->0.1],
	c_. i[a_,{b1p1_,b2_,b3_}][p1_,p2_,p3_] + cmp1_. i[a_,{b1m1_,b2_,b3_}][p1_,p2_,p3_] :>
		2 c (b1p1-1) i[a-1,{b1p1-1,b2,b3}][p1,p2,p3] /; Simplify[b1p1-b1m1==2 && cmp1==-c p1^2, assume, TimeConstraint->0.1],
	c_. i[a_,{b1_,b2p1_,b3_}][p1_,p2_,p3_] + cmp2_. i[a_,{b1_,b2m1_,b3_}][p1_,p2_,p3_] :>
		2 c (b2p1-1) i[a-1,{b1,b2p1-1,b3}][p1,p2,p3] /; Simplify[b2p1-b2m1==2 && cmp2==-c p2^2, assume, TimeConstraint->0.1],
	c_. i[a_,{b1_,b2_,b3p1_}][p1_,p2_,p3_] + cmp3_. i[a_,{b1_,b2_,b3m1_}][p1_,p2_,p3_] :>
		2 c (b3p1-1) i[a-1,{b1,b2,b3p1-1}][p1,p2,p3] /; Simplify[b3p1-b3m1==2 && cmp3==-c p3^2, assume, TimeConstraint->0.1],
	c1_. i[ap1_,{b1m1_,b2_,b3_}][p1_,p2_,p3_] + c2_. i[ap1_,{b1_,b2m1_,b3_}][p1_,p2_,p3_] + c3_. i[ap1_,{b1_,b2_,b3m1_}][p1_,p2_,p3_] :>
		c1/p[1]^2 * (ap1-b1-b2-b3) * i[ap1-1, {b1,b2,b3}][p1,p2,p3] 
		/; Simplify[b1-b1m1==1 && b2-b2m1==1 && b3-b3m1==1 && c1 p[2]^2==c2 p[1]^2 && c1 p[3]^2==c3 p[1]^2 && !TrueQ[IsDivergent[i[ap1-1,{b1,b2,b3}], assume]], assume, TimeConstraint->0.1]
};


suMomentum3Identities[assume___] := {
	c3_. p[3][\[Mu]_] + c1_. p[1][\[Mu]_] :> -c1 p[2][\[Mu]]/; Simplify[c1==c3,assume,TimeConstraint->0.1],
	c3_. p[3][\[Mu]_] + c2_. p[2][\[Mu]_] :> -c2 p[1][\[Mu]]/; Simplify[c2==c3,assume,TimeConstraint->0.1],
	c2_. p[2][\[Mu]_] + c1_. p[1][\[Mu]_] :> -c1 p[3][\[Mu]]/; Simplify[c1==c2,assume,TimeConstraint->0.1]
	 };


SimplifyFunctionKIdentities[assume_] := # /. Union[
	suKIdentities[assume], 
	suKNegative[assume], 
	suKOrder[assume], 
	suMomentum3Identities[assume]]&;


suMomentum = {
	(x_+y_)[j_] :> x[j]+y[j],
	(a_*x_)[j_] :> a*x[j] /; NumericQ[a],
	(-x_)[j_] :> -x[j],
	0[j_] :> 0,
	x_\[CenterDot](y_+z_) :> x\[CenterDot]y + x\[CenterDot]z,
	x_\[CenterDot](a_*y_) :> a(x\[CenterDot]y) /; NumericQ[a],
	x_\[CenterDot](-y_) :> -(x\[CenterDot]y),
	x_\[CenterDot]0 :> 0 };


suMomentumEx[vector_] := { 
	(x_+y_)[j_] :> x[j]+y[j],
	(a_*x_)[j_] :> a*x[j] /; FreeQ[a,p] && FreeQ[a,vector],
	(-x_)[j_] :> -x[j],
	0[j_] :> 0,
	x_\[CenterDot](y_+z_) :> x\[CenterDot]y + x\[CenterDot]z,
	x_\[CenterDot](a_*y_) :> a(x\[CenterDot]y) /; FreeQ[a,p] && FreeQ[a,vector],
	x_\[CenterDot](-y_) :> -(x\[CenterDot]y),
	x_\[CenterDot]0 :> 0 };


suMomentumEx[vectors_List] := { 
	(x_+y_)[j_] :> x[j]+y[j],
	(a_*x_)[j_] :> a*x[j] /; FreeQ[a,p] && AllTrue[vectors, FreeQ[a,#]&],
	(-x_)[j_] :> -x[j],
	0[j_] :> 0,
	x_\[CenterDot](y_+z_) :> x\[CenterDot]y + x\[CenterDot]z,
	x_\[CenterDot](a_*y_) :> a(x\[CenterDot]y) /; FreeQ[a,p] && AllTrue[vectors, FreeQ[a,#]&],
	x_\[CenterDot](-y_) :> -(x\[CenterDot]y),
	x_\[CenterDot]0 :> 0 };


suProductsToMagnitudes = { 
	p[2]\[CenterDot]p[3] -> 1/2 (p[1]^2 - p[2]^2 - p[3]^2),
	p[1]\[CenterDot]p[3] -> 1/2 (p[2]^2 - p[1]^2 - p[3]^2),
	p[1]\[CenterDot]p[2] -> 1/2 (p[3]^2 - p[2]^2 - p[1]^2),
	Sqrt[p[n_]\[CenterDot]p[n_]] :> p[n],
	Sqrt[p\[CenterDot]p] -> p,
	Sqrt[p[n_]^2] :> p[n],
	Sqrt[p^2] -> p,
	p[n_]\[CenterDot]p[n_] :> p[n]^2,
	p\[CenterDot]p -> p^2 };


suMomentum3 = {
	p[3][\[Mu]_] :> -p[1][\[Mu]]-p[2][\[Mu]],
	p[3]\[CenterDot]x_ :> -(p[1]\[CenterDot]x)-(p[2]\[CenterDot]x) };


suLoopLinearity = {
	LoopIntegral[x__][X_ + Y_, y__] :> LoopIntegral[x][X,y] + LoopIntegral[x][Y,y],
	LoopIntegral[x__][-X_, y__] :> -LoopIntegral[x][X,y],
	LoopIntegral[x__][0, y__] :> 0,
	LoopIntegral[x__][a_*X_, k_, y___] :> a*LoopIntegral[x][X,k,y] /; FreeQ[a,k] && a =!= 1,
	LoopIntegral[x__][a_, k_, y___] :> a*LoopIntegral[x][1,k,y] /; FreeQ[a,k] && a =!= 1
};


suLoopReduction = {
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][x_.*(k_\[CenterDot]k_)^n_.,k_,p_] :> 
		LoopIntegral[d, {\[Delta]1-n//Simplify,\[Delta]2}][x, k, p],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][x_.*(k_\[CenterDot]p_)^n_.,k_,p_] :> 
		(p\[CenterDot]p /. suProductsToMagnitudes)/2*LoopIntegral[d, {\[Delta]1,\[Delta]2}][x*(k\[CenterDot]p)^(n-1), k,p]
		-1/2*LoopIntegral[d, {\[Delta]1,\[Delta]2-1}][x*(k\[CenterDot]p)^(n-1), k,p]
		+1/2*LoopIntegral[d, {\[Delta]1-1,\[Delta]2}][x*(k\[CenterDot]p)^(n-1), k,p],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][x_.*(k_\[CenterDot]k_)^n_.,k_] :> 
		LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3-n//Simplify}][x, k],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][x_.*(k_\[CenterDot]p[1])^n_.,k_] :> 
		p[1]^2/2*LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][x*(k\[CenterDot]p[1])^(n-1), k]
		-1/2*LoopIntegral[d, {\[Delta]1,\[Delta]2-1,\[Delta]3}][x*(k\[CenterDot]p[1])^(n-1), k]
		+1/2*LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3-1}][x*(k\[CenterDot]p[1])^(n-1),k],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][x_.*(k_\[CenterDot]p[2])^n_., k_] :>
		-p[2]^2/2*LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][x*(k\[CenterDot]p[2])^(n-1), k]
		+1/2*LoopIntegral[d,{\[Delta]1-1,\[Delta]2,\[Delta]3}][x*(k\[CenterDot]p[2])^(n-1),k]
		-1/2*LoopIntegral[d,{\[Delta]1,\[Delta]2,\[Delta]3-1}][x*(k\[CenterDot]p[2])^(n-1),k],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][x_.*(k_\[CenterDot]k_+2 k_\[CenterDot]p[2]+p[2]^2)^n_.,k_] :> 
		LoopIntegral[d, {\[Delta]1-n//Simplify,\[Delta]2,\[Delta]3}][x, k],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][x_.*(k_\[CenterDot]k_-2 k_\[CenterDot]p[1]+p[1]^2)^n_.,k_] :> 
		LoopIntegral[d, {\[Delta]1,\[Delta]2-n//Simplify,\[Delta]3}][x, k],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][x_.*(k_\[CenterDot]k_-2 k_\[CenterDot]p_+p_^2)^n_.,k_,p_] :> 
		LoopIntegral[d, {\[Delta]1,\[Delta]2-n//Simplify}][x, k, p],
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][x_.*(k_\[CenterDot]k_+2 k_\[CenterDot]p_+p_^2)^n_.,k_, -p_] :> 
		LoopIntegral[d, {\[Delta]1,\[Delta]2-n//Simplify}][x, k, -p]		
};


suLoopZero = {
	LoopIntegral[d_, {n_Integer,\[Delta]2_}][x_,k_,p_] /; n<=0 && FreeQ[x, k\[CenterDot]_]:> 0,
	LoopIntegral[d_, {\[Delta]1_,n_Integer}][x_,k_,p_] /; n<=0 && FreeQ[x, k\[CenterDot]_]:> 0,
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,n_Integer}][x_,k_] /; n<=0 && FreeQ[x, k\[CenterDot]_]:> 
		LoopIntegral[d, {\[Delta]2,\[Delta]1}][x*(k\[CenterDot]k)^(-n) /. {k -> k+p[1]} //. suMomentumEx[k] /. suProductsToMagnitudes, k, p[3]],
	LoopIntegral[d_, {\[Delta]1_,n_Integer,\[Delta]3_}][x_,k_] /; n<=0 && FreeQ[x, k\[CenterDot]_]:> 
		LoopIntegral[d, {\[Delta]3,\[Delta]1}][x*(k\[CenterDot]k - 2 k\[CenterDot]p[1]+p[1]^2)^(-n), k, -p[2]],
	LoopIntegral[d_, {n_Integer,\[Delta]2_,\[Delta]3_}][x_,k_] /; n<=0 && FreeQ[x, k\[CenterDot]_] :> 
		LoopIntegral[d, {\[Delta]3,\[Delta]2}][x*(k\[CenterDot]k + 2 k\[CenterDot]p[2]+p[2]^2)^(-n), k, p[1]]
};


suLoopZeroRecursive = {
	LoopIntegral[d_, {n_Integer,\[Delta]2_}][x_,k_,p_] /; n<=0 && FreeQ[x, LoopIntegral] && FreeQ[x, k\[CenterDot]_] :> 0,
	LoopIntegral[d_, {\[Delta]1_,n_Integer}][x_,k_,p_] /; n<=0 && FreeQ[x, LoopIntegral] && FreeQ[x, k\[CenterDot]_] :> 0,
	LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,n_Integer}][x_,k_] /; n<=0 && FreeQ[x, LoopIntegral] && FreeQ[x, k\[CenterDot]_] :> 
		LoopIntegral[d, {\[Delta]2,\[Delta]1}][x*(k\[CenterDot]k)^(-n) /. {k -> k+p[1]} //. suMomentumEx[k] /. suProductsToMagnitudes, k, p[3]],
	LoopIntegral[d_, {\[Delta]1_,n_Integer,\[Delta]3_}][x_,k_] /; n<=0 && FreeQ[x, LoopIntegral] && FreeQ[x, k\[CenterDot]_] :> 
		LoopIntegral[d, {\[Delta]3,\[Delta]1}][x*(k\[CenterDot]k - 2 k\[CenterDot]p[1]+p[1]^2)^(-n), k, -p[2]],
	LoopIntegral[d_, {n_Integer,\[Delta]2_,\[Delta]3_}][x_,k_] /; n<=0 && FreeQ[x, LoopIntegral] && FreeQ[x, k\[CenterDot]_] :> 
		LoopIntegral[d, {\[Delta]3,\[Delta]2}][x*(k\[CenterDot]k + 2 k\[CenterDot]p[2]+p[2]^2)^(-n), k, p[1]]
};


suLoopReductionZero = Union[suLoopReduction, suLoopZero, suLoopLinearity];
suLoopReductionZeroRecursive = Union[suLoopReduction, suLoopZeroRecursive, suLoopLinearity];


suNLDiff = { 
	Derivative[i_,j_,k_][NL][a_,b_,c_] :> D[NL[1][a,b,c],{c,k},{b,j},{a,i-1}] /; i>=1,
	Derivative[0, j_,k_][NL][a_,b_,c_] :> D[NL[2][a,b,c],{c,k},{b,j-1}] /; j>=1,
	Derivative[0, 0, k_][NL][a_,b_,c_] :> D[NL[3][a,b,c],{c,k-1}] /; k>=1,
	Derivative[i_,j_,k_][NL[1]][a_,b_,c_] :> D[2/(a Sqrt[\[Lambda][a,b,c]])*(a^2 Log[a^2]+1/2(c^2-a^2-b^2) Log[b^2]+1/2(b^2-a^2-c^2) Log[c^2]),{a,i},{b,j},{c,k}],
	Derivative[i_,j_,k_][NL[2]][a_,b_,c_] :> D[2/(b Sqrt[\[Lambda][a,b,c]])*(b^2 Log[b^2]+1/2(c^2-a^2-b^2) Log[a^2]+1/2(a^2-c^2-b^2) Log[c^2]),{a,i},{b,j},{c,k}],
	Derivative[i_,j_,k_][NL[3]][a_,b_,c_] :> D[2/(c Sqrt[\[Lambda][a,b,c]])*(c^2 Log[c^2]+1/2(b^2-a^2-c^2) Log[a^2]+1/2(a^2-c^2-b^2)Log[b^2]),{a,i},{b,j},{c,k}],
	NL[2][a_,b_,c_]:>2/(b Sqrt[\[Lambda][a,b,c]])*(b^2 Log[b^2]+1/2(c^2-a^2-b^2) Log[a^2]+1/2(a^2-c^2-b^2) Log[c^2]),
	NL[3][a_,b_,c_]:>2/(c Sqrt[\[Lambda][a,b,c]])*(c^2 Log[c^2]+1/2(b^2-a^2-c^2) Log[a^2]+1/2(a^2-c^2-b^2) Log[b^2]),
	NL[1][a_,b_,c_]:>2/(a Sqrt[\[Lambda][a,b,c]])*(a^2 Log[a^2]+1/2(c^2-a^2-b^2) Log[b^2]+1/2(b^2-a^2-c^2) Log[c^2]),
	Derivative[i_,j_,k_][\[Lambda]][a_,b_,c_] :> D[\[Lambda]fun[a,b,c],{a,i},{b,j},{c,k}]
};


suPolyGammaToHarmonic[assume___] := { PolyGamma[0, z_] :> HarmonicNumber[z-1] - EulerGamma 
		/; Simplify[z>=1 && Element[z,Integers], assume, TimeConstraint->0.1],
	PolyGamma[n_Integer, z_] :> (-1)^n *n! * (HarmonicNumber[z-1,n+1] - Zeta[n+1]) 
		/; Simplify[n > 0 && z>=1 && Element[z,Integers], assume, TimeConstraint->0.1],
	PolyGamma[0, z_] :> -EulerGamma - 2Log[2] + 2 HarmonicNumber[2z-1] - HarmonicNumber[z-1/2]
		/; Simplify[z>=1/2 && Divisible[2z+1,2], assume, TimeConstraint->0.1],
	PolyGamma[n_Integer, z_] :> (-1)^n *n! * ( 2^(n+1)HarmonicNumber[2z-1, n+1]
			-HarmonicNumber[z-1/2, n+1] - (2^(n+1)-1) Zeta[n+1] )
		/; Simplify[n>0 && z>=1/2 && Divisible[2z+1,2], assume, TimeConstraint->0.1],
	PolyGamma[n_Integer?EvenQ, z_] :> PolyGamma[n, 1-z] 
		/; Simplify[n>=0 && z<=-1/2 && Divisible[2z+1,2], assume, TimeConstraint->0.1],
	PolyGamma[n_Integer?OddQ, z_] :> -PolyGamma[n, 1-z] + (2 Pi)^(n+1) (-1)^((n-1)/2) * (2^(n+1)-1) * BernoulliB[n+1]/(n+1)
		/; Simplify[n>=0 && z<=-1/2 && Divisible[2z+1,2], assume, TimeConstraint->0.1]
};


fcall:PolyGammaToHarmonic[exp_, opts___?OptionQ] := Module[{ValidOpts, assume},
	ValidOpts = First /@ Options[PolyGammaToHarmonic];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[PolyGammaToHarmonic::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	assume = {Assumptions} /. Flatten[{opts}] /. Options[PolyGammaToHarmonic] // Flatten;
	assume = Assumptions -> If[$Assumptions === True, assume, Union[$Assumptions, assume]];

	exp /. suPolyGammaToHarmonic[assume]
];


fcall:PolyGammaToHarmonic[exp_, k___] := Null /; Message[PolyGammaToHarmonic::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
PolyGammaToHarmonic[] := Null /; Message[PolyGammaToHarmonic::argrx, "PolyGammaToHarmonic", 0, 1];


KExpandOnce[exp_, ks_] := 
	Expand[exp] /. If[!FreeQ[exp, LoopIntegral], {
			LoopIntegral[d_,\[Delta]s_][X_,k_,p_] :> LoopIntegral[d,\[Delta]s][KExpandOnce[X, Union[ks, {k}]], k, p],
			LoopIntegral[d_,\[Delta]s_][X_,k_] :> LoopIntegral[d,\[Delta]s][KExpandOnce[X, Union[ks, {k}]], k] 
		}, {}] /. 
	{ i[a_,b_] :> i[Expand[a], Expand[b]] } //. 
	suMomentumEx[ks] //. 
	suLoopLinearity;


KExpandFunctions[exp_, inlevel_] := Module[{level, res, Pow},
	level = If[MemberQ[{0,1,2,3,4, D,Diff,\[Lambda],NL,Integer}, inlevel],
		inlevel,
		Message[KExpand::lev, inlevel]; 1];
	level = Switch[level, D, 1, Diff, 1, Integer, 2, \[Lambda], 3, NL, 4, _, level];

	res = exp
		 /. suPolyGammaToHarmonic[]
		//. If[1 <= level, suNLDiff, {}]
		 /. If[3 <= level, { \[Lambda]->\[Lambda]fun }, {}]
		 /. If[4 <= level, { NL->NLfun }, {}];
	If[level === 2, 
		res /. { Power[\[Lambda][p__],n_] :> Pow[{p},n] /; !IntegerQ[n] }
			/. { \[Lambda][p__]^(n_.) :> \[Lambda]fun[p]^n /; IntegerQ[n] }
			/. { Pow[p_,n_] :> (\[Lambda]@@p)^n },
		res]
];


fcall:KExpand[exp_, opts___?OptionQ] := Module[{ValidOpts, level},
	ValidOpts = First /@ Options[KExpand];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[KExpand::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	level = Level /. Flatten[{opts}] /. Options[KExpand];

	KExpandFunctions[KExpandOnce[Unformat[exp], {}] //. suKDiff, level]
]; 
KExpand[exp_List, opts___?OptionQ] := KExpand[#,opts]& /@ exp;
KExpand[exp_SeriesData, opts___?OptionQ] :=
	SeriesData[exp[[1]],exp[[2]], KExpand[#,opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];


fcall:KExpand[exp_, k___] := Null /; Message[KExpand::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
KExpand[] := Null /; Message[KExpand::argrx, "KExpand", 0, 1];


KFullExpand[exp_] := KExpand[exp, Level -> NL];
KFullExpand[x___] := Null /; Message[KFullExpand::argrx, "KFullExpand", Length[{x}], 1];
KFullExpand[exp_List] := KFullExpand /@ exp;
KFullExpand[exp_SeriesData] :=
	SeriesData[exp[[1]],exp[[2]], KFullExpand /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];


KSimplifyOnce[exp_, ks_, assume_] := 
	Contract[	
		exp /. If[!FreeQ[exp, LoopIntegral], {
			LoopIntegral[d_,\[Delta]s_][X_,k_,p_] :> LoopIntegral[d,\[Delta]s][KSimplifyOnce[X, Union[ks, {k}], assume], k, p],
			LoopIntegral[d_,\[Delta]s_][X_,k_] :> LoopIntegral[d,\[Delta]s][KSimplifyOnce[X, Union[ks, {k}], assume], k]
		}, {}]
		/. { i[a_,b_] :> i[Simplify[a,assume, TimeConstraint->0.1], Simplify[b,assume, TimeConstraint->0.1]] }
		//.suMomentumEx[ks]
		//.{ a_*LoopIntegral[x__][X_, k_, y___] :> LoopIntegral[x][a*X,k,y] /; FreeQ[a,k] && a =!= 1 }
	] //. suLoopReductionZeroRecursive;


fcall:KSimplify[exp_, opts___?OptionQ] := Module[{ii,LI, ValidOpts, defu,defv,assume},
	ValidOpts = First /@ Options[KSimplify];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[KSimplify::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	assume = {Assumptions} /. Flatten[{opts}] /. Options[KSimplify] // Flatten;
	assume = Assumptions -> If[$Assumptions === True, assume, Union[$Assumptions, assume]];

	Simplify[KSimplifyOnce[Unformat[exp], {}, assume] 
			//. suKDiff 
			//. suKNegative[assume] 
			/.  suPolyGammaToHarmonic[assume], 
		assume, 
		TransformationFunctions -> {Automatic, SimplifyFunctionKIdentities[assume]} ] 
	/.  suKInverseOrder 
	//. suNLDiff
];


fcall:KSimplify[exp_, k___] := Null /; Message[KSimplify::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
KSimplify[] := Null /; Message[KSimplify::argrx, "KSimplify", 0, 1];


KSimplify[exp_List, opts___?OptionQ] := KSimplify[#,opts]& /@ exp;


KSimplify[exp_SeriesData, opts___?OptionQ] :=
	SeriesData[exp[[1]],exp[[2]], KSimplify[#,opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];


PrepareIntegral[a_,{b1_,b2_,b3_}, defu_,defv_, assume___] := Module[{acf0,acf1,bcf0,bcf1},
	If[FreeQ[a,\[Epsilon]]&&FreeQ[b1,\[Epsilon]]&&FreeQ[b2,\[Epsilon]]&&FreeQ[b3,\[Epsilon]],
		If[TrueQ[IsDivergent[i[a,{b1,b2,b3}], assume]],
			i[a,{b1,b2,b3}, defu,defv],
			i[a,{b1,b2,b3}, 0,{0,0,0}]
		],
		Quiet[
			acf0 = a /. \[Epsilon]->0;
			bcf0 = {b1,b2,b3} /. \[Epsilon]->0;
			acf1 = D[a,\[Epsilon]] /. \[Epsilon]->0;
			bcf1 = D[{b1,b2,b3},\[Epsilon]] /. \[Epsilon]->0;
		];
		If[!TrueQ[Simplify[a-acf0-acf1 \[Epsilon] == 0, assume, TimeConstraint->0.1]],
			Message[TripleK::nonlin, "\[Alpha]", i[a,{b1,b2,b3}]]];
		If[!TrueQ[Simplify[{b1,b2,b3}-bcf0-bcf1 \[Epsilon] == 0, assume, TimeConstraint->0.1]],
			Message[TripleK::nonlin, "\[Beta]", i[a,{b1,b2,b3}]]];		
		i[acf0,bcf0, acf1,bcf1]
	]
];


PrepareIntegral[a_,{b1_,b2_}, defu_,defv_, assume___] := Module[{acf0,acf1,bcf0,bcf1},
	If[FreeQ[a,\[Epsilon]]&&FreeQ[b1,\[Epsilon]]&&FreeQ[b2,\[Epsilon]],
		If[TrueQ[IsDivergent[i[a,{b1,b2}], assume]],
			i[a,{b1,b2}, defu, {defv[[1]],defv[[2]]}],
			i[a,{b1,b2}, 0,{0,0}]
		],
		Quiet[
			acf0 = a /. \[Epsilon]->0;
			bcf0 = {b1,b2} /. \[Epsilon]->0;
			acf1 = D[a,\[Epsilon]] /. \[Epsilon]->0;
			bcf1 = D[{b1,b2},\[Epsilon]] /. \[Epsilon]->0;
		];
		If[TrueQ[Simplify[a-acf0-acf1 \[Epsilon] != 0, assume, TimeConstraint->0.1]],
			Message[TripleK::nonlin, "\[Alpha]", i[a,{b1,b2}]]];
		If[TrueQ[Simplify[{b1,b2}-bcf0-bcf1 \[Epsilon] != 0, assume, TimeConstraint->0.1]],
			Message[TripleK::nonlin, "\[Beta]", i[a,{b1,b2}]]];		
		i[acf0,bcf0, acf1,bcf1]
	]
];


PrepareExpression[exp_, defu_,defv_, assume___] :=
	Unformat[exp] //. suKDiff /.
		{ i[a_,{b1_,b2_,b3_}] :> PrepareIntegral[a,{b1,b2,b3}, defu,defv,assume],
		  i[a_,{b1_,b2_}] :> PrepareIntegral[a,{b1,b2}, defu,defv,assume] } //. 
	 	suKNegativeReg //. 
	 	suKOrderReg;


suPostpareExpression = { i[a_,b_,u_,v_] :> i[a + u \[Epsilon],b + v \[Epsilon]] };


(* ::Section:: *)
(*Momentum manipulations*)


suContract[d_] := {
		\[Delta][i_,j_]\[Delta][j_,k_] :> \[Delta][i,k] /; !NumericQ[j],
		\[Delta][i_,j_]\[Delta][i_,j_] :> d /; !NumericQ[i]&&!NumericQ[j],
		\[Delta][i_,i_] :> d /; !NumericQ[i],
		x_[i_]y_[i_] :> x\[CenterDot]y /; !NumericQ[i],
		\[Delta][i_,j_]x_[j_] :> x[i] /; !NumericQ[j],
		x_[i_]^2 :> x\[CenterDot]x /; !NumericQ[i]
			};


suContractVectors[d_, vects_] := {
		\[Delta][i_,j_]\[Delta][j_,k_] :> \[Delta][i,k] /; !NumericQ[j],
		\[Delta][i_,j_]\[Delta][i_,j_] :> d /; !NumericQ[i]&&!NumericQ[j],
		\[Delta][i_,i_] :> d /; !NumericQ[i],
		x_[i_]y_[i_] :> x\[CenterDot]y /; !NumericQ[i] && MemberQ[vects, x] && MemberQ[vects, y],
		\[Delta][i_,j_]x_[j_] :> x[i] /; !NumericQ[j] && MemberQ[vects, x],
		x_[i_]^2 :> x\[CenterDot]x /; !NumericQ[i] && MemberQ[vects, x]
			};


suContractIndices[d_, js_List] := {
		\[Delta][i_,j_]\[Delta][j_,k_] :> \[Delta][i,k] /; MemberQ[js,j],
		\[Delta][i_,j_]^2 :> \[Delta][i,i] /; MemberQ[js,j],
		\[Delta][j_,j_] :> d /; MemberQ[js,j],
		x_[j_]y_[j_] :> x\[CenterDot]y /; MemberQ[js,j],
		\[Delta][i_,j_]x_[j_] :> x[i] /; MemberQ[js,j],
		x_[j_]^2 :> x\[CenterDot]x  /; MemberQ[js,j]
			};


suContractIndices[d_, j_] := {
		\[Delta][i_,j]\[Delta][j,k_] :> \[Delta][i,k],
		\[Delta][i_,j]^2 :> \[Delta][i,i],
		\[Delta][j,j] :> d,
		x_[j]y_[j] :> x\[CenterDot]y,
		\[Delta][i_,j]x_[j] :> x[i],
		x_[j]^2 :> x\[CenterDot]x
			};


suContractVectorsIndices[d_, vects_, js_List] := {
		\[Delta][i_,j_]\[Delta][j_,k_] :> \[Delta][i,k] /; MemberQ[js,j],
		\[Delta][i_,j_]^2 :> \[Delta][i,i] /; MemberQ[js,j],
		\[Delta][j_,j_] :> d /; MemberQ[js,j],
		x_[j_]y_[j_] :> x\[CenterDot]y /; MemberQ[js,j] && MemberQ[vects, x] && MemberQ[vects, y],
		\[Delta][i_,j_]x_[j_] :> x[i] /; MemberQ[js,j] && MemberQ[vects, x],
		x_[j_]^2 :> x\[CenterDot]x  /; MemberQ[js,j] && MemberQ[vects, x]
			};


suContractVectorsIndices[d_, vects_, j_] := {
		\[Delta][i_,j]\[Delta][j,k_] :> \[Delta][i,k],
		\[Delta][i_,j]^2 :> \[Delta][i,i],
		\[Delta][j,j] :> d,
		x_[j]y_[j] :> x\[CenterDot]y /; MemberQ[vects, x] && MemberQ[vects, y],
		\[Delta][i_,j]x_[j] :> x[i] /; MemberQ[vects, x],
		x_[j]^2 :> x\[CenterDot]x /; MemberQ[vects, x]
			};


ContractOnce[exp_, vects_, dim_, idx_] :=
	Expand[exp] /. 
		If[!FreeQ[exp, LoopIntegral], {
			LoopIntegral[d_,\[Delta]s_][num_,k_] :> LoopIntegral[d,\[Delta]s][
				ContractOnce[num, Union[vects,{k}], dim, idx], k],
			LoopIntegral[d_,\[Delta]s_][num_,k_,p_] :> LoopIntegral[d,\[Delta]s][
				ContractOnce[num, Union[vects,{k}], dim, idx], k, p]
		}, {}] //.
		If[idx === All,
			suContractVectors[dim, vects],
			suContractVectorsIndices[dim, vects, idx]
		];


DoContract[exp_, dim_, vects_, idx_] := 
	If[vects === All,
		If[idx === All,
			Expand[Unformat[exp]] //. suContract[dim],
			Expand[Unformat[exp]] //. suContractIndices[dim, idx]
		],
		ContractOnce[
			Unformat[exp], 
			If[ListQ[vects], Union[vects, {p, p[1],p[2],p[3]}], {vects,p, p[1],p[2],p[3]}], 
			dim, 
			idx
		]
	] /. suProductsToMagnitudes;


fcall:Contract[exp_, opts___?OptionQ] := Module[{ValidOpts, dim, vects, idx},
	ValidOpts = First /@ Options[Contract];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[Contract::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{dim, vects, idx} = {Dimension, Vectors, Indices} /. Flatten[{opts}] /. Options[Contract];
	DoContract[exp, dim, If[vects === Automatic, If[idx === All, {}, All], vects], idx]
];


fcall:Contract[exp_, \[Mu]_, opts___?OptionQ] := Module[{ValidOpts, dim, vects},
	ValidOpts = First /@ Options[Contract];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[Contract::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{dim, vects} = {Dimension, Vectors} /. Flatten[{opts}] /. Options[Contract];
	DoContract[exp, dim, If[vects === Automatic, All, vects], \[Mu]]
];


fcall:Contract[exp_, \[Mu]_, \[Nu]_, opts___?OptionQ] := Module[{ValidOpts, dim, vects},
	ValidOpts = First /@ Options[Contract];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[Contract::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{dim, vects} = {Dimension, Vectors} /. Flatten[{opts}] /. Options[Contract];
	DoContract[exp /. \[Nu]->\[Mu], dim, If[vects === Automatic, All, vects], \[Mu]]
];


Contract[exp_List, opts___?OptionQ] := Contract[#,opts]& /@ exp;
Contract[exp_List, \[Mu]_, opts___?OptionQ] := Contract[#,\[Mu],opts]& /@ exp;
Contract[exp_List, \[Mu]_, \[Nu]_, opts___?OptionQ] := Contract[#,\[Mu],\[Nu],opts]& /@ exp;


Contract[exp_SeriesData, opts___?OptionQ] :=
	SeriesData[exp[[1]],exp[[2]], Contract[#,opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];
Contract[exp_SeriesData, \[Mu]_, opts___?OptionQ] :=
	SeriesData[exp[[1]],exp[[2]], Contract[#,\[Mu],opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];
Contract[exp_SeriesData, \[Mu]_, \[Nu]_, opts___?OptionQ] :=
	SeriesData[exp[[1]],exp[[2]], Contract[#,\[Mu],\[Nu],opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];


fcall:Contract[exp_, k___] := Null /; Message[Contract::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
Contract[] := Null /; Message[Contract::argbt, "Contract", 0, 1, 3];


OldDiff[exp_, p[mom_], idx_] /; mom===1 || mom===2 :=
Module[{Vtemp, Ltemp, L3temp},
	D[Unformat[exp]
		//.suMomentum3
		//.suMomentum
		/. { p[mom][i_] :> Vtemp[i] }
		/. { p[mom] -> Ltemp, p[3] -> L3temp }
		//.{ Ltemp\[CenterDot]x_ :> p[mom]\[CenterDot]x }
		/. { Ltemp -> Ltemp[p[mom]], 
			 Vtemp[i_] :> Vtemp[p[mom],i],
			 L3temp -> L3temp[p[mom]] },
		p[mom]]
		/. { 
			Derivative[1][Ltemp][p[mom]] -> p[mom][idx] / p[mom],
			Derivative[1][L3temp][p[mom]] -> (p[1][idx] + p[2][idx]) / p[3],
			Derivative[1,0][Vtemp][p[mom],i_] :> \[Delta][idx,i],
			Derivative[0,1][CenterDot][x_, y_] :> x[idx],
			Derivative[1,0][CenterDot][x_, y_] :> y[idx],
			Ltemp[p[mom]] -> p[mom],
			L3temp[p[mom]] -> p[3],
			Vtemp[p[mom],i_] :> p[mom][i]
			 }
];


OldDiff[exp_, mom_, idx_] := Module[{Vtemp, Ltemp},
	D[exp
		//.suMomentumEx[mom]
		/. { mom[i_] :> Vtemp[mom,i] }
		/. { mom :> Ltemp[mom] }
		//.{ Vtemp[Ltemp[mom],i_] :> Vtemp[mom,i], Ltemp[mom]\[CenterDot]x_ :> mom\[CenterDot]x },
		mom]
		/. { 
			Derivative[1][Ltemp][mom] :> mom[idx] / mom,
			Derivative[1,0][Vtemp][mom,i_] :> \[Delta][idx,i],
			Derivative[0,1][CenterDot][x_, y_] :> x[idx],
			Derivative[1,0][CenterDot][x_, y_] :> y[idx],
			Ltemp[mom] :> mom,
			Vtemp[mom,i_] :> mom[i]
			 }
];


DiffOnce[exp_, p[mom_], idx_, ks_] /; mom===1 || mom===2 := 
Module[{loopint, Vtemp, Ltemp, L3temp, ptemp, itemp},
	loopint = !FreeQ[exp, LoopIntegral];
	D[exp
		/. If[loopint,  { 
				LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_,k_] :> 
					LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][num /. p->ptemp, k][Vtemp[itemp]],
				LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][num_,k_,pp_] :> 
					LoopIntegral[d, {\[Delta]1,\[Delta]2}][num /. p->ptemp, k, pp /. p->ptemp][Vtemp[itemp]]
			}, {}]
		//.suMomentum3
		//.suMomentumEx[ks]
		/. { p[mom][i_] :> Vtemp[i] }
		/. { p[mom] -> Ltemp, p[3] -> L3temp }
		//.{ Ltemp\[CenterDot]x_ :> p[mom]\[CenterDot]x }
		/. { Ltemp -> Ltemp[p[mom]], 
			 Vtemp[i_] :> Vtemp[p[mom],i],
			 L3temp -> L3temp[p[mom]] }, 
		p[mom]]
		/. { 
			Derivative[1][Ltemp][p[mom]] -> p[mom][idx] / p[mom],
			Derivative[1][L3temp][p[mom]] -> (p[1][idx] + p[2][idx]) / p[3],
			Derivative[1,0][Vtemp][p[mom],i_] :> \[Delta][idx,i],
			Derivative[0,1][CenterDot][x_, y_] :> x[idx],
			Derivative[1,0][CenterDot][x_, y_] :> y[idx],
			Ltemp[p[mom]] -> p[mom],
			L3temp[p[mom]] -> p[3],
			Vtemp[p[mom],i_] :> p[mom][i]
			 }
		/. If[loopint, {
 			\[Delta][idx, itemp] -> 1,
				Derivative[1][LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_,k_]][p[mom][itemp]] :>
					  If[mom===1, +2\[Delta]2 LoopIntegral[d, {\[Delta]1,\[Delta]2+1,\[Delta]3}][num*(k[idx]-p[1][idx]) /. ptemp->p, k], 0]
					+ If[mom===2, -2\[Delta]1 LoopIntegral[d, {\[Delta]1+1,\[Delta]2,\[Delta]3}][num*(p[2][idx]+k[idx]) /. ptemp->p, k], 0]
					+ LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][DiffOnce[num /. ptemp->p, p[mom], idx, Union[ks, {k}]], k],
				Derivative[1][LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][num_,k_,pp_]][p[mom][itemp]] :>
					+ 2\[Delta]2 LoopIntegral[d, {\[Delta]1,\[Delta]2+1}][num*(k[idx]-pp[idx]) 
									/. ptemp->p 
									/. suMomentum3
									//.suMomentumEx[Union[ks, {k}]], 
							k, pp /. ptemp->p] * D[pp /. ptemp->p /. p[3]->-p[1]-p[2], p[mom]]
					+ LoopIntegral[d, {\[Delta]1,\[Delta]2}][DiffOnce[num /. ptemp->p, p[mom], idx, Union[ks, {k}]], k, pp /. ptemp->p],
				LoopIntegral[d_,{\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_, k_][p[mom][itemp]] :> 
					LoopIntegral[d,{\[Delta]1,\[Delta]2,\[Delta]3}][num /. ptemp->p, k],
				LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][num_, k_, pp_][p[mom][itemp]] :> 
					LoopIntegral[d,{\[Delta]1,\[Delta]2}][num /. ptemp->p, k, pp /. ptemp->p]
			}, {}]
		/. { LoopIntegral[x__][0, y__] :> 0 }
];


DiffOnce[exp_, mom_, idx_, ks_] := Module[{loopint, Vtemp, Ltemp, ptemp, itemp},
	loopint = !FreeQ[exp, LoopIntegral];
	D[exp
		//.suMomentumEx[Flatten[{mom,ks}]]
		/. If[loopint, {
				LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_,k_] :> 
					LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][num /. mom->ptemp, k][Vtemp[mom, itemp]],
				LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][num_,k_,pp_] :> 
					LoopIntegral[d, {\[Delta]1,\[Delta]2}][num /. mom->ptemp, k, pp /. mom->ptemp][Vtemp[mom, itemp]]
			}, {}]
		/. { mom[i_] :> Vtemp[mom,i] }
		/. { mom :> Ltemp[mom] }
		//.{ Vtemp[Ltemp[mom],i_] :> Vtemp[mom,i], Ltemp[mom]\[CenterDot]x_ :> mom\[CenterDot]x },
		mom]
		/. { 
			Derivative[1][Ltemp][mom] :> mom[idx] / mom,
			Derivative[1,0][Vtemp][mom,i_] :> \[Delta][idx,i],
			Derivative[0,1][CenterDot][x_, y_] :> x[idx],
			Derivative[1,0][CenterDot][x_, y_] :> y[idx],
			Ltemp[mom] :> mom,
			Vtemp[mom,i_] :> mom[i]
			 }
		/. If[loopint, {
				\[Delta][idx, itemp] -> 1,
				Derivative[1][LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_,k_]][mom[itemp]] :>
					LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][DiffOnce[num /. ptemp->mom, mom, idx, Union[ks, {k}]], k],
				Derivative[1][LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][num_,k_,pp_]][mom[itemp]] :>
					+ 2\[Delta]2 LoopIntegral[d, {\[Delta]1,\[Delta]2+1}][num*(k[idx]-pp[idx]) 
										/. ptemp->mom 
										//.suMomentumEx[Union[ks, {k}]], 
							k, pp /. ptemp->mom] * D[pp /. ptemp->mom, mom]
					+ LoopIntegral[d, {\[Delta]1,\[Delta]2}][DiffOnce[num /. ptemp->mom, mom, idx, Union[ks, {k}]], k, pp /. ptemp->mom],
				LoopIntegral[d_,{\[Delta]1_,\[Delta]2_,\[Delta]3_}][num_, k_][mom[itemp]] :> 
					LoopIntegral[d,{\[Delta]1,\[Delta]2,\[Delta]3}][num /. ptemp->mom, k],
				LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][num_, k_, pp_][mom[itemp]] :> 
					LoopIntegral[d,{\[Delta]1,\[Delta]2}][num /. ptemp->mom, k, pp /. ptemp->mom]
			}, {}]
		/. { LoopIntegral[x__][0, y__] :> 0 }
];


Diff[exp_, mom_, idx_] := Module[{pos},
	If[FreeQ[exp, LoopIntegral],
		OldDiff[exp, mom /. suMomentumUnformat, idx],
		pos = FirstPosition[{exp}, LoopIntegral[__][_,mom,___]];
		If[Head[pos]===Missing,
			DiffOnce[Unformat[exp], mom /. suMomentumUnformat, idx, {}],
			Message[Diff::darg, mom, Extract[{exp}, pos]]; Null
		]
	] 
	/. suProductsToMagnitudes
];


Diff[exp_List, mom_, idx_] := Diff[#, mom, idx]& /@ exp;
Diff[exp_SeriesData, mom_, idx_] := 
	SeriesData[exp[[1]],exp[[2]], Diff[#, mom, idx]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]];


Diff[x___] := Null /; Message[Diff::argrx, "Diff", Length[{x}], 3];


(* ::Section:: *)
(*Momentum integrals to multiple-K*)


MomentumIntConstant[d_, m_, a_] := 1/( (4 Pi)^(d/2) 2^m ) * a^(-d/2-m);
SchwingerIntegralValue[at_, {\[Delta]1_,\[Delta]2_}][p_] := 
	2^(\[Delta]1+\[Delta]2+2at+3)/( Gamma[-3at - 2\[Delta]1 - 2\[Delta]2] ) * 
		i[-2at - \[Delta]1 - \[Delta]2 - 1, {-2at - \[Delta]1 - 2\[Delta]2, -2at - 2\[Delta]1 - \[Delta]2}][
			Simplify[Sqrt[p\[CenterDot]p] //. suMomentum /. suProductsToMagnitudes, TimeConstraint->0.1]];
SchwingerIntegralValue[at_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}] := 
	2^(at+4)/( Gamma[-2at - \[Delta]1 - \[Delta]2 - \[Delta]3] ) * 
		i[-at - 1, {-at - \[Delta]2 - \[Delta]3, -at - \[Delta]1 - \[Delta]3, -at - \[Delta]1 - \[Delta]2}][
			p[1],p[2],p[3]];


LoopToSchwinger[exp_] := Module[{l, s,st, b1,b2,b3,bt, counter=1,idx,j},
	Expand[Unformat[Expand[exp]]
		/. suProductsToMagnitudes
		/. suMomentum3
		/. { LoopIntegral[x__][X_,y__] :> LoopIntegral[x][Expand[X],y] }
		//.suMomentum
		//.suLoopReductionZero
		/. { 
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][X_,k_,p_] :>
				LoopIntegral[d, {\[Delta]1,\[Delta]2}][
				Expand[X /. { k -> l + s[2]p/st } //. suMomentumEx[l]], 
				l, p],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][X_,k_] :>
				LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][
				Expand[X /. { k -> l + (s[2]*p[1])/st - (s[1]*p[2])/st } 
						//. suMomentumEx[l]], 
				l]
			}
		//. suLoopLinearity
		//. { (l\[CenterDot]x_)^n_. :> l[idx[counter]]x[idx[counter++]] * (l\[CenterDot]x)^(n-1) /; TrueQ[n>0] }
		//. suLoopLinearity
		//.{ l[\[Mu]_]^n_ :> l@@ConstantArray[\[Mu],n], l[\[Mu]___]l[\[Nu]___] :> l[\[Mu],\[Nu]] }
		/. { 
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][l[x___],l,p_] :> 0 /; OddQ[Length[{x}]],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][l[x___],l,p_] :> 
				MomentumIntConstant[d, Length[{x}]/2, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2]) *
					l[x] s[1]^b1 * s[2]^b2 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2}][p],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][1,l,p_] :> 
				MomentumIntConstant[d, 0, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2]) * 
					s[1]^b1 * s[2]^b2 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2}][p],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][l[x___],l] :> 0 /; OddQ[Length[{x}]],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][l[x___],l] :> 
				MomentumIntConstant[d, Length[{x}]/2, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2] Gamma[\[Delta]3]) *
				l[x] s[1]^b1 * s[2]^b2 * s[3]^b3 * st^bt *
				SchwingerIntegral[0, {\[Delta]1, \[Delta]2, \[Delta]3}],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][1,l] :> 
				MomentumIntConstant[d, 0, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2] Gamma[\[Delta]3]) * 
					s[1]^b1 * s[2]^b2 * s[3]^b3 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2, \[Delta]3}]
			}
		//.{ l[x__] :> Sum[\[Delta][First[{x}],{x}[[j]]]l@@Delete[{x},{{1},{j}}], {j,2,Length[{x}]}] }
		/. { l[] -> 1 }
	]
		/. { s[1]^(a1_) s[2]^(a2_)s[3]^(a3_) st^(at_) SchwingerIntegral[At_, {A1_,A2_,A3_}] :>
				SchwingerIntegral[At+at, {A1+a1, A2+a2, A3+a3}],
			s[1]^(a1_) s[2]^(a2_) st^(at_) SchwingerIntegral[At_, {A1_,A2_}][p_] :>
				SchwingerIntegral[At+at, {A1+a1, A2+a2}][p] }
		/. { b1->0, b2->0, b3->0, bt->0 }
	//. suContractIndices[d, idx/@Range[counter-1]]
];


LoopToSchwinger[exp_, ks_] := Module[{l, s,st, b1,b2,b3,bt, counter=1,idx,j, SI},
	Expand[Expand[exp]
		/. { LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][X_,k_,p_] :> LoopIntegral[d,{\[Delta]1,\[Delta]2}][
					Expand[X] //. suMomentumEx[Union[ks,{k}]]
							  /. suProductsToMagnitudes, k, p],
			 LoopIntegral[d_,{\[Delta]1_,\[Delta]2_,\[Delta]3_}][X_,k_] :>  
				LoopIntegral[d,{\[Delta]1,\[Delta]2,\[Delta]3}][X /.  suMomentum3 
											//. suMomentumEx[Union[ks,{k}]] 
											/. suProductsToMagnitudes, k] }
		//. suMomentumEx[ks]
		//. suLoopReductionZeroRecursive
		/. { 
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][X_,k_,p_] :>
				SI[LoopIntegral[d, {\[Delta]1,\[Delta]2}][
				Expand[X /. { k -> l + s[2]p/st } //. suMomentumEx[Union[ks, {l}]]], 
				l, p]],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][X_,k_] :>
				SI[LoopIntegral[d, {\[Delta]1,\[Delta]2,\[Delta]3}][
				Expand[X /. { k -> l + (s[2]*p[1])/st - (s[1]*p[2])/st } 
						//. suMomentumEx[Union[ks, {l}]]], 
				l]]
			}
		//. suLoopLinearity
		//. { (l\[CenterDot]x_)^n_. :> l[idx[counter]]x[idx[counter++]] * (l\[CenterDot]x)^(n-1) /; TrueQ[n>0] }
		//. suLoopLinearity
		//.{ l[\[Mu]_]^n_ :> l@@ConstantArray[\[Mu],n], l[\[Mu]___]l[\[Nu]___] :> l[\[Mu],\[Nu]] }
		/. { 
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][l[x___],l,p_] :> 0 /; OddQ[Length[{x}]],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][l[x___],l,p_] :> 
				MomentumIntConstant[d, Length[{x}]/2, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2]) *
					l[x] s[1]^b1 * s[2]^b2 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2}][p],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_}][1,l,p_] :> 
				MomentumIntConstant[d, 0, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2]) * 
					s[1]^b1 * s[2]^b2 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2}][p],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][l[x___],l] :> 0 /; OddQ[Length[{x}]],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][l[x___],l] :> 
				MomentumIntConstant[d, Length[{x}]/2, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2] Gamma[\[Delta]3]) *
				l[x] s[1]^b1 * s[2]^b2 * s[3]^b3 * st^bt *
				SchwingerIntegral[0, {\[Delta]1, \[Delta]2, \[Delta]3}],
			LoopIntegral[d_, {\[Delta]1_,\[Delta]2_,\[Delta]3_}][1,l] :> 
				MomentumIntConstant[d, 0, st] / (Gamma[\[Delta]1] Gamma[\[Delta]2] Gamma[\[Delta]3]) * 
					s[1]^b1 * s[2]^b2 * s[3]^b3 * st^bt *
					SchwingerIntegral[0, {\[Delta]1, \[Delta]2, \[Delta]3}]
			}
		/. { s[1]^(a1_) s[2]^(a2_)s[3]^(a3_) st^(at_) SchwingerIntegral[At_, {A1_,A2_,A3_}] :>
				SchwingerIntegral[At+at, {A1+a1, A2+a2, A3+a3}],
			s[1]^(a1_) s[2]^(a2_) st^(at_) SchwingerIntegral[At_, {A1_,A2_}][p_] :>
				SchwingerIntegral[At+at, {A1+a1, A2+a2}][p] }
		/. { b1->0, b2->0, b3->0, bt->0 }
		/. { SI[x_] :> x }
		//.{ l[x__] :> Sum[\[Delta][First[{x}],{x}[[j]]]l@@Delete[{x},{{1},{j}}], {j,2,Length[{x}]}] }
		/. { l[] -> 1 }
	]
	//. suContractIndices[d, idx/@Range[counter-1]]
];


SchwingerToK[exp_] := exp /. { SchwingerIntegral -> SchwingerIntegralValue } //. 
		suMomentum /. suProductsToMagnitudes;


SchwingerToK[exp_, ks_] := Simplify[exp /. 
	{ SchwingerIntegral -> SchwingerIntegralValue } /. 
	{ i[a_,{b1_,b2_}] :> DoubleKValue[a,{b1,b2}] }//. 
	suMomentumEx[ks] /. 
	suProductsToMagnitudes];


LoopToKOnce[exp_, ks_] := Module[{res},
	res = LoopToSchwinger[
		exp /. { LoopIntegral[d_,{\[Delta]1_,\[Delta]2_}][X_,k_,p_] :> LoopIntegral[d,{\[Delta]1,\[Delta]2}][
		 			LoopToKOnce[X, Union[ks, {k}]], k,p] /; !FreeQ[X, LoopIntegral],
		 		LoopIntegral[d_,{\[Delta]1_,\[Delta]2_,\[Delta]3_}][X_,k_] :> LoopIntegral[d,{\[Delta]1,\[Delta]2,\[Delta]3}][
		 			LoopToKOnce[X, Union[ks, {k}]], k] /; !FreeQ[X, LoopIntegral] },
		ks];
	If[!FreeQ[res, LoopIntegral],
		Message[LoopToK::loopfail, exp]];
	SchwingerToK[res, ks]
];


fcall:LoopToK[exp_, opts___?OptionQ] := Module[{ValidOpts, recursive},
	ValidOpts = First /@ Options[LoopToK];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[LoopToK::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	recursive = Recursive /. Flatten[{opts}] /. Options[LoopToK]; 
	If [!BooleanQ[recursive], Message[LoopToK::barg, "ExpandDoubleK", recursive]; recursive = ExpandDoubleK /. Options[LoopToK]];

	If[recursive, 
		LoopToKOnce[Unformat[exp] /. suProductsToMagnitudes /.  suMomentum3, {}],
		SchwingerToK[LoopToSchwinger[exp]]
	]
];


fcall:LoopToK[exp_, k___] := Null /; Message[LoopToK::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
LoopToK[] := Null /; Message[LoopToK::argrx, "LoopToK", 0, 1];


(* ::Section:: *)
(*Divergences*)


IsDivergent[i[a_?NumericQ, {b1_?NumericQ, b2_?NumericQ}]] := (
	(IntegerQ[(1+a-b1-b2)/2] && 1+a-b1-b2 <= 0) || 
	(IntegerQ[(1+a+b1-b2)/2] && 1+a+b1-b2 <= 0) ||
	(IntegerQ[(1+a-b1+b2)/2] && 1+a-b1+b2 <= 0) ||
	(IntegerQ[(1+a+b1+b2)/2] && 1+a+b1+b2 <= 0) );
IsDivergent[i[a_, {b1_, b2_}], assume___] := Simplify[
	(Element[(1+a-b1-b2)/2, Integers] && 1+a-b1-b2 <= 0) || 
	(Element[(1+a+b1-b2)/2, Integers] && 1+a+b1-b2 <= 0) ||
	(Element[(1+a-b1+b2)/2, Integers] && 1+a-b1+b2 <= 0) ||
	(Element[(1+a+b1+b2)/2, Integers] && 1+a+b1+b2 <= 0), assume, TimeConstraint->0.1 ];
IsDivergent[i[a_, {b1_, b2_}][p_], assume___] := IsDivergent[i[a, {b1,b2}], assume];
IsDivergent[i[a_, {b1_, b2_}, u_,v_], assume___] := IsDivergent[i[a, {b1,b2}], assume];
IsDivergent[i[a_, {b1_, b2_}, u_,v_][p_], assume___] := IsDivergent[i[a, {b1,b2}], assume];


IsDivergent[i[a_?NumericQ, {b1_?NumericQ, b2_?NumericQ, b3_?NumericQ}]] := (
	(IntegerQ[(1+a-b1-b2-b3)/2] && 1+a-b1-b2-b3 <= 0) || 
	(IntegerQ[(1+a+b1-b2-b3)/2] && 1+a+b1-b2-b3 <= 0) ||
	(IntegerQ[(1+a-b1+b2-b3)/2] && 1+a-b1+b2-b3 <= 0) ||
	(IntegerQ[(1+a+b1+b2-b3)/2] && 1+a+b1+b2-b3 <= 0) ||
	(IntegerQ[(1+a-b1-b2+b3)/2] && 1+a-b1-b2+b3 <= 0) || 
	(IntegerQ[(1+a+b1-b2+b3)/2] && 1+a+b1-b2+b3 <= 0) ||
	(IntegerQ[(1+a-b1+b2+b3)/2] && 1+a-b1+b2+b3 <= 0) ||
	(IntegerQ[(1+a+b1+b2+b3)/2] && 1+a+b1+b2+b3 <= 0) );
IsDivergent[i[a_, {b1_, b2_, b3_}], assume___] := Simplify[ 
	(Element[(1+a-b1-b2-b3)/2, Integers] && 1+a-b1-b2-b3 <= 0) || 
	(Element[(1+a+b1-b2-b3)/2, Integers] && 1+a+b1-b2-b3 <= 0) ||
	(Element[(1+a-b1+b2-b3)/2, Integers] && 1+a-b1+b2-b3 <= 0) ||
	(Element[(1+a+b1+b2-b3)/2, Integers] && 1+a+b1+b2-b3 <= 0) ||
	(Element[(1+a-b1-b2+b3)/2, Integers] && 1+a-b1-b2+b3 <= 0) || 
	(Element[(1+a+b1-b2+b3)/2, Integers] && 1+a+b1-b2+b3 <= 0) ||
	(Element[(1+a-b1+b2+b3)/2, Integers] && 1+a-b1+b2+b3 <= 0) ||
	(Element[(1+a+b1+b2+b3)/2, Integers] && 1+a+b1+b2+b3 <= 0), assume, TimeConstraint->0.1 ];
IsDivergent[i[a_, {b1_, b2_, b3_}][p1_,p2_,p3_], assume___] := IsDivergent[i[a, {b1,b2,b3}], assume];
IsDivergent[i[a_, {b1_, b2_, b3_}, u_, v_], assume___] := IsDivergent[i[a, {b1,b2,b3}], assume];
IsDivergent[i[a_, {b1_, b2_, b3_}, u_, v_][p1_,p2_,p3_], assume___] := IsDivergent[i[a, {b1,b2,b3}], assume];


BesselKCoeff[\[Sigma]_,\[Nu]_][j_]:=(-1)^(j) Gamma[-\[Sigma] \[Nu]-j]/(2^(\[Sigma] \[Nu]+2j+1) j!);
BesselKCoeff[1,n_,v_][j_] := (-1/2)^(1 + 2*j + n)/(v*\[Epsilon]*j!*(j + n)!) + 
 ((-1)^(2*j + n)*2^(-1 - 2*j - n)*
   (Log[2] + PolyGamma[0, 1 + j + n]))/(j!*(j + n)!) - 
 ((-1)^(2*j + n)*2^(-2 - 2*j - n)*v*\[Epsilon]*(Pi^2 + 3*Log[2]^2 + 
    6*Log[2]*PolyGamma[0, 1 + j + n] + 
    3*PolyGamma[0, 1 + j + n]^2 - 3*PolyGamma[1, 1 + j + n]))/
  (3*j!*(j + n)!) + ((-1)^(2*j + n)*2^(-2 - 2*j - n)*v^2*\[Epsilon]^2*
   (Pi^2*Log[2] + Log[2]^3 + Pi^2*PolyGamma[0, 1 + j + n] + 
    3*Log[2]^2*PolyGamma[0, 1 + j + n] + 
    3*Log[2]*PolyGamma[0, 1 + j + n]^2 + PolyGamma[0, 1 + j + n]^
     3 - 3*Log[2]*PolyGamma[1, 1 + j + n] - 
    3*PolyGamma[0, 1 + j + n]*PolyGamma[1, 1 + j + n] + 
    PolyGamma[2, 1 + j + n]))/(3*j!*(j + n)!) - 
 ((-1)^(2*j + n)*2^(-4 - 2*j - n)*v^3*\[Epsilon]^3*
   (7*Pi^4 + 30*Pi^2*Log[2]^2 + 15*Log[2]^4 + 
    60*Pi^2*Log[2]*PolyGamma[0, 1 + j + n] + 
    60*Log[2]^3*PolyGamma[0, 1 + j + n] + 
    30*Pi^2*PolyGamma[0, 1 + j + n]^2 + 
    90*Log[2]^2*PolyGamma[0, 1 + j + n]^2 + 
    60*Log[2]*PolyGamma[0, 1 + j + n]^3 + 
    15*PolyGamma[0, 1 + j + n]^4 - 
    30*Pi^2*PolyGamma[1, 1 + j + n] - 
    90*Log[2]^2*PolyGamma[1, 1 + j + n] - 
    180*Log[2]*PolyGamma[0, 1 + j + n]*PolyGamma[1, 1 + j + n] - 
    90*PolyGamma[0, 1 + j + n]^2*PolyGamma[1, 1 + j + n] + 
    45*PolyGamma[1, 1 + j + n]^2 + 
    60*Log[2]*PolyGamma[2, 1 + j + n] + 
    60*PolyGamma[0, 1 + j + n]*PolyGamma[2, 1 + j + n] - 
    15*PolyGamma[3, 1 + j + n]))/(45*j!*(j + n)!);
BesselKCoeff[-1,n_,v_][j_] := ((-1)^(2*j + n)*2^(-1 - 2*j + n))/(v*\[Epsilon]*j!*(j - n)!) + 
 ((-1)^(2*j + n)*2^(-1 - 2*j + n)*
   (Log[2] + PolyGamma[0, 1 + j - n]))/(j!*(j - n)!) + 
 ((-1)^(2*j + n)*2^(-2 - 2*j + n)*v*\[Epsilon]*(Pi^2 + 3*Log[2]^2 + 
    6*Log[2]*PolyGamma[0, 1 + j - n] + 
    3*PolyGamma[0, 1 + j - n]^2 - 3*PolyGamma[1, 1 + j - n]))/
  (3*j!*(j - n)!) + ((-1)^(2*j + n)*2^(-2 - 2*j + n)*v^2*\[Epsilon]^2*
   (Pi^2*Log[2] + Log[2]^3 + Pi^2*PolyGamma[0, 1 + j - n] + 
    3*Log[2]^2*PolyGamma[0, 1 + j - n] + 
    3*Log[2]*PolyGamma[0, 1 + j - n]^2 + PolyGamma[0, 1 + j - n]^
     3 - 3*Log[2]*PolyGamma[1, 1 + j - n] - 
    3*PolyGamma[0, 1 + j - n]*PolyGamma[1, 1 + j - n] + 
    PolyGamma[2, 1 + j - n]))/(3*j!*(j - n)!) + 
 ((-1)^(2*j + n)*2^(-4 - 2*j + n)*v^3*\[Epsilon]^3*
   (7*Pi^4 + 30*Pi^2*Log[2]^2 + 15*Log[2]^4 + 
    60*Pi^2*Log[2]*PolyGamma[0, 1 + j - n] + 
    60*Log[2]^3*PolyGamma[0, 1 + j - n] + 
    30*Pi^2*PolyGamma[0, 1 + j - n]^2 + 
    90*Log[2]^2*PolyGamma[0, 1 + j - n]^2 + 
    60*Log[2]*PolyGamma[0, 1 + j - n]^3 + 
    15*PolyGamma[0, 1 + j - n]^4 - 
    30*Pi^2*PolyGamma[1, 1 + j - n] - 
    90*Log[2]^2*PolyGamma[1, 1 + j - n] - 
    180*Log[2]*PolyGamma[0, 1 + j - n]*PolyGamma[1, 1 + j - n] - 
    90*PolyGamma[0, 1 + j - n]^2*PolyGamma[1, 1 + j - n] + 
    45*PolyGamma[1, 1 + j - n]^2 + 
    60*Log[2]*PolyGamma[2, 1 + j - n] + 
    60*PolyGamma[0, 1 + j - n]*PolyGamma[2, 1 + j - n] - 
    15*PolyGamma[3, 1 + j - n]))/(45*j!*(j - n)!);


KDivergenceValue[a_, {b1_,b2_}, u_, {v1_,v2_}, pp_, order_, assume___] /; order <= 0 := 
	Series[DoubleKValue[a+u \[Epsilon],{b1+v1 \[Epsilon],b2+v2 \[Epsilon]}][pp], {\[Epsilon],0,order}, assume];


KDivergenceValue[type_, a_?NumericQ, {b1_?NumericQ,b2_?NumericQ,b3_?NumericQ}, u_, {v1_,v2_,v3_}, order_] /; order <= 0 := 
	KDivergenceValue[type, a, {b1,b2,b3}, u, {v1,v2,v3}, order] =
			Module[{cutoff, nn, i1,i2,s},
	Sum[
		nn = (-s . {b1,b2,b3} - a - 1)/2;
		If[TrueQ[0 <= nn && IntegerQ[nn]],
			Sum[Series[cutoff^(-\[Epsilon](u+s . {v1,v2,v3}))/(\[Epsilon](u+s . {v1,v2,v3})) 
					   * p[1]^(2i1 + If[s[[1]]<0,0,2(b1+v1 \[Epsilon])]) 
				       * p[2]^(2i2 + If[s[[2]]<0,0,2(b2+v2 \[Epsilon])]) 
				       * p[3]^(2(nn-i1-i2) + If[s[[3]]<0,0,2(b3+v3 \[Epsilon])])
				       * BesselKCoeff[s[[1]],b1+v1 \[Epsilon]][i1] 
				       * BesselKCoeff[s[[2]],b2+v2 \[Epsilon]][i2] 
				       * BesselKCoeff[s[[3]],b3+v3 \[Epsilon]][nn-i1-i2], {\[Epsilon],0,order}],
				 {i1,0,nn}, {i2,0,nn-i1}],
			0],
		{s, type}
	] /. { Log[cutoff]->0 }
];


SingularQ[n_, assume___] := TrueQ[Simplify[Element[n, Integers] && 0<=n, assume, TimeConstraint->0.1]];


KDivergenceValue[type_, a_, {b1_,b2_,b3_}, u_, {v1_,v2_,v3_}, order_, assume___] /; order <= 0 := 
	KDivergenceValue[type, a, {b1,b2,b3}, u, {v1,v2,v3}, order, assume] =
			Module[{cutoff, nn, s, s1,s2,s3, sing0,sing1,sing2,sing3, 
					j,k, j0,j1,k0,k1, ass, res, factor },

	ass = Assumptions -> Union[Assumptions /. assume, 
			{k>=0, k<=k1, j>=0, j<=j1, Element[k,Integers], Element[j,Integers]}];

	sing0 = SingularQ[(b1+b2+b3-a-1)/2, assume];
	sing1 = SingularQ[#, assume]& /@ { (-b1+b2+b3-a-1)/2, (b1-b2+b3-a-1)/2, (b1+b2-b3-a-1)/2 };
	sing2 = SingularQ[#, assume]& /@ { (-b1-b2+b3-a-1)/2, (-b1+b2-b3-a-1)/2, (b1-b2-b3-a-1)/2 };
	sing3 = SingularQ[(-b1-b2-b3-a-1)/2, assume];

	factor[special_, ss_] := cutoff^(-\[Epsilon](u+ss . {v1,v2,v3}))/(\[Epsilon](u+ss . {v1,v2,v3})) 
					   * p[1]^(2Subscript[K, 2] + If[ss[[1]]<0,0,2(b1+v1 \[Epsilon])]) 
				       * p[2]^(2(nn-Subscript[K, 1]-Subscript[K, 2]) + If[ss[[2]]<0,0,2(b2+v2 \[Epsilon])]) 
				       * p[3]^(2Subscript[K, 1] + If[ss[[3]]<0,0,2(b3+v3 \[Epsilon])])
				       * If[special == 1, BesselKCoeff[ss[[1]],b1,v1][Subscript[K, 2]], BesselKCoeff[ss[[1]],b1+v1 \[Epsilon]][Subscript[K, 2]]]
				       * If[special == 2, BesselKCoeff[ss[[2]],b2,v2][nn-Subscript[K, 1]-Subscript[K, 2]], BesselKCoeff[ss[[2]],b2+v2 \[Epsilon]][nn-Subscript[K, 1]-Subscript[K, 2]]]
				       * If[special == 3, BesselKCoeff[ss[[3]],b3,v3][Subscript[K, 1]], BesselKCoeff[ss[[3]],b3+v3 \[Epsilon]][Subscript[K, 1]]];
	
	If[sing3, Message[KDivergence::ksing, "(+++)"]];
	If[sing2[[1]], Message[KDivergence::ksing, "(++-)"]];
	If[sing2[[2]], Message[KDivergence::ksing, "(+-+)"]];
	If[sing2[[3]], Message[KDivergence::ksing, "(-++)"]];

	Sum[
		nn = Simplify[(-s . {b1,b2,b3} - a - 1)/2, assume, TimeConstraint->0.1];
		If[SingularQ[nn], Sum[Series[factor[0, s], {\[Epsilon],0,order}, assume],
				 {Subscript[K, 1],0,nn}, {Subscript[K, 2],0,nn-Subscript[K, 1]}], 0],
		{s, Complement[type, {{-1,-1,-1}}]}
	] 
	+ If[MemberQ[type, {-1,-1,-1}],
		s = {-1,-1,-1};
		nn = Simplify[(-s . {b1,b2,b3} - a - 1)/2, assume, TimeConstraint->0.1];
		If[sing0 && sing1[[3]], HoldForm[Sum][Series[factor[3, s], {\[Epsilon],0,order}, assume],
			{Subscript[K, 1],b3,nn}, {Subscript[K, 2],0,nn-Subscript[K, 1]}], 0]
		+ If[sing0 && sing1[[1]], HoldForm[Sum][Series[factor[1, s], {\[Epsilon],0,order}, assume], 
			{Subscript[K, 1], 0, Simplify[If[sing0 && sing1[[3]], Min[nn, b3-1], nn], assume, TimeConstraint->0.1]}, 
			{Subscript[K, 2],b1,nn-Subscript[K, 1]}], 0]
		+ If[sing0 && sing1[[2]], HoldForm[Sum][Series[factor[2, s], {\[Epsilon],0,order}, assume], 
			{Subscript[K, 1], 0, Simplify[If[sing0 && sing1[[3]], Min[nn-b2, b3-1], nn-b2], assume, TimeConstraint->0.1]}, 
			{Subscript[K, 2], 0, Simplify[If[sing0 && sing1[[1]], Min[nn-Subscript[K, 1]-b2,b1-1], nn-Subscript[K, 1]-b2], assume, TimeConstraint->0.1]}], 0]
		+ If[sing0, HoldForm[Sum][Series[factor[0, s], {\[Epsilon],0,order}, assume], 
			{Subscript[K, 1], 0, Simplify[If[sing0 && sing1[[3]], Min[nn, b3-1], nn], assume, TimeConstraint->0.1]},
			{Subscript[K, 2], Simplify[If[sing0 && sing1[[2]], Max[0, nn-b2+1-Subscript[K, 1]], 0], assume, TimeConstraint->0.1], 
				 Simplify[If[sing0 && sing1[[1]], Min[nn-Subscript[K, 1],b1-1], nn-Subscript[K, 1]], assume, TimeConstraint->0.1]}], 0]
		/. { HoldForm[Sum][x_, {j_,j0_,j1_}, {k_,k0_,k1_}] :> 
				HoldForm[Sum][x /. {j->j+j0, k->k+k0},{ j,0, Simplify[j1-j0/.{k->k+k0},ass,TimeConstraint->0.1]},
								{k,0,Simplify[k1-k0 /. {j->j+j0},ass,TimeConstraint->0.1]}] }
		/. { HoldForm[Sum][x_, {j_,0,j1_}, {k_,0,k1_}] :> 0 /; Simplify[j1<0 || k1<0, ass,TimeConstraint->0.1] }
		//.{ HoldForm[Sum][x_, {j_,0,j1_}, {k_,0,k1_}] :> 
				HoldForm[Sum][x /. {k->0}, {j,0,j1/.{k->0}}] /; Simplify[k1==0, ass,TimeConstraint->0.1],
			 HoldForm[Sum][x_, {j_,0,j1_}, {k_,0,k1_}] :>
			    HoldForm[Sum][x /. {j->0}, {k,0,k1/.{j->0}}] /; Simplify[j1==0,ass,TimeConstraint->0.1],
			 HoldForm[Sum][x_, {j_,0,j1_}] :> (x /. {j->0}) /; Simplify[j1==0,ass,TimeConstraint->0.1] }, 0]
	/. { Log[cutoff]->0 }
];


AllDivs = Flatten[Table[{s1,s2,s3}, {s1,{-1,1}}, {s2,{-1,1}}, {s3,{-1,1}}],2];


KAllDivergenceValue[a_, {b1_,b2_,b3_}, u_, {v1_,v2_,v3_}, order_, assume___] :=
	KDivergenceValue[AllDivs, a,{b1,b2,b3}, u,{v1,v2,v3}, order, assume];


fcall:KDivergence[exp_, opts___?OptionQ] := Module[{ValidOpts,type,order,defu,defv,assume},
	ValidOpts = First /@ Options[KDivergence];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[KDivergence::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{type,order,defu,defv,assume} = {Type, ExpansionOrder, uParameter, vParameters, {Assumptions}} 
		/. Flatten[{opts}] /. Options[KDivergence];
	If [!IntegerQ[order], Message[KDivergence::iarg, order]];
	assume = Assumptions -> If[$Assumptions === True, Flatten[assume], Union[$Assumptions,Flatten[assume]]];
	type = DeleteDuplicates @ Cases[#, {s1_,s2_,s3_} /; 
							(s1===1||s1===-1) &&
							(s2===1||s2===-1) &&
							(s3===1||s3===-1)]& @
		If[type === All, 
			AllDivs,
			If[!ListQ[type],
				Message[KDivergence::targ, type]; {},
				If[Length[type] > 0,
					If[ListQ[First[type]],
						type,
						{type}]]]];
		
	PrepareExpression[exp, defu,defv, assume]
		/. { i[x__][p__] :> 0 /; !IsDivergent[i@@{x}] }
		/. { i[a_?NumericQ, {b1_?NumericQ,b2_?NumericQ,b3_?NumericQ}, u_, {v1_,v2_,v3_}][p1_,p2_,p3_] :> 
				( KDivergenceValue[type, a,{b1,b2,b3},u,{v1,v2,v3},order] /. {p[1]->p1, p[2]->p2, p[3]->p3} )
			}
		/. { i[a_, {b1_,b2_}, u_, {v1_,v2_}][pp_] :> 
				KDivergenceValue[a,{b1,b2},u,{v1,v2},pp,order, assume],
			 i[a_, {b1_,b2_,b3_}, u_, {v1_,v2_,v3_}][p1_,p2_,p3_] :> 
				( KDivergenceValue[type, a,{b1,b2,b3},u,{v1,v2,v3},order, assume] /. {p[1]->p1, p[2]->p2, p[3]->p3} )
			}
		/. suPolyGammaToHarmonic[assume]
		/. suPostpareExpression
];


KDivergence[exp_SeriesData, opts___?OptionQ] := 
	If[exp[[1]] === \[Epsilon],
		KDivergence[Normal[exp], opts],
		SeriesData[exp[[1]],exp[[2]], KDivergence[#, opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]]];


fcall:KDivergence[exp_, k___] := Null /; Message[KDivergence::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
KDivergence[] := Null /; Message[KDivergence::argrx, "KDivergence", 0, 1];


(* ::Section:: *)
(*Reduction scheme*)


\[Lambda]fun[a_,b_,c_]:=a^4+b^4+c^4 - 2a^2 b^2 - 2 a^2 c^2 - 2 b^2 c^2;
Xfun[a_,b_,c_]:=(-a^2+b^2+c^2-Sqrt[\[Lambda]fun[a,b,c]])/(2c^2);
Yfun[a_,b_,c_]:=(-b^2+a^2+c^2-Sqrt[\[Lambda]fun[a,b,c]])/(2c^2);
NLfun[a_,b_,c_] := Pi^2/6-2Log[a/c]Log[b/c]+Log[Xfun[a,b,c]]Log[Yfun[a,b,c]]-PolyLog[2,Xfun[a,b,c]]-PolyLog[2,Yfun[a,b,c]];


\[Lambda]val = p[1]^4+p[2]^4+p[3]^4 - 2 p[1]^2 p[2]^2-2 p[1]^2 p[3]^2-2 p[2]^2 p[3]^2;
Xval = (-p[1]^2+p[2]^2+p[3]^2-Sqrt[\[Lambda]val])/(2p[3]^2);
Yval = (-p[2]^2+p[1]^2+p[3]^2-Sqrt[\[Lambda]val])/(2p[3]^2);
NLval = Pi^2/6 - 2 Log[p[1]/p[3]]Log[p[2]/p[3]]+Log[Xval]Log[Yval]-PolyLog[2,Xval]-PolyLog[2,Yval];


vt[v_] := v[[1]]+v[[2]]+v[[3]];
vtt[v_] := v[[1]]^2+v[[2]]^2+v[[3]]^2;


i1000 = 1/(2 Sqrt[\[Lambda][p[1],p[2],p[3]]])*NL[p[1], p[2], p[3]] + SeriesData[\[Epsilon],0,{0},0,1,1];


i2111fin = 2 p[1]^2 p[2]^2 p[3]^2 / (\[Lambda][p[1],p[2],p[3]])^(3/2) NL[p[1], p[2], p[3]] + 1/(2 \[Lambda][p[1],p[2],p[3]])( p[1]^2(p[2]^2+p[3]^2-p[1]^2) Log[p[1]^2]+p[2]^2(p[1]^2+p[3]^2-p[2]^2) Log[p[2]^2]+p[3]^2(p[1]^2+p[2]^2-p[3]^2) Log[p[3]^2] );
i2111div[u_,v_] := 1/((u - vt[v]) \[Epsilon])+u / (u - vt[v])(Log[2]-EulerGamma);
i2111[u_,v_] := i2111fin + i2111div[u,v] + SeriesData[\[Epsilon],0,{0},0,1,1];


i0111m2[u_,v_] := 1/(2(vt[v] - u))Sum[p[j]^2/(u - vt[v] + 2 v[[j]]), {j,1,3}];
i0111m1[u_,v_] := 1/4 Sum[p[j]^2 Log[p[j]^2]/(u - vt[v] + 2 v[[j]]), {j,1,3}] + (u (1-2EulerGamma+2Log[2])-vt[v]) / (4(vt[v]-u))Sum[p[j]^2/(u-vt[v]+2v[[j]]), {j,1,3}];
i0111sch[u_,v_] := 1/8 Sum[v[[j]]/(u-vt[v]+2v[[j]])p[j]^2Log[p[j]^2]^2,{j,1,3}]+1/8(u(1-2EulerGamma+2Log[2])-vt[v])Sum[p[j]^2Log[p[j]^2]/(u-vt[v]+2v[[j]]),{j,1,3}]+((vt[v]-u(1-EulerGamma+Log[2]))^2+u^2(EulerGamma-Log[2])^2+1/3 Pi^2 vtt[v])/(8(vt[v]-u))Sum[p[j]^2/(u-vt[v]+2v[[j]]),{j,1,3}];
i0111fin = -1/8 Sqrt[\[Lambda][p[1],p[2],p[3]]]*NL[p[1],p[2],p[3]];
i0111scviol = 1/16( (p[3]^2-p[1]^2-p[2]^2)Log[p[1]^2]Log[p[2]^2]+(p[2]^2-p[1]^2-p[3]^2)Log[p[1]^2]Log[p[3]^2]+(p[1]^2-p[2]^2-p[3]^2)Log[p[2]^2]Log[p[3]^2] );


i0111div[u_,v_] := i0111m2[u,v] / \[Epsilon]^2 + i0111m1[u,v] / \[Epsilon] + i0111sch[u,v];
i0111[u_,v_] := i0111div[u,v] + i0111fin + i0111scviol + SeriesData[\[Epsilon],0,{0},0,1,1];


iahalfhalfhalf[a_] := (Pi/2)^(3/2) * Gamma[a-1/2] (p[1] + p[2] + p[3])^(1/2 - a) + SeriesData[\[Epsilon],0,{0},0,1,1];
iahalfhalfhalfreg[a_,u_,v_] := Module[{div, div0},
	div  = KAllDivergenceValue[a, {1/2,1/2,1/2}, u,v, 0];
	div0 = KAllDivergenceValue[a, {1/2,1/2,1/2}, 1,{0,0,0}, 0];
	iahalfhalfhalf[a + \[Epsilon]] + div - div0 + SeriesData[\[Epsilon],0,{0},0,1,1]
 ];


halfBessel[j_,\[Beta]_] := (Abs[\[Beta]]-1/2+j)!/(2^j * j! * (Abs[\[Beta]]-1/2-j)!);


iahalfhalfhalf[a_,{b1_,b2_,b3_}] := Module[{k1,k2,k3},
	(Pi/2)^(3/2) * Sum[Gamma[a-1/2-k1-k2-k3] * (p[1] + p[2] + p[3])^(1/2-a+k1+k2+k3) 
			* p[1]^(b1-k1-1/2) p[2]^(b2-k2-1/2) p[3]^(b3-k3-1/2)
			* halfBessel[k1,b1] halfBessel[k2,b2] halfBessel[k3,b3],
		{k1,0,Abs[b1]-1/2}, {k2,0,Abs[b2]-1/2}, {k3,0,Abs[b3]-1/2}] + SeriesData[\[Epsilon],0,{0},0,1,1]
];
iahalfhalfhalfreg[a_,{b1_,b2_,b3_}, u_,v_] := Module[{div, div0},
	div  = KAllDivergenceValue[a, {b1,b2,b3}, u,v, 0];
	div0 = KAllDivergenceValue[a, {b1,b2,b3}, 1,{0,0,0}, 0];
	iahalfhalfhalf[a + \[Epsilon], {b1,b2,b3}] + div - div0 + SeriesData[\[Epsilon],0,{0},0,1,1]
];


IntExpK[\[Mu]_,\[Alpha]_,\[Nu]_,\[Beta]_] := 1/(\[Alpha]+\[Beta])^(\[Mu]+\[Nu])Gamma[\[Mu]+\[Nu]]Gamma[\[Mu]-\[Nu]]/Gamma[\[Mu]+1/2]Hypergeometric2F1[\[Mu]+\[Nu],\[Nu]+1/2,\[Mu]+1/2,(\[Alpha]-\[Beta])/(\[Alpha]+\[Beta])];


iaBhalfhalf[a_,{B_,b2_,b3_}] := Module[{k2,k3},
	2^(B-1)*Pi^(3/2)*p[1]^(2B) * Sum[p[2]^(b2-k2-1/2) p[3]^(b3-k3-1/2)
		* halfBessel[k2,b2] halfBessel[k3,b3] * IntExpK[a-k2-k3, p[2]+p[3], B, p[1]],
	{k2,0,Abs[b2]-1/2}, {k3,0,Abs[b3]-1/2}] + SeriesData[\[Epsilon],0,{0},0,1,1]
];
iaBhalfhalfreg[a_,{B_,b2_,b3_}, u_,v_] := Module[{div, div0},
	div  = KAllDivergenceValue[a, {B,b2,b3}, u,v, 0];
	div0 = KAllDivergenceValue[a, {B,b2,b3}, 1,{0,0,0}, 0];
	iaBhalfhalf[a + \[Epsilon], {B,b2,b3}] + div - div0 + SeriesData[\[Epsilon],0,{0},0,1,1]
];


TripleKHalf[a_,{1/2,1/2,1/2}, u_,v_] /; OddQ[2a] && 2a<=1 := 
	TripleKHalf[a,{1/2,1/2,1/2}, u,v] = iahalfhalfhalfreg[a, u, v];
TripleKHalf[a_,{1/2,1/2,1/2}, u_,v_] := TripleKHalf[a,{1/2,1/2,1/2}, u,v] = iahalfhalfhalf[a];
TripleKHalf[a_,{b1_,b2_,b3_}, u_,v_] /; OddQ[2b1] && OddQ[2b2] && OddQ[2b3] && IsDivergent[i[a,{b1,b2,b3}]] := 
	TripleKHalf[a,{b1,b2,b3}, u,v] = iahalfhalfhalfreg[a, {b1,b2,b3}, u, v];
TripleKHalf[a_,{b1_,b2_,b3_}, u_,v_] /; OddQ[2b1] && OddQ[2b2] && OddQ[2b3] := 
	TripleKHalf[a,{b1,b2,b3}, u,v] = Series[iahalfhalfhalf[a + \[Epsilon], {b1,b2,b3}],{\[Epsilon],0,0}];


TripleKHalf[a_,{b1_,b2_,b3_}, u_,v_] /; !OddQ[2b1] && OddQ[2b2] && OddQ[2b3] && IsDivergent[i[a,{b1,b2,b3}]] := 
	TripleKHalf[a,{b1,b2,b3}, u,v] = iaBhalfhalfreg[a, {b1,b2,b3}, u, v];
TripleKHalf[a_,{b1_,b2_,b3_}, u_,v_] /; !OddQ[2b1] &&OddQ[2b2] && OddQ[2b3] := 
	TripleKHalf[a,{b1,b2,b3}, u,v] = Series[iaBhalfhalf[a + \[Epsilon], {b1,b2,b3}],{\[Epsilon],0,0}];
TripleKHalf[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}] /; OddQ[2b1] && !OddQ[2b2] && OddQ[2b3] := 
	TripleKHalf[a,{b1,b2,b3}, u,{v1,v2,v3}] = Swap[iaBhalfhalfreg[a, {b2,b1,b3}, u, {v2,v1,v3}], p[2],p[1]];
TripleKHalf[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}] /; OddQ[2b1] && OddQ[2b2] && !OddQ[2b3] := 
	TripleKHalf[a,{b1,b2,b3}, u,{v1,v2,v3}] = Swap[iaBhalfhalfreg[a, {b3,b2,b1}, u, {v3,v2,v1}], p[3],p[1]];


DoubleKValue[a_,{b1_,b2_}][p_] := 2^(a-2)*p^(-1-a+b1+b2) / Gamma[1+a] * 
	Gamma[(1+a-b1-b2)/2]*Gamma[(1+a+b1-b2)/2]*Gamma[(1+a-b1+b2)/2]*Gamma[(1+a+b1+b2)/2];
DoubleKValue[a_,{b1_,b2_}, u_,{v1_,v2_}][p_] := 
	Series[DoubleKValue[a+u \[Epsilon],{b1+v1 \[Epsilon],b2+v2 \[Epsilon]}][p], {\[Epsilon],0,0}];
(* DoubleKValue[a_,{b1_,b2_}, u_,{v1_,v2_}][p_] := 
	DoubleKValue[a+u \[Epsilon],{b1+v1 \[Epsilon],b2+v2 \[Epsilon]}][p]; *)


MinPos[b_] := Ordering[b,1][[1]];
MaxPos[b_] := Ordering[b,-1][[1]];


n0[a_,b_] := ( Abs[b[[1]]]+Abs[b[[2]]]+Abs[b[[3]]]-a-1 )/2;
n1[a_,b_] := ( Abs[b[[1]]]+Abs[b[[2]]]+Abs[b[[3]]]-2Min[Abs[b]]-a-1 )/2;
n2[a_,b_] := (-Abs[b[[1]]]-Abs[b[[2]]]-Abs[b[[3]]]+2Max[Abs[b]]-a-1 )/2;
n3[a_,b_] := (-Abs[b[[1]]]-Abs[b[[2]]]-Abs[b[[3]]]-a-1 )/2;
GetNs[i[a_,b_]] := { n0[a,b], n1[a,b], n2[a,b], n3[a,b] };


IsSolvable[i[a_,{b1_,b2_}]] = True;
IsSolvable[i[a_,{b1_,b2_}][p_]] = True;
IsSolvable[i[a_Integer,{b1_Integer,b2_Integer,b3_Integer}]] := Module[{ns},
	ns = GetNs[i[a,{b1,b2,b3}]];
	AllTrue[ns, IntegerQ] && ns[[-2]]<0
];
IsSolvable[i[a_,{b1_,b2_,b3_}]] :=
	(OddQ[2 b1] && OddQ[2 b2]) || (OddQ[2 b1] && OddQ[2 b3]) || (OddQ[2 b2] && OddQ[2 b3]);
IsSolvable[i[a_,{b1_,b2_,b3_}][p1_,p2_,p3_]] := IsSolvable[i[a,{b1,b2,b3}]];


IOperator[j_,\[Beta]_][F_] := p[j]^(2 \[Beta]) F;
MOperator[j_,\[Beta]_][F_] := 2\[Beta] F-(p[j] D[F, p[j]]);
LOperator[j_][F_] := -1/p[j] D[F, p[j]];
KOperator[j_,\[Beta]_][F_] := D[F, p[j],p[j]] - (2\[Beta] - 1)/p[j] D[F, p[j]];
BOperator[a_,b_][F_] := PrefBOperator[a,b][ 
	p[1]^2 MOperator[2, b[[2]]-1][MOperator[3, b[[3]]-1][F]] + 
	p[2]^2 MOperator[3, b[[3]]-1][MOperator[1, b[[1]]-1][F]] + 
	p[3]^2 MOperator[1, b[[1]]-1][MOperator[2, b[[2]]-1][F]] ];
PrefBOperator[a_,b_][F_] := 1/(a+1-b[[1]]-b[[2]]-b[[3]]) * F;


TripleKValue[1,{0,0,0}, u_,v_] = i1000;
TripleKValue[2,{1,1,1}, u_,v_] := TripleKValue[2,{1,1,1}, u,v] = i2111[u,v];
TripleKValue[0,{1,1,1}, u_,v_] := TripleKValue[0,{1,1,1}, u,v] = i0111[u,v];


TripleKValue[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}] /; b2<b3 := Swap[TripleKValue[a,{b1,b3,b2}, u,{v1,v3,v2}],p[2],p[3]];
TripleKValue[a_,{b1_,b2_,b3_}, u_,{v1_,v2_,v3_}] /; b1<b2 := Swap[TripleKValue[a,{b2,b1,b3}, u,{v2,v1,v3}],p[1],p[2]];


TripleKValue[a_,b_, u_,v_] /; Min[b]<0 := TripleKValue[a,b, u,v] =
	IOperator[MinPos[b], Min[b]+\[Epsilon]*v[[MinPos[b]]]][
		TripleKValue[a, b*(1-2UnitVector[3, MinPos[b]]), u, v*(1-2UnitVector[3, MinPos[b]])]];


TripleKValue[a_,b_, u_,v_] /; ( (n0[a,b]<0 && 1<a && 0<Max[b]) 
		|| (n1[a,b]<0 && 0<=n0[a,b] && 2<a && Min[b]<a-1) 
		|| (n1[a,b]==0 && 0<n0[a,b] && 0<a && Min[b]<a+1) ) := 
TripleKValue[a,b, u,v] =
	MOperator[MaxPos[b], Max[b]-1+\[Epsilon]*v[[MaxPos[b]]]][
		TripleKValue[a-1, b-UnitVector[3, MaxPos[b]], u, v]];


TripleKValue[a_,b_, u_,v_] /; ( 0<=Min[b] && 0<n1[a,b] ) := 
TripleKValue[a,b, u,v] =
	PrefBOperator[a+u*\[Epsilon],b+v*\[Epsilon]][
		p[1]^2 TripleKValue[a+1, {b[[1]]-1, b[[2]], b[[3]]}, u,v] + 
		p[2]^2 TripleKValue[a+1, {b[[1]], b[[2]]-1, b[[3]]}, u,v] + 
		p[3]^2 TripleKValue[a+1, {b[[1]], b[[2]], b[[3]]-1}, u,v] ];


TripleKValue[a_,b_, u_,v_] /; ( 0<a && n1[a,b]==0 && n0[a,b]==0 ) := 
TripleKValue[a,b, u,v] =
	LOperator[MinPos[b]][TripleKValue[a-1, b+UnitVector[3,MinPos[b]],u,v]];


TripleKValue[a_,{b_,b_,b_}, u_,v_] /; ( (a-1==b && 2<a) || (a+1==b && 0<a) ) := 
TripleKValue[a,b, u,v] =
	BOperator[a+u*\[Epsilon], {b,b,b}+\[Epsilon]*v][TripleKValue[a-1, {b-1,b-1,b-1}, u, v]];


TripleKValue[a_,{0,0,0}, u_,v_] /; 1<a := 
TripleKValue[a,{0,0,0}, u,v] =
	KOperator[1, \[Epsilon]*v[[1]]][TripleKValue[a-2, {0,0,0},u,v]];


fcall:KEvaluate[exp_, opts___?OptionQ] := Module[{ValidOpts,defu,defv, level},
	ValidOpts = First /@ Options[KEvaluate];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[KEvaluate::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{defu,defv, level} = {uParameter, vParameters, ExpansionLevel} 
		/. Flatten[{opts}] /. Options[KEvaluate];
		
	PrepareExpression[exp, defu,defv]
		/. { i[a_?NumericQ, {b1_?NumericQ,b2_?NumericQ}, u_, {v1_,v2_}][p_] :> 
				DoubleKValue[a,{b1,b2},u,{v1,v2}][p],
			 i[a_, {b1_,b2_}, u_, {v1_,v2_}][p_] :> 
				DoubleKValue[a + u \[Epsilon],{b1 + v1 \[Epsilon], b2 + v2 \[Epsilon]}][p],
			 i[a_Integer, {b1_Integer,b2_Integer,b3_Integer}, u_, {v1_,v2_,v3_}][p1_,p2_,p3_] :> 
				( TripleKValue[a,{b1,b2,b3},u,{v1,v2,v3}] /. {p[1]->p1, p[2]->p2, p[3]->p3} ) 
				/; TrueQ[IsSolvable[i[a,{b1,b2,b3}]]],
			 i[a_, {b1_,b2_,b3_}, u_, {v1_,v2_,v3_}][p1_,p2_,p3_] :>
				( TripleKHalf[a,{b1,b2,b3},u,{v1,v2,v3}] /. {p[1]->p1, p[2]->p2, p[3]->p3} )
				/; (OddQ[2 b1] && OddQ[2 b2]) || (OddQ[2 b1] && OddQ[2 b3]) || (OddQ[2 b2] && OddQ[2 b3])
			}
		/. suPolyGammaToHarmonic[]
		/. suPostpareExpression
		// KExpandFunctions[#, level]&
];


KEvaluate[exp_SeriesData, opts___?OptionQ] := 
	If[exp[[1]] === \[Epsilon],
		KEvaluate[Normal[exp], opts],
		SeriesData[exp[[1]],exp[[2]], KEvaluate[#, opts]& /@ (exp[[3]]), exp[[4]],exp[[5]],exp[[6]]]];


fcall:KEvaluate[exp_, k___] := Null /; Message[KEvaluate::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
KEvaluate[] := Null /; Message[KEvaluate::argrx, "KEvaluate", 0, 1];


fcall:LoopEvaluate[exp_, opts___?OptionQ] := Module[{ValidOpts, defu,defv,level,recursive},
	ValidOpts = First /@ Options[LoopEvaluate];
	Scan[If[!MemberQ[ValidOpts, First[#]],
		Message[LoopEvaluate::optx, ToString[First[#]], ToString[Unevaluated[fcall]]]]&, Flatten[{opts}]];
	{recursive, defu,defv, level} = {Recursive, uParameter, vParameters, ExpansionLevel} 
		/. Flatten[{opts}] /. Options[LoopEvaluate];	

	KEvaluate[LoopToK[exp, Recursive -> recursive],
		{ uParameter -> defu, vParameters -> defv, ExpansionLevel -> level }]
];


fcall:LoopEvaluate[exp_, k___] := Null /; Message[LoopEvaluate::nonopt, Last[{k}], 1, ToString[Unevaluated[fcall]]];
LoopEvaluate[] := Null /; Message[LoopEvaluate::argrx, "LoopEvaluate", 0, 1];


(* ::Section:: *)
(*Conformal operators*)


ScalarKOp[F_, p_, \[Beta]_] := D[F,p,p] - (2\[Beta]-1)/p * D[F,p];
ScalarKKOp[F_,p1_,p2_,\[Beta]1_,\[Beta]2_] := ScalarKOp[F, p1, \[Beta]1] - ScalarKOp[F, p2, \[Beta]2];


ScalarKOp[x___] := Null /; Message[ScalarKOp::argrx, "KOp", Length[{x}], 3];
ScalarKKOp[x___] := Null /; Message[ScalarKOp::argrx, "KOp", Length[{x}], 5];


(* ::Section:: *)
(*End*)


End[];
Protect @@ Names["TripleK`*"];
EndPackage[];

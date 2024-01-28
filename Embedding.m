(* ::Package:: *)

(* ::Chapter:: *)
(*Front End*)


(* ::Input::Initialization:: *)
SortNumber2[n_Integer]:=Module[{a1,a2},
a1=Max[x/.Solve[1/2 x (1+ x)==n,x,Reals]];
a1=Ceiling[a1];
a2=n-1/2 a1 (a1-1);
{a1,a2}]
SortNumber3[n_Integer]:=Module[{a1,a2,a3},
a1=Max[x/.Solve[1/6 x (2+3 x+x^2)==n,x,Reals]];
a1=Ceiling[a1];
{a2,a3}=SortNumber2[n-1/6 a1 (a1^2-1)];
{a1,a2,a3}]
SortNumber4[n_Integer]:=Module[{a1,a2,a3,a4},
a1=Max[x/.Solve[1/24 x (6+11 x+6 x^2+x^3)==n,x,Reals]];
a1=Ceiling[a1];
{a2,a3,a4}=SortNumber3[n-1/24 *(-1+a1) a1 (2+3 a1+a1^2)];
{a1,a2,a3,a4}]


(* ::Section:: *)
(*p6mm*)


(* ::Input::Initialization:: *)
p6m1dIrrep={{1,1,1,1,0},{-1,1,1,1,c},{1,-1,1,1,r},{-1,-1,1,1,c+r}};
p6m2dIrrep={{({
 {1, 0},
 {0, -1}
}),({
 {1, 0},
 {0, 1}
}),({
 {Cos[2Pi/3], Sin[2Pi/3]},
 {-Sin[2Pi/3], Cos[2Pi /3]}
}),({
 {Cos[2Pi/3], Sin[2Pi/3]},
 {-Sin[2Pi/3], Cos[2Pi /3]}
}),c},{({
 {1, 0},
 {0, -1}
}),({
 {-1, 0},
 {0, -1}
}),({
 {Cos[2Pi/3], Sin[2Pi/3]},
 {-Sin[2Pi/3], Cos[2Pi /3]}
}),({
 {Cos[2Pi/3], Sin[2Pi/3]},
 {-Sin[2Pi/3], Cos[2Pi /3]}
}),c},{({
 {Cos[Pi/3], Sin[Pi/3]},
 {-Sin[Pi/3], Cos[Pi /3]}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),r},{({
 {Cos[2Pi/3], Sin[2Pi/3]},
 {-Sin[2Pi/3], Cos[2Pi /3]}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),r}};
p6m3dIrrep={{({
 {0, 1, 0},
 {0, 0, 1},
 {1, 0, 0}
}),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", "1", "0"},
{"1", "0", "0"},
{"0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {1, 0, 0},
 {0, -1, 0},
 {0, 0, -1}
}),({
 {-1, 0, 0},
 {0, 1, 0},
 {0, 0, -1}
}),r},{({
 {0, 1, 0},
 {0, 0, 1},
 {-1, 0, 0}
}),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "1"}], "0"},
{
RowBox[{"-", "1"}], "0", "0"},
{"0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {1, 0, 0},
 {0, -1, 0},
 {0, 0, -1}
}),({
 {-1, 0, 0},
 {0, 1, 0},
 {0, 0, -1}
}),c+r},{({
 {0, 1, 0},
 {0, 0, 1},
 {1, 0, 0}
}),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "1"}], "0"},
{
RowBox[{"-", "1"}], "0", "0"},
{"0", "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {1, 0, 0},
 {0, -1, 0},
 {0, 0, -1}
}),({
 {-1, 0, 0},
 {0, 1, 0},
 {0, 0, -1}
}),0},{({
 {0, 1, 0},
 {0, 0, 1},
 {-1, 0, 0}
}),\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", "1", "0"},
{"1", "0", "0"},
{"0", "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {1, 0, 0},
 {0, -1, 0},
 {0, 0, -1}
}),({
 {-1, 0, 0},
 {0, 1, 0},
 {0, 0, -1}
}),c}};
p6m4dIrrep={({
 {1/2, 0, 0, Sqrt[3]/2},
 {0, -(1/2), -(Sqrt[3]/2), 0},
 {0, Sqrt[3]/2, -(1/2), 0},
 {-(Sqrt[3]/2), 0, 0, 1/2}
}),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"-", "1"}], "0", "0", "0"},
{"0", 
RowBox[{"-", "1"}], "0", "0"},
{"0", "0", "1", "0"},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{
RowBox[{"-", 
FractionBox["1", "2"]}], 
FractionBox[
SqrtBox["3"], "2"], "0", "0"},
{
RowBox[{"-", 
FractionBox[
SqrtBox["3"], "2"]}], 
RowBox[{"-", 
FractionBox["1", "2"]}], "0", "0"},
{"0", "0", 
RowBox[{"-", 
FractionBox["1", "2"]}], 
RowBox[{"-", 
FractionBox[
SqrtBox["3"], "2"]}]},
{"0", "0", 
FractionBox[
SqrtBox["3"], "2"], 
RowBox[{"-", 
FractionBox["1", "2"]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {-(1/2), Sqrt[3]/2, 0, 0},
 {-(Sqrt[3]/2), -(1/2), 0, 0},
 {0, 0, -(1/2), -(Sqrt[3]/2)},
 {0, 0, Sqrt[3]/2, -(1/2)}
})};


(* ::Input::Initialization:: *)
p6mO2MatrixForm[i_Integer]:=Module[{a1,a2,t1,t2,m,tab},
tab=Table[0,{m,1,6}];

If[i<=36,{a1,a2}=SortNumber2[i];
t1=Quotient[a1-1,4];
t2=Quotient[a2-1,4];
a1=Mod[a1-1,4]+1;
a2=Mod[a2-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p6m1dIrrep[[a1]][[m]],p6m1dIrrep[[a2]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2}];
tab[[6]]=PolynomialMod[p6m1dIrrep[[a1]][[5]]+p6m1dIrrep[[a2]][[5]]+t*(t1+t2),2];

];

If[i>36,a1=i-36;
t1=Quotient[a1-1,4];
a1=Mod[a1-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=p6m2dIrrep[[a1]][[m]];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1}];
tab[[6]]=p6m2dIrrep[[a1]][[5]];];

tab
]

p6mO3MatrixForm[i_Integer]:=Module[{a1,a2,a3,t1,t2,t3,m,tab},
tab=Table[0,{m,1,6}];

If[i<=120,{a1,a2,a3}=SortNumber3[i];
t1=Quotient[a1-1,4];
t2=Quotient[a2-1,4];
t3=Quotient[a3-1,4];
a1=Mod[a1-1,4]+1;
a2=Mod[a2-1,4]+1;
a3=Mod[a3-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p6m1dIrrep[[a1]][[m]],p6m1dIrrep[[a2]][[m]],p6m1dIrrep[[a3]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2,(-1)^t3}];
tab[[6]]=PolynomialMod[p6m1dIrrep[[a1]][[5]]+p6m1dIrrep[[a2]][[5]]+p6m1dIrrep[[a3]][[5]]+t*(t1+t2+t3),2];];

If[i>120&&i<=184,a1=i-120;
a2=Quotient[a1-1,8]+1;
a1=Mod[a1-1,8]+1;

t1=Quotient[a1-1,4];
a1=Mod[a1-1,4]+1;
t2=Quotient[a2-1,4];
a2=Mod[a2-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p6m2dIrrep[[a1]][[m]],0},{0,p6m1dIrrep[[a2]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2}];
tab[[6]]=PolynomialMod[p6m2dIrrep[[a1]][[5]]+p6m1dIrrep[[a2]][[5]]+t*t2,2];];

If[i>184,a1=i-184;
t1=Quotient[a1-1,4];
a1=Mod[a1-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=p6m3dIrrep[[a1]][[m]];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t1}];
tab[[6]]=p6m3dIrrep[[a1]][[5]]+t*t1;];

tab
]

p6mO4MatrixForm[i_Integer]:=Module[{a1,a2,a3,a4,t1,t2,t3,t4,m,tab},
tab=Table[0,{m,1,6}];

If[i<=330,{a1,a2,a3,a4}=SortNumber4[i];
t1=Quotient[a1-1,4];
t2=Quotient[a2-1,4];
t3=Quotient[a3-1,4];
t4=Quotient[a4-1,4];
a1=Mod[a1-1,4]+1;
a2=Mod[a2-1,4]+1;
a3=Mod[a3-1,4]+1;
a4=Mod[a4-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p6m1dIrrep[[a1]][[m]],p6m1dIrrep[[a2]][[m]],p6m1dIrrep[[a3]][[m]],p6m1dIrrep[[a4]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2,(-1)^t3,(-1)^t4}];
tab[[6]]=PolynomialMod[p6m1dIrrep[[a1]][[5]]+p6m1dIrrep[[a2]][[5]]+p6m1dIrrep[[a3]][[5]]+p6m1dIrrep[[a4]][[5]]+t*(t1+t2+t3+t4),2];];

If[i>330&&i<=618,a1=i-330;
a2=Quotient[a1-1,8]+1;
a1=Mod[a1-1,8]+1;
{a2,a3}=SortNumber2[a2];

t1=Quotient[a1-1,4];
a1=Mod[a1-1,4]+1;
t2=Quotient[a2-1,4];
a2=Mod[a2-1,4]+1;
t3=Quotient[a3-1,4];
a3=Mod[a3-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p6m2dIrrep[[a1]][[m]],0,0},{0,p6m1dIrrep[[a2]][[m]],0},{0,0,p6m1dIrrep[[a3]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2,(-1)^t3}];
tab[[6]]=PolynomialMod[p6m2dIrrep[[a1]][[5]]+p6m1dIrrep[[a2]][[5]]+p6m1dIrrep[[a3]][[5]]+t*(t2+t3),2];];

If[i>618&&i<=654,a1=i-618;
{a1,a2}=SortNumber2[a1];
t1=Quotient[a1-1,4];
t2=Quotient[a2-1,4];
a1=Mod[a1-1,4]+1;
a2=Mod[a2-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p6m2dIrrep[[a1]][[m]],0},{0,p6m2dIrrep[[a2]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2,(-1)^t2}];
tab[[6]]=PolynomialMod[p6m2dIrrep[[a1]][[5]]+p6m2dIrrep[[a2]][[5]],2];];

If[i>654&&i<=718,a1=i-654;
a2=Quotient[a1-1,8]+1;
a1=Mod[a1-1,8]+1;

t1=Quotient[a1-1,4];
a1=Mod[a1-1,4]+1;
t2=Quotient[a2-1,4];
a2=Mod[a2-1,4]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p6m1dIrrep[[a1]][[m]],0},{0,p6m3dIrrep[[a2]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2,(-1)^t2,(-1)^t2}];
tab[[6]]=PolynomialMod[p6m3dIrrep[[a2]][[5]]+p6m1dIrrep[[a1]][[5]]+t*(t1+t2),2];];

If[i>718&&i<=720,a1=i-718;
t1=Quotient[a1-1,2];

For[m=1,m<=4,m++,
tab[[m]]=p6m4dIrrep[[m]];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t1,(-1)^t1}];
tab[[6]]=0;];

tab
]


(* ::Input::Initialization:: *)
p6mDQCPEmbedding[i_Integer]:=Module[{O2Mat,m,w1,O3Sign,O3T,tab},
O2Mat = p6mO2MatrixForm[i];

w1=PolynomialMod[t+r+O2Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=4,m++,If[p6m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p6m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p6m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O2Mat[[m]]}})];];
tab[[5]]=ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O2Mat[[5]]}})];

tab
]


(* ::Input::Initialization:: *)
p6mDSLEmbedding[i_Integer,j_Integer]:=Module[{O2Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O2Mat =p6mO2MatrixForm[j];
O3Mat= p6mO3MatrixForm[i];

w1=PolynomialMod[t+r+O2Mat[[6]]+O3Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=4,m++,If[p6m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p6m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p6m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O3Mat[[m]]}})],O2Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O3Mat[[5]]}})],O2Mat[[5]]};

tab
]


(* ::Input::Initialization:: *)
p6mSLEmbedding[i_Integer,j_Integer]:=Module[{O4Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O3Mat =p6mO3MatrixForm[j];
O4Mat= p6mO4MatrixForm[i];

w1=PolynomialMod[t+r+O3Mat[[6]]+O4Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=4,m++,If[p6m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p6m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p6m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O4Mat[[m]]}})],O3Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O4Mat[[5]]}})],O3Mat[[5]]};

tab
]


(* ::Input::Initialization:: *)
p6mEmbedding[n_Integer,lsm_Integer,p_Integer]:=Module[{a1,a2,tab},
If[n==5,
a1=data[[1]][[1]][[lsm]][[p]];
tab = p6mDQCPEmbedding[a1];];
If[n==6,
a1=data[[2]][[1]][[lsm]][[p]][[1]];
a2=data[[2]][[1]][[lsm]][[p]][[2]];
tab = p6mDSLEmbedding[a1,a2];];
If[n==7,
a1=data[[3]][[1]][[lsm]][[p]][[1]];
a2=data[[3]][[1]][[lsm]][[p]][[2]];
tab = p6mSLEmbedding[a1,a2];];
tab]


p6mStability[n_Integer,lsm_Integer,p_Integer]:=Module[{check},
check=False;
If[n==5,
If[lsm==1,If[(p==9)||(p==10),check=True;];];
If[lsm==4,If[(p==2),check=True;];];];
If[n==6,
If[lsm==1,If[(p==48),check=True;];];
If[lsm==2,check=True;];
If[lsm==3,check=True;];
If[lsm==4,If[(p==3)||(p==12)||(p==13)||(p==15)||(p==18),check=True;];];];
If[n==7,
If[lsm==2,If[(p==1)||(p==2)||(p==3)||(p==5)||(p==8)||(p==9)||(p==32)||(p==34),check=True;];];
If[lsm==3,If[(p==1)||(p==2)||(p==3)||(p==5)||(p==7)||(p==10)||(p==12)||(p==22)||(p==26),check=True;];];
If[lsm==4,If[(p==99)||(p==164),check=True;];];];
check]


p6mPrintEmbedding[n_Integer,lsm_Integer,p_Integer]:=Module[{flag},
If[n==5,
Print["C6:n->",p6mEmbedding[5,lsm,p][[1]]//MatrixForm,"n"];
Print["M:n->",p6mEmbedding[5,lsm,p][[2]]//MatrixForm,"n"];
Print["T1:n->",p6mEmbedding[5,lsm,p][[3]]//MatrixForm,"n"];
Print["T2:n->",p6mEmbedding[5,lsm,p][[4]]//MatrixForm,"n"];
Print["T:n->",p6mEmbedding[5,lsm,p][[5]]//MatrixForm,"n"];];
If[n==6,
Print["C6:n->",p6mEmbedding[6,lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[6,lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p6mEmbedding[6,lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[6,lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p6mEmbedding[6,lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[6,lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p6mEmbedding[6,lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[6,lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p6mEmbedding[6,lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[6,lsm,p][[5]][[2]]]//MatrixForm];];
If[n==7,
Print["C6:n->",p6mEmbedding[7,lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[7,lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p6mEmbedding[7,lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[7,lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p6mEmbedding[7,lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[7,lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p6mEmbedding[7,lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[7,lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p6mEmbedding[7,lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p6mEmbedding[7,lsm,p][[5]][[2]]]//MatrixForm];];
flag=True;
flag]


p6mSL5dEmbedding[i_Integer,j_Integer]:=Module[{O2Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O3Mat =p6mO3MatrixForm[j];
O2Mat= p6mO2MatrixForm[i];

w1=PolynomialMod[t+r+O3Mat[[6]]+O2Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=4,m++,If[p6m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p6m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p6m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[5],0},{0,O2Mat[[m]]}})],O3Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[5],0},{0,O2Mat[[5]]}})],O3Mat[[5]]};

tab
]


p6m5dEmbedding[lsm_Integer,p_Integer]:=Module[{a1,a2,tab},
tab=0;

a1=dataSL5Rep[[1]][[lsm]][[p]][[1]];
a2=dataSL5Rep[[1]][[lsm]][[p]][[2]];
tab = p6mSL5dEmbedding[a1,a2];

tab
]


p6m5dStability[lsm_Integer,p_Integer]:=Module[{check},
check=False;
If[lsm==1,If[(p==3)||(p==11),check=True;];];
check]


p6m5dPrintEmbedding[lsm_Integer,p_Integer]:=Module[{flag},
Print["C6:n->",p6m5dEmbedding[lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p6m5dEmbedding[lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p6m5dEmbedding[lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p6m5dEmbedding[lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p6m5dEmbedding[lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p6m5dEmbedding[lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p6m5dEmbedding[lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p6m5dEmbedding[lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p6m5dEmbedding[lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p6m5dEmbedding[lsm,p][[5]][[2]]]//MatrixForm];
flag=True;
flag]


(* ::Section:: *)
(*p4mm*)


(* ::Input::Initialization:: *)
p4m1dIrrep={{1,1,1,1,0},{-1,1,1,1,c},{1,-1,1,1,r},{1,1,-1,-1,x},{-1,-1,1,1,c+r},{-1,1,-1,-1,c+x},{1,-1,-1,-1,r+x},{-1,-1,-1,-1,c+r+x}};
p4m2dIrrep={{({
 {0, 1},
 {1, 0}
}),({
 {1, 0},
 {0, 1}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, -1}
}),c+x},{({
 {0, 1},
 {-1, 0}
}),({
 {-1, 0},
 {0, 1}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, -1}
}),r+x},{({
 {0, 1},
 {1, 0}
}),({
 {-1, 0},
 {0, -1}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, -1}
}),c+x},{({
 {0, 1},
 {-1, 0}
}),({
 {1, 0},
 {0, -1}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, -1}
}),r+x},{({
 {0, 1},
 {-1, 0}
}),({
 {-1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),({
 {1, 0},
 {0, 1}
}),r},{({
 {0, 1},
 {-1, 0}
}),({
 {-1, 0},
 {0, 1}
}),({
 {-1, 0},
 {0, -1}
}),({
 {-1, 0},
 {0, -1}
}),r}};
p4m4dIrrep={{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"0", "1", "0", "0"},
{"1", "0", "0", "0"},
{"0", "0", "1", "0"},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
}),({
 {1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c+r},{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "1"}], "0", "0"},
{
RowBox[{"-", "1"}], "0", "0", "0"},
{"0", "0", 
RowBox[{"-", "1"}], "0"},
{"0", "0", "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
}),({
 {1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c+r},
{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),({
 {0, 0, 0, 1},
 {0, 0, 1, 0},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, Cos[2Pi  p/n], Sin[2Pi  p/n]},
 {0, 0, -Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c},{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),({
 {0, 0, 0, -1},
 {0, 0, -1, 0},
 {0, -1, 0, 0},
 {-1, 0, 0, 0}
}),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, Cos[2Pi  p/n], Sin[2Pi  p/n]},
 {0, 0, -Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c},{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"0", "1", "0", "0"},
{"1", "0", "0", "0"},
{"0", "0", "1", "0"},
{"0", "0", "0", "1"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, -1, 0},
 {0, 0, 0, -1}
}),({
 {-1, 0, 0, 0},
 {0, -1, 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c+r},{({
 {0, 0, 1, 0},
 {0, 0, 0, 1},
 {0, 1, 0, 0},
 {1, 0, 0, 0}
}),\!\(\*
TagBox[
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", "1"}], "0", "0"},
{
RowBox[{"-", "1"}], "0", "0", "0"},
{"0", "0", 
RowBox[{"-", "1"}], "0"},
{"0", "0", "0", 
RowBox[{"-", "1"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),({
 {Cos[2Pi  p/n], Sin[2Pi  p/n], 0, 0},
 {-Sin[2Pi  p/n], Cos[2Pi  p/n], 0, 0},
 {0, 0, -1, 0},
 {0, 0, 0, -1}
}),({
 {-1, 0, 0, 0},
 {0, -1, 0, 0},
 {0, 0, Cos[2Pi  p/n], -Sin[2Pi  p/n]},
 {0, 0, Sin[2Pi  p/n], Cos[2Pi  p/n]}
}),c+r}};


(* ::Input::Initialization:: *)
p4mO2MatrixForm[i_Integer]:=Module[{a1,a2,t1,t2,m,tab},
tab=Table[0,{m,1,6}];

If[i<=136,{a1,a2}=SortNumber2[i];
t1=Quotient[a1-1,8];
t2=Quotient[a2-1,8];
a1=Mod[a1-1,8]+1;
a2=Mod[a2-1,8]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p4m1dIrrep[[a1]][[m]],p4m1dIrrep[[a2]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2}];
tab[[6]]=PolynomialMod[p4m1dIrrep[[a1]][[5]]+p4m1dIrrep[[a2]][[5]]+t*(t1+t2),2];];

If[i>136,a1=i-136;
t1=Quotient[a1-1,6];
a1=Mod[a1-1,6]+1;

For[m=1,m<=4,m++,
tab[[m]]=p4m2dIrrep[[a1]][[m]];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1}];
tab[[6]]=p4m2dIrrep[[a1]][[5]];];

tab
]

p4mO3MatrixForm[i_Integer]:=Module[{a1,a2,a3,t1,t2,t3,m,tab},
tab=Table[0,{m,1,6}];

If[i<=816,{a1,a2,a3}=SortNumber3[i];
t1=Quotient[a1-1,8];
t2=Quotient[a2-1,8];
t3=Quotient[a3-1,8];
a1=Mod[a1-1,8]+1;
a2=Mod[a2-1,8]+1;
a3=Mod[a3-1,8]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p4m1dIrrep[[a1]][[m]],p4m1dIrrep[[a2]][[m]],p4m1dIrrep[[a3]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2,(-1)^t3}];
tab[[6]]=PolynomialMod[p4m1dIrrep[[a1]][[5]]+p4m1dIrrep[[a2]][[5]]+p4m1dIrrep[[a3]][[5]]+t*(t1+t2+t3),2];];

If[i>816,a1=i-816;
a2=Quotient[a1-1,12]+1;
a1=Mod[a1-1,12]+1;

t1=Quotient[a1-1,6];
a1=Mod[a1-1,6]+1;
t2=Quotient[a2-1,8];
a2=Mod[a2-1,8]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p4m2dIrrep[[a1]][[m]],0},{0,p4m1dIrrep[[a2]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2}];
tab[[6]]=PolynomialMod[p4m2dIrrep[[a1]][[5]]+p4m1dIrrep[[a2]][[5]]+t*t2,2];];

tab
]

p4mO4MatrixForm[i_Integer]:=Module[{a1,a2,a3,a4,t1,t2,t3,t4,m,tab},
tab=Table[0,{m,1,6}];

If[i<=3876,{a1,a2,a3,a4}=SortNumber4[i];
t1=Quotient[a1-1,8];
t2=Quotient[a2-1,8];
t3=Quotient[a3-1,8];
t4=Quotient[a4-1,8];
a1=Mod[a1-1,8]+1;
a2=Mod[a2-1,8]+1;
a3=Mod[a3-1,8]+1;
a4=Mod[a4-1,8]+1;

For[m=1,m<=4,m++,
tab[[m]]=DiagonalMatrix[{p4m1dIrrep[[a1]][[m]],p4m1dIrrep[[a2]][[m]],p4m1dIrrep[[a3]][[m]],p4m1dIrrep[[a4]][[m]]}];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t2,(-1)^t3,(-1)^t4}];
tab[[6]]=PolynomialMod[p4m1dIrrep[[a1]][[5]]+p4m1dIrrep[[a2]][[5]]+p4m1dIrrep[[a3]][[5]]+p4m1dIrrep[[a4]][[5]]+t*(t1+t2+t3+t4),2];];

If[i>3876&&i<=5508,a1=i-3876;
a2=Quotient[a1-1,12]+1;
a1=Mod[a1-1,12]+1;
{a2,a3}=SortNumber2[a2];

t1=Quotient[a1-1,6];
a1=Mod[a1-1,6]+1;
t2=Quotient[a2-1,8];
a2=Mod[a2-1,8]+1;
t3=Quotient[a3-1,8];
a3=Mod[a3-1,8]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p4m2dIrrep[[a1]][[m]],0,0},{0,p4m1dIrrep[[a2]][[m]],0},{0,0,p4m1dIrrep[[a3]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2,(-1)^t3}];
tab[[6]]=PolynomialMod[p4m2dIrrep[[a1]][[5]]+p4m1dIrrep[[a2]][[5]]+p4m1dIrrep[[a3]][[5]]+t*(t2+t3),2];];

If[i>5508,a1=i-5508;
{a1,a2}=SortNumber2[a1];
t1=Quotient[a1-1,6];
t2=Quotient[a2-1,6];
a1=Mod[a1-1,6]+1;
a2=Mod[a2-1,6]+1;

For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p4m2dIrrep[[a1]][[m]],0},{0,p4m2dIrrep[[a2]][[m]]}})];];
tab[[5]]=DiagonalMatrix[{(-1)^t1,(-1)^t1,(-1)^t2,(-1)^t2}];
tab[[6]]=PolynomialMod[p4m2dIrrep[[a1]][[5]]+p4m2dIrrep[[a2]][[5]],2];];

tab
]


(* ::Input::Initialization:: *)
p4mDQCPEmbedding[i_Integer]:=Module[{O2Mat,m,w1,O3Sign,O3T,tab},
O2Mat = p4mO2MatrixForm[i];

w1=PolynomialMod[t+r+O2Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=8,m++,If[p4m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p4m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]=ArrayFlatten[({{p4m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O2Mat[[m]]}})];];
tab[[5]]=ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O2Mat[[5]]}})];

tab
]


(* ::Input::Initialization:: *)
p4mDSLEmbedding[i_Integer,j_Integer]:=Module[{O2Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O2Mat =p4mO2MatrixForm[j];
O3Mat= p4mO3MatrixForm[i];

w1=PolynomialMod[t+r+O2Mat[[6]]+O3Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=8,m++,If[p4m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p4m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p4m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O3Mat[[m]]}})],O2Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O3Mat[[5]]}})],O2Mat[[5]]};

tab
]


(* ::Input::Initialization:: *)
p4mSLEmbedding[i_Integer,j_Integer]:=Module[{O4Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O3Mat =p4mO3MatrixForm[j];
O4Mat= p4mO4MatrixForm[i];

w1=PolynomialMod[t+r+O3Mat[[6]]+O4Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=8,m++,If[p4m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p4m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p4m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O4Mat[[m]]}})],O3Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,O4Mat[[5]]}})],O3Mat[[5]]};

tab
]


(* ::Input::Initialization:: *)
p4mEmbedding[n_Integer,lsm_Integer,p_Integer]:=Module[{a1,a2,tab},
If[n==5,
a1=data[[1]][[2]][[lsm]][[p]];
tab = p4mDQCPEmbedding[a1];];
If[n==6,
a1=data[[2]][[2]][[lsm]][[p]][[1]];
a2=data[[2]][[2]][[lsm]][[p]][[2]];
tab = p4mDSLEmbedding[a1,a2];];
If[n==7,
a1=data[[3]][[2]][[lsm]][[p]][[1]];
a2=data[[3]][[2]][[lsm]][[p]][[2]];
tab = p4mSLEmbedding[a1,a2];];
tab]


p4mStability[n_Integer,lsm_Integer,p_Integer]:=Module[{check},
check=False;
If[n==5,
If[(lsm==2)||(lsm==3)||(lsm==5)||(lsm==8),check=True;];];
If[n==6,
If[lsm==2,If[(p==3)||(p==9)||(p==14),check=True;];];
If[lsm==3,If[(p==2)||(p==9)||(p==11),check=True;];];
If[lsm==4,If[(p==12)||(p==16)||(p==18)||(p==26),check=True;];];
If[(lsm==6)||(lsm==7),check=True;];
If[lsm==8,If[(p==7)||(p==10),check=True;];];];
If[n==7,
If[lsm==2,If[(p==19)||(p==30),check=True;];];
If[lsm==3,If[(p==20)||(p==105),check=True;];];
If[lsm==6,If[(p==5)||(p==14)||(p==19)||(p==37)||(p==51)||(p==54),check=True;];];
If[lsm==7,If[(p==2)||(p==5)||(p==30)||(p==38)||(p==54)||(p==64),check=True;];];];
check]


p4mPrintEmbedding[n_Integer,lsm_Integer,p_Integer]:=Module[{flag},
If[n==5,
Print["C4:n->",p4mEmbedding[5,lsm,p][[1]]//MatrixForm,"n"];
Print["M:n->",p4mEmbedding[5,lsm,p][[2]]//MatrixForm,"n"];
Print["T1:n->",p4mEmbedding[5,lsm,p][[3]]//MatrixForm,"n"];
Print["T2:n->",p4mEmbedding[5,lsm,p][[4]]//MatrixForm,"n"];
Print["T:n->",p4mEmbedding[5,lsm,p][[5]]//MatrixForm,"n"];];
If[n==6,
Print["C4:n->",p4mEmbedding[6,lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[6,lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p4mEmbedding[6,lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[6,lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p4mEmbedding[6,lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[6,lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p4mEmbedding[6,lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[6,lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p4mEmbedding[6,lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[6,lsm,p][[5]][[2]]]//MatrixForm];];
If[n==7,
Print["C4:n->",p4mEmbedding[7,lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[7,lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p4mEmbedding[7,lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[7,lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p4mEmbedding[7,lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[7,lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p4mEmbedding[7,lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[7,lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p4mEmbedding[7,lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p4mEmbedding[7,lsm,p][[5]][[2]]]//MatrixForm];];
flag=True;
flag]


p4mSL5dEmbedding[i_Integer,j_Integer]:=Module[{O2Mat,O3Mat,m,w1,O3Sign,O3T,tab},
O3Mat =p4mO3MatrixForm[j];
O2Mat= p4mO2MatrixForm[i];

w1=PolynomialMod[t+r+O3Mat[[6]]+O2Mat[[6]],2];
O3Sign=0;O3T=0;
For[m=1,m<=8,m++,If[p4m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p4m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p4m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[5],0},{0,O2Mat[[m]]}})],O3Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[5],0},{0,O2Mat[[5]]}})],O3Mat[[5]]};

tab
]


p4m5dEmbedding[lsm_Integer,p_Integer]:=Module[{a1,a2,tab},
tab=0;

a1=dataSL5Rep[[2]][[lsm]][[p]][[1]];
a2=dataSL5Rep[[2]][[lsm]][[p]][[2]];
tab = p4mSL5dEmbedding[a1,a2];

tab
]


p4m5dStability[lsm_Integer,p_Integer]:=Module[{check},
check=False;
If[lsm==1,If[(p==7)||(p==10)||(p==19)||(p==21),check=True;];];
If[(lsm==2)||(lsm==3),check=True;];
If[lsm==4,If[p==2,check=True;];];
check]


p4m5dPrintEmbedding[lsm_Integer,p_Integer]:=Module[{flag},
Print["C4:n->",p4m5dEmbedding[lsm,p][[1]][[1]]//MatrixForm,"n",Transpose[p4m5dEmbedding[lsm,p][[1]][[2]]]//MatrixForm];
Print["M:n->",p4m5dEmbedding[lsm,p][[2]][[1]]//MatrixForm,"n",Transpose[p4m5dEmbedding[lsm,p][[2]][[2]]]//MatrixForm];
Print["T1:n->",p4m5dEmbedding[lsm,p][[3]][[1]]//MatrixForm,"n",Transpose[p4m5dEmbedding[lsm,p][[3]][[2]]]//MatrixForm];
Print["T2:n->",p4m5dEmbedding[lsm,p][[4]][[1]]//MatrixForm,"n",Transpose[p4m5dEmbedding[lsm,p][[4]][[2]]]//MatrixForm];
Print["T:n->",p4m5dEmbedding[lsm,p][[5]][[1]]//MatrixForm,"n",Transpose[p4m5dEmbedding[lsm,p][[5]][[2]]]//MatrixForm];
flag=True;
flag]


(* ::Chapter:: *)
(*Incommensurate*)


(* ::Input::Initialization:: *)
IncommensurateCheck[i_Integer]:=Module[{flag1,flag2},
flag1=0;flag2=0;

If[(i==3893)||(i==5510),flag1=1;flag2=0;];
If[(i==4523)||(i==5543),flag1=1;flag2=1;];
If[(i==4025)||(i==5517),flag1=2;flag2=0;];
If[(i==4943)||(i==5562),flag1=2;flag2=1;];
If[(i==4001)||(i==4254),flag1=3;flag2=0;];
If[(i==4919)||(i==5460),flag1=3;flag2=1;];
If[(i==3929)||(i==4194),flag1=4;flag2=0;];
If[(i==4655)||(i==5304),flag1=4;flag2=1;];
If[(i==4098)||(i==5515),flag1=5;flag2=0;];
If[(i==5112)||(i==5560),flag1=5;flag2=1;];
If[(i==4290)||(i==5513),flag1=6;flag2=0;];
If[(i==5496)||(i==5552),flag1=6;flag2=1;];

{flag1,flag2}]


(* ::Input::Initialization:: *)
Incommensurate[i_Integer,j_Integer,k_Integer]:=Module[{O3Mat, O4Mat, w1,m,O3Sign,O3T,tab},
O4Mat=p4m4dIrrep[[i]];
O3Mat=p4mO3MatrixForm[k];

w1=PolynomialMod[t+r+O3Mat[[6]]+O4Mat[[5]]+t*j,2];
O3Sign=0;O3T=0;
For[m=1,m<=8,m++,If[p4m1dIrrep[[m]][[5]]==w1,O3Sign=m;];If[p4m1dIrrep[[m]][[5]]+t==w1,O3Sign=m;O3T=1;];];

tab=Table[0,{m,5}];
For[m=1,m<=4,m++,
tab[[m]]={ArrayFlatten[({{p4m1dIrrep[[O3Sign]][[m]]*IdentityMatrix[3],0},{0,O4Mat[[m]]}})],O3Mat[[m]]};];
tab[[5]]={ArrayFlatten[({{(-1)^O3T*IdentityMatrix[3],0},{0,(-1)^j*IdentityMatrix[4]}})],O3Mat[[5]]};

tab]


(* ::Input::Initialization:: *)
IncommensurateEmbedding[lsm_Integer,p_Integer]:=Module[{a1,a2,f1,f2,flag,tab},
flag=False;tab=0;

a1=data[[3]][[2]][[lsm]][[p]][[1]];
a2=data[[3]][[2]][[lsm]][[p]][[2]];
{f1,f2}=IncommensurateCheck[a1];

If[f1>0,tab = Incommensurate[f1,f2,a2];flag=True;];

{flag,tab}]


IncommensuratePrintEmbedding[lsm_Integer,p_Integer]:=Module[{flag, tab},
{flag,tab}=IncommensurateEmbedding[lsm,p];
If[flag,
Print["C4:n->",tab[[1]][[1]]//MatrixForm,"n",Transpose[tab[[1]][[2]]]//MatrixForm];
Print["M:n->",tab[[2]][[1]]//MatrixForm,"n",Transpose[tab[[2]][[2]]]//MatrixForm];
Print["T1:n->",tab[[3]][[1]]//MatrixForm,"n",Transpose[tab[[3]][[2]]]//MatrixForm];
Print["T2:n->",tab[[4]][[1]]//MatrixForm,"n",Transpose[tab[[4]][[2]]]//MatrixForm];
Print["T:n->",tab[[5]][[1]]//MatrixForm,"n",Transpose[tab[[5]][[2]]]//MatrixForm];
, Print[p, " does not lie at any incommensurate realizations"];];
flag]

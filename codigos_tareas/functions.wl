(* ::Package:: *)

(* functions.wl *)
(* Author: Jos\[EAcute] Alfredo de Le\[OAcute]n *)
(* First created: February 18, 2022 *)
(* Last update: March, 4, 2022 *)
BeginPackage["functions`"]
PivotIndex::usage="PivotIndex[] finds the pivot element's index in the submatrix \[ScriptM] of matrix M, formed by deleting rows and columns above and to the left of \!\(\*SubscriptBox[\(M\), \(i, i\)]\). PivotIndex[] implements the search of the largest matrix element of the absolute value of M. It takes into account that the input matrix M could be the augmented matrix of a linear system of equations, then the last column is left out of pivoting procedure."
Pivot::usage="Pivot[] does the pivoting in the submatrix \[ScriptM] of matrix M, where Subscript[\[ScriptM], 1,1]=Subscript[M, i,i]."
UpperDiagonalIndices::usage="Function to get all indices above the diagonal of an N \[Times] N matrix"
LowerDiagonalIndices::usage="Function to get all indices under the diagonal of an N \[Times] N matrix"
MatrixToZeroUnderPivot::usage="Function to construct the matrix that zeroes all elements under a pivoting element in position (R,R) during Gaussian elimination procedure of matrix M. The idea behind the function is to construct the appropiated matrix to multiply it by the left and make the necessary operations."
MatrixToZeroAbovePivot::usage="Function to construct the matrix that zeroes all elements above a pivoting element in position (R,R) during Gauss-Jordan elimination procedure of matrix M. The idea behind the function is to construct the appropiated matrix to multiply it by the left and make the necessary operations."
GaussianElimination::usage="Function to implement the Gaussian elimination procedure of matrix M."
BackElimination::usage="Function to implement the back elimination procedure of the uppper-triangular matrix Udummy."
GaussJordan::usage="Function to implement Gauss Jordan."
Assignationk::usage="Assignationk[M,N,n] computes the value of k for the Fock state n of N bosons and M sites"
HilbertSpaceDim::usage="HilbertSpaceDim[N, M] computes the dimension of Hilbert Space of N bosons and M sites"
FockBasis::usage="FockBasis[N, M] computes the basis of Fock states for N bosons and M sites"



Begin["`Private`"]
(* PivotIndex[] finds the pivot element's index in the submatrix \[ScriptM] of matrix M, formed by deleting rows and columns above and to the left of Subscript[M, i,i]. PivotIndex[] implements the search of the largest matrix element of the absolute value of M. It takes into account that the input matrix M could be the augmented matrix of a linear system of equations, then the last column is left out of pivoting procedure. *)
PivotIndex[M_,i_]:=(i-1)+FirstPosition[#,Max[#]]&[Abs[If[SquareMatrixQ[#],#[[i;;,i;;]],#[[i;;,i;;-2]]]&[M]]]


(* Pivot[] does the pivoting in the submatrix \[ScriptM] of matrix M, where Subscript[\[ScriptM], 1,1]=Subscript[M, i,i]. *)
Pivot[M_,i_,colSwaps_]:=Module[{pivotIndex,Mdummy,colSwps},Mdummy=M;
pivotIndex=PivotIndex[M,i];colSwps=colSwaps;
(* Interchange rows *)
Mdummy[[{i,pivotIndex[[1]]}]]=Mdummy[[{pivotIndex[[1]],i}]];
(* Interchange columns *)
Mdummy[[All,{i,pivotIndex[[2]]}]]=Mdummy[[All,{pivotIndex[[2]],i}]];
(* Multiply ith row by 1/a, where a is the pivot element *)
Mdummy[[{i}]]*=1/Mdummy[[i,i]];
(* Record of the column interchanges *)
colSwps[[{i,pivotIndex[[2]]}]]=colSwps[[{pivotIndex[[2]],i}]];
{Mdummy,colSwps}]


(* Function to get all indices above the diagonal of an N \[Times] N matrix *)
UpperDiagonalIndices[N_]:=Table[{i,j},{j,Range[N]},{i,j-1}]


(* Function to get all indices under the diagonal of an N \[Times] N matrix *)
LowerDiagonalIndices[N_]:=Table[{i,j},{j,N},{i,Range[j+1,N]}]


(* Function to construct the matrix that zeroes all elements under a pivoting element in position (R,R) during Gaussian elimination procedure of matrix M. The idea behind the function is to construct the appropiated matrix to multiply it by the left and make the necessary operations. *)
MatrixToZeroUnderPivot[M_,R_]:=Module[{lowDiagInd,N},
N=Length[M];
lowDiagInd=LowerDiagonalIndices[N];
IdentityMatrix[N]+SparseArray[lowDiagInd[[R]]->(-M[[#[[1]],#[[2]]]]&/@lowDiagInd[[R]]),{N,N}]
]


(* Function to construct the matrix that zeroes all elements above a pivoting element in position (R,R) during Gauss-Jordan elimination procedure of matrix M. The idea behind the function is to construct the appropiated matrix to multiply it by the left and make the necessary operations. *)
MatrixToZeroAbovePivot[U_,R_]:=Module[{uppDiagInd,N},
N=Length[U];
uppDiagInd=UpperDiagonalIndices[N];
IdentityMatrix[N]+SparseArray[uppDiagInd[[R]]->(-U[[#[[1]],#[[2]]]]&/@uppDiagInd[[R]]),{N,N}]
]


(* Function to implement the Gaussian elimination procedure of matrix M *)
GaussianElimination[M_]:=Module[{N,Mdummy,pivotIndex,colSwaps},
Mdummy=M;N=Length[Mdummy];
(* Create a list to keep record of column swaps. Consider that the matrix M may not be square. *)
colSwaps=Range[If[SquareMatrixQ[#],N,Dimensions[#][[2]]-1]&[Mdummy]];
(* If the submatrix formed by deleting the rows and columns above and to the left of position (i,i) is not the zero matrix, then call Pivot[] and make zeroes under the pivot element. *)
Table[
If[AnyTrue[Flatten[If[SquareMatrixQ[Mdummy],Mdummy[[i;;,i;;]],Mdummy[[i;;,i;;-2]]]],#!=0&],
{Mdummy,colSwaps}=Pivot[Mdummy,i,colSwaps];
Mdummy[[i,i]]=Round[Mdummy[[i,i]]];
Mdummy=MatrixToZeroUnderPivot[Mdummy,i] . Mdummy;];
,{i,N}];
{Mdummy,colSwaps}
]

(* Function to implement the back elimination procedure of the uppper-triangular matrix Udummy *)
BackElimination[Udummy_]:=Module[{N,U},
U=Udummy;
(* Find the last 1 on the diagonal *)
N=If[AnyTrue[#,#==0&],FromDigits[FirstPosition[#,0]]-1,Length[U]]&[Table[U[[i,i]],{i,Length[U]}]];
(* Iterate to make zeroes above diagonal 1's. *)
Table[U=MatrixToZeroAbovePivot[U,i] . U;
,{i,2,N}];
U
]


(* Function to implement Gauss Jordan *)
GaussJordan[M_]:={BackElimination[#[[1]]],#[[2]]}&[GaussianElimination[M]]


(* Assignationk[M,N,n] computes the value of k for the Fock state n of N bosons and M sites *)
Assignationk[M_,N_,fockState_]:=If[fockState[[1;;M-1]]==ConstantArray[0,M-1],M-1,FromDigits[Last[Position[Normal[fockState[[1;;M-1]]],x_ /;x!=0]]]]


(* HilbertSpaceDim[N, M] computes the dimension of Hilbert Space of N bosons and M sites *)
HilbertSpaceDim[N_,M_]:=(N+M-1)!/(N!(M-1)!)


(* FockBasis[N, M] computes the basis of Fock states for N bosons and M sites *)
FockBasis[N_,M_]:=Module[{k,fockState},
k=1;
Join[
(* Compute the first lexycographical Fock state *)
{fockState=SparseArray[{1->N},{N}]},
(* Construct the following elements of Fock basis *)
Table[
(* With \[Eta] the new Fock state and n the previous one, assign Subscript[\[Eta], i]=Subscript[n, i] (1<=i<=k-1), Subscript[\[Eta], k]=Subscript[n, k]-1 y Subscript[\[Eta], i]=0 (i>=k+2) *)
fockState=SparseArray[Join[Table[i->fockState[[i]],{i,k-1}],{k->fockState[[k]]-1}],{N}];
(* With \[Eta] the new Fock state and n the previous one, assign Subscript[\[Eta], k+1]=N-\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(k\)]
\*SubscriptBox[\(\[Eta]\), \(i\)]\) *)
fockState[[k+1]]=N-Total[fockState[[1;;k]]];
(* Compute next value of k *)
k=Assignationk[M,N,fockState];
fockState
,HilbertSpaceDim[N,M]-1]
]
]


End[];
EndPackage[]




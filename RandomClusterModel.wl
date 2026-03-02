(* ::Package:: *)

ClearAll["Global`*"]
LaunchKernels[$ProcessorCount];


(*ZacharyKarateClub (Monte-Carlo simulation).*)
graph=ExampleData[{"NetworkGraph","ZacharyKarateClub"}];
vertexList=VertexList[graph];
vertexListOrdering=AssociationThread[vertexList->Range[Length@vertexList]];
edgeList=EdgeList[graph];


maxIterations=1000000;
\[CapitalDelta]p=0.05;\[CapitalDelta]q=0.05;


Print["Time: ", AbsoluteTiming[dataRCConnectedness=Table[
Module[{dataRC=ParallelTable[
Module[{subGraph=Graph[vertexList,Pick[edgeList,Thread[RandomReal[{0,1},Length@edgeList]<p]]],
connectedComponents},
connectedComponents=ConnectedComponents[subGraph]/.vertexListOrdering;
{SparseArray[Thread[Flatten[Tuples[#,2]&/@connectedComponents,1]->1],{Length@vertexList,Length@vertexList}],Length@connectedComponents}
],
maxIterations],dataRC2},
dataRC2=dataRC[[;;,2]]-Min[dataRC[[;;,2]]];
Table[
Module[{weightsRC=q^dataRC2,totalWeightsRC},
totalWeightsRC=Total[weightsRC];
weightsRC . dataRC[[;;,1]]/totalWeightsRC
]
,{q,1+\[CapitalDelta]q,2.,\[CapitalDelta]q}]

]
,{p,\[CapitalDelta]p,1-\[CapitalDelta]p,\[CapitalDelta]p}];
]]


DumpSave[FileNameJoin[{NotebookDirectory[],"RandomClusterModel_karate-club_"<>ToString[maxIterations]<>"_p={0.05,0.95,0.05}_q={1.05,2.,0.05}.mx"}],dataRCConnectedness];






(*ZacharyKarateClub (this the theoretical exact result using the Tutte polynomial - only works on small graphs).*)


(*Tutte polynomial.*)
graph=ExampleData[{"NetworkGraph","ZacharyKarateClub"}];
edgeList=EdgeList@AdjacencyGraph[AdjacencyMatrix[graph]];
graph=Graph[Range[VertexCount[edgeList]],edgeList];
partition[edgeList_,p_,q_]:=q p^(VertexCount[edgeList]-1) (1-p)^(EdgeCount[edgeList]-VertexCount[edgeList]+1) TuttePolynomial[edgeList,{1+(q(1-p))/p,1+p/(1-p)}];
z0=partition[edgeList,p,q];
\[Mu]Connectedness[edgeList_,p_,q_,i_,j_]:=Module[{zContract=partition[edgeList/.{i->j},p,q]},(*FullSimplify[1+q/(q-1)((1-p)zContract-(1-p)(z0-p zContract)/(1-p))/z0]*)(1-q zContract/z0)/(1-q)];


\[Mu]ConnectednessData=ParallelTable[If[i>=j,0,\[Mu]Connectedness[edgeList,p,q,i,j]],{i,1,VertexCount[graph]},{j,1,VertexCount[graph]}];
(*The output is a matrix (row: node i; column: node j). Each entry is a polynomial of p and q.*)
DumpSave[FileNameJoin[{NotebookDirectory[],"RandomClusterModel_karate-club_theoretical.mx"}],\[Mu]ConnectednessData];






(* ::Package:: *)

ClearAll["Global`*"]

(* --- Helper: resolve output directory (notebook or script) --- *)
outputDirectory = Quiet[Check[NotebookDirectory[], $Failed]];
If[FailureQ[outputDirectory] || !DirectoryQ[outputDirectory],
  outputDirectory = Directory[];
  Print["Warning: NotebookDirectory[] unavailable; falling back to working directory: ", outputDirectory];
];

(* --- Launch parallel kernels with validation --- *)
Module[{launched},
  launched = LaunchKernels[$ProcessorCount];
  If[launched === $Failed || Length[Kernels[]] === 0,
    Print["Error: Failed to launch parallel kernels. Aborting."];
    Abort[];
  ];
  Print["Launched ", Length[Kernels[]], " parallel kernels."];
];


(*ZacharyKarateClub (Monte-Carlo simulation).*)
graph = ExampleData[{"NetworkGraph", "ZacharyKarateClub"}];
If[!GraphQ[graph],
  Print["Error: Failed to load ZacharyKarateClub graph from ExampleData. Aborting."];
  Abort[];
];
vertexList = VertexList[graph];
vertexListOrdering = AssociationThread[vertexList -> Range[Length@vertexList]];
edgeList = EdgeList[graph];


maxIterations = 1000000;
\[CapitalDelta]p = 0.05; \[CapitalDelta]q = 0.05;


Print["Time: ", AbsoluteTiming[dataRCConnectedness = Table[
  Module[{dataRC = ParallelTable[
    Module[{subGraph = Graph[vertexList, Pick[edgeList, Thread[RandomReal[{0, 1}, Length@edgeList] < p]]],
      connectedComponents},
      connectedComponents = ConnectedComponents[subGraph] /. vertexListOrdering;
      {SparseArray[Thread[Flatten[Tuples[#, 2] & /@ connectedComponents, 1] -> 1], {Length@vertexList, Length@vertexList}], Length@connectedComponents}
    ],
    maxIterations], dataRC2},

    (* Validate parallel results before proceeding *)
    If[Length[dataRC] === 0,
      Print["Error: ParallelTable returned empty results for p=", p, ". Skipping."];
      Nothing,
      dataRC2 = dataRC[[;; , 2]] - Min[dataRC[[;; , 2]]];
      Table[
        Module[{weightsRC = q^dataRC2, totalWeightsRC},
          totalWeightsRC = Total[weightsRC];
          If[totalWeightsRC == 0 || !NumericQ[totalWeightsRC],
            Print["Warning: Zero or non-numeric total weight for p=", p, ", q=", q, ". Returning zero matrix."];
            SparseArray[{}, {Length@vertexList, Length@vertexList}],
            weightsRC . dataRC[[;; , 1]] / totalWeightsRC
          ]
        ]
        , {q, 1 + \[CapitalDelta]q, 2., \[CapitalDelta]q}]
    ]
  ]
  , {p, \[CapitalDelta]p, 1 - \[CapitalDelta]p, \[CapitalDelta]p}];
]]


(* --- Save Monte Carlo results --- *)
Module[{mcPath, saveResult},
  mcPath = FileNameJoin[{outputDirectory, "RandomClusterModel_karate-club_" <> ToString[maxIterations] <> "_p={0.05,0.95,0.05}_q={1.05,2.,0.05}.mx"}];
  saveResult = DumpSave[mcPath, dataRCConnectedness];
  If[saveResult === $Failed || saveResult === {},
    Print["Error: Failed to save Monte Carlo results to ", mcPath];,
    Print["Saved Monte Carlo results to ", mcPath];
  ];
];


(* ===================================================================== *)
(*ZacharyKarateClub (this the theoretical exact result using the Tutte polynomial - only works on small graphs).*)


(*Tutte polynomial.*)
graph = ExampleData[{"NetworkGraph", "ZacharyKarateClub"}];
If[!GraphQ[graph],
  Print["Error: Failed to reload ZacharyKarateClub graph. Aborting theoretical calculation."];
  Abort[];
];
edgeList = EdgeList@AdjacencyGraph[AdjacencyMatrix[graph]];
graph = Graph[Range[VertexCount[graph]], edgeList];

partition[g_, p_, q_] := Module[{vc, ec},
  vc = VertexCount[g];
  ec = EdgeCount[g];
  If[vc === 0,
    Print["Error: partition called on graph with 0 vertices."];
    Return[$Failed];
  ];
  q p^(vc - 1) (1 - p)^(ec - vc + 1) TuttePolynomial[g, {1 + (q (1 - p)) / p, 1 + p / (1 - p)}]
];

z0 = partition[graph, p, q];
If[z0 === $Failed,
  Print["Error: Failed to compute partition function z0. Aborting."];
  Abort[];
];

\[Mu]Connectedness[g_, p_, q_, i_, j_] := Module[{contracted, zContract},
  contracted = g /. {i -> j};
  zContract = partition[contracted, p, q];
  If[zContract === $Failed,
    Print["Error: partition failed for contracted graph (i=", i, ", j=", j, ")."];
    Return[$Failed];
  ];
  (* The (1-q) denominator is removable via L'Hopital at q=1, but flag if hit numerically *)
  (1 - q zContract / z0) / (1 - q)
];


\[Mu]ConnectednessData = ParallelTable[
  If[i >= j, 0,
    Module[{result},
      result = \[Mu]Connectedness[graph, p, q, i, j];
      If[result === $Failed,
        Print["Warning: Connectedness computation failed for (i=", i, ", j=", j, "). Using 0."];
        0,
        result
      ]
    ]
  ],
  {i, 1, VertexCount[graph]}, {j, 1, VertexCount[graph]}
];
(*The output is a matrix (row: node i; column: node j). Each entry is a polynomial of p and q.*)

(* --- Save theoretical results --- *)
Module[{theoPath, saveResult},
  theoPath = FileNameJoin[{outputDirectory, "RandomClusterModel_karate-club_theoretical.mx"}];
  saveResult = DumpSave[theoPath, \[Mu]ConnectednessData];
  If[saveResult === $Failed || saveResult === {},
    Print["Error: Failed to save theoretical results to ", theoPath];,
    Print["Saved theoretical results to ", theoPath];
  ];
];

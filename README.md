# Random Cluster Model

This repository contains Mathematica code to calculate the pairwise Random Cluster (RC) connectedness of nodes in, e.g., the Zachary Karate Club network under the RC model. 

The codebase provides two distinct approaches to evaluating the connectedness probability between nodes: a high-throughput Monte Carlo simulation for numerical estimation and an exact analytical calculation using the Tutte polynomial. 

## Overview 

The Random Cluster model is a generalization of percolation and the Potts model. It assigns a probability weight *p* to edges and a parameter *q* to the number of connected components. This makes it an excellent framework for studying phase transitions and cluster formations in complex networks.

This script executes the model on the standard Zachary Karate Club graph (34 nodes) and outputs `.mx` files containing the connectedness matrices.

### 1. Monte Carlo Simulation
Because exact calculations scale poorly with network size, the first half of the code uses a highly parallelized Monte Carlo approach. 

* **Parameters:** It sweeps through edge probabilities *p* and cluster weights *q* in increments of Δp = 0.05 and Δq = 0.05. 
* **Method:** It generates 1,000,000 subgraph configurations per *p*-value. To maintain numerical stability when calculating the statistical weights (*q* raised to the number of components), it shifts the exponent relative to the minimum number of components observed.
* **Output:** Saves an `.mx` file containing numerical matrices of the pairwise connectedness probabilities.

### 2. Exact Theoretical Calculation
For small graphs, it is possible to compute the exact partition function of the RC model by mapping it to the Tutte polynomial.

* **Method:** The script evaluates the exact partition function using Mathematica's built-in `TuttePolynomial` function. The pairwise connectedness is then determined exactly by contracting the edge between node *i* and node *j* and computing the ratio of the contracted partition function to the full partition function.
* **Output:** Saves an `.mx` file containing a symmetric matrix where each entry (*i*, *j*) is a symbolic polynomial of *p* and *q* representing the exact connectedness probability.

## Prerequisites

* **Mathematica:** The code relies on built-in Wolfram Language functions (`TuttePolynomial`, `ConnectedComponents`, `ParallelTable`).
* **Multi-core Processor:** The script automatically launches kernels equal to your `$ProcessorCount` to run the Monte Carlo loop and the pairwise matrix generation in parallel.

## Usage

1. Open the file in Mathematica.
2. Evaluate the notebook. 
3. The script will automatically fetch the Zachary Karate Club graph via `ExampleData`.
4. The simulation will run and drop two data dumps into the same directory as the notebook:
   * `RandomClusterModel_karate-club_1000000_p={0.05,0.95,0.05}_q={1.05,2.,0.05}.mx`
   * `RandomClusterModel_karate-club_theoretical.mx`

*Note: Calculating the Tutte polynomial is #P-hard. While it evaluates relatively quickly on the 34-node Karate Club graph, this theoretical exact method will bottleneck significantly if applied to much larger empirical networks.*

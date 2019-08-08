# LaplaceGillespie
C/C++ codes for the Laplace Gillespie algorithm, one to generate positively correlated inter-event times etc.

Usage:

1. Download all codes in the same folder.

2. Download [mt19937ar.c](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html)
in the same folder.
 
3. Compile and run.

NB: corr-LGA.cc and nMGA-LGA-no-binary-tree.cc just generate event sequences as done in the first part of the numerical simulations in the paper. sir-node-centric-gillespie.cc carries out simulations of the SIR model in either the well-mixed population (i.e., complete graph) or given networks. So, most users would be interested in using sir-node-centric-gillespie.cc or its modifications.

Original paper: [Naoki Masuda & Luis E. C. Rocha, SIAM Review, 60, 95-115 (2018).](https://epubs.siam.org/doi/10.1137/16M1055876)

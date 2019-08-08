/* Laplace Gillespie algorithm implemented using binary tree data structure

   Copyright (C) 2016, Naoki Masuda, All rights reserved.                          Please cite the following paper when you use this code.
   Naoki Masuda, Luis E. C. Rocha. A Gillespie algorithm for non-Markovian stochastic processes. arXiv (2016).

   A power-law distribution of inter-event times assumed.
   p(lambda) is distributed according to a gamma distribution.

   Usage: a.out Np

     Np: number of renewal processes running in parallel

   -DNOOUT: suppress outputting the distribution of inter-event times to perform speed comparison (Fig. 1(e)).
   If NOOUT is undefined, the distribution of inter-event times is passed to cout (used in Fig. 1(a)-1(d)). 
*/
#include <iostream>
using namespace std;
#include <cstdlib> // atof
#include <cstring>
#include <ctime>

#include "mt19937ar.c"
/* Mersenne Twister to generate random variates on [0,1].
   mt19937ar.c is available at 
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c

   If one prefers to use the in-built random number generator,
   (genrand_int32()+0.5)/4294967296.0 should be replaced by (double)rand()/RAND_MAX
   and
   init_genrand(time(NULL)) should be replaced by, e.g., srand(time(NULL)) */

#include "gamma_marsaglia_tsang.cc"
/* Generate gamma distributed random variables using the Marsaglia-Tsang algorithm. It requires mt19937ar.c */

int compare_gillespie(const void *a, const void *b); // for qsort

int main(int argc, char **argv) {

  if (argc != 2) {
    cerr << "LGA-binary-tree.out Np" << endl;
    cerr << "Np: number of processes running in parallel" << endl;
    exit(8);
  }

  init_genrand(time(NULL)); // initialize the random number generator
  
  int i,j,ii; // counters
  double alpha[3]; // shape parameter of the gamma distribution, psi(tau) = alpha/(1+tau)^{alpha+1}
  alpha[0] = 1.0;
  alpha[1] = 1.5;
  alpha[2] = 2.0;
  int Np = atoi(argv[1]); // # processes

  // preparation to construct a binary tree
  int log2_Np_upper = (int)(log(Np-1e-6)/log(2)) + 1;
  // 2^{log2_Np_upper-1} < Np <= 2^{log2_Np_upper}
  int pow2[log2_Np_upper+1]; // power[i] = 2^i.
  pow2[0]=1;
  for (i=0 ; i<log2_Np_upper ; i++)
    pow2[i+1] = 2 * pow2[i];

  // event rates and the binary tree
  double lambda[2*pow2[log2_Np_upper]]; // rate of Poisson processes and associated binary tree
  double Dlambda; // working var for updating the binary tree
  int start[log2_Np_upper+1]; // start[i]=j if the i-th layer of the binary tree starts with lambda[j]

  /* --- layer 0 ---
   lambda[0], ..., lambda[Np-1] contain original lambda.
   lambda[Np], ..., lambda[2^{log2_Np_upper}-1] are set to 0 (dummy).
   These elements of lambda[] defines layer 0, so start[0] = 0.

   --- layer 1 ---  
   lambda[2^{log2_Np_upper}], .., lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}-1] contain the sum of two adjacent lambda[i]'s in layer 0,
   where 0 <= i < 2^{log2_Np_upper}.
   Therefore, start[1] = 2^{log2_Np_upper}.

   --- layer 2 ---
   lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}], ..., lambda[2^{log2_Np_upper}+2^{log2_Np_upper-1}+2^{log2_Np_upper-2}-1] contain the sum of two adjacent lambda[i]'s in layer 1,
   where 2^{log2_Np_upper} <= i < 2^{log2_Np_upper}+2^{log2_Np_upper-1}).
   Therefore, start[2] = 2^{log2_Np_upper}+2^{log2_Np_upper-1}.
   Similar for layers 3, 4, ... */
  start[0] = 0;
  for (i=0 ; i<log2_Np_upper ; i++)
    start[i+1] = start[i] + pow2[log2_Np_upper-i];
  for (i=0 ; i<Np ; i++)
    lambda[i] = rgamma(alpha[i%3],1); // gamma-distributed random var
  for (i=Np ; i<2*pow2[log2_Np_upper] ; i++)
    lambda[i] = 0.0; // dummy
  for (i = 0 ; i < log2_Np_upper ; i++) { // construct binary tree
    for (j = 0 ; j < pow2[log2_Np_upper-i] ; j++)
      lambda[start[i+1]+j/2] += lambda[start[i]+j];
  }

  double t_prev[Np]; // time of the previous event for each process
  int Nevent[Np]; // # events
  int Np_rec = (Np > 3)? 3 : Np; // # processes whose inter-event times are recorded (Fig. 1(a)-1(d))
  int Ntau = 1000000; // # inter-event times generated for one process
#ifndef NOOUT  
  double **tau = new double*[Np_rec]; // inter-event time for the first Np_rec processes
  for (i=0 ; i<Np_rec ; i++)
    tau[i] = new double[Ntau]; // initialization
#endif // ifndef NOOUT
  // initialization
  for (i=0 ; i<Np ; i++) {
    t_prev[i] = 0.0;
    Nevent[i] = 0;
  }

  double t = 0.0;
  double dt;
  double ra;
  i = 0; // initialization (dummy)

  // renewal process dynamics start here
  while (Nevent[i] < Ntau) { // till any process generates Ntau events

    // determine the process that has fired an event
    ra = (genrand_int32()+0.5)/4294967296.0 * lambda[2*pow2[log2_Np_upper]-2];
    // ra is uniformly distributed between 0 and lambda[2*pow2[log2_Np_upper]-2]
    i = 0;
    for (j = log2_Np_upper-1 ; j >= 0 ; j--) { // fast search using binary tree
      i *= 2;
      if (ra > lambda[start[j]+i]) {
	ra -= lambda[start[j]+i]; 
	i++;
      }
    }
    if (i >= Np) {
      cerr << "selection of the process failed" << endl;
      exit(8);
    }
    // the process that fired the event (i.e., i) determined

    // spend another uniformly generated random variate 
    dt = -log((genrand_int32()+0.5)/4294967296.0)/lambda[2*pow2[log2_Np_upper]-2];

    // avoid spending another random variate (Appendix A)
    // dt = -log(ra/lambda[i])/lambda[2*pow2[log2_Np_upper]-2];

    t += dt; // advance time
   
#ifndef NOOUT // record the inter-event time
    if (i < Np_rec)
      tau[i][Nevent[i]] = t - t_prev[i];
#endif // ifndef NOOUT
    t_prev[i] = t;
    Nevent[i]++;

    // redraw the event rate for the i-th process and update the binary tree
    Dlambda = lambda[i];
    lambda[i] = rgamma(alpha[i%3],1);
    Dlambda = lambda[i] - Dlambda; // increment in lambda[i]
    ii = i;
    for (j = 0 ; j < log2_Np_upper ; j++) {
      ii /= 2;
      lambda[start[j+1]+ii] += Dlambda; // update the binary tree
    }

  } // renewal process dynamics completed

// output the distribution of inter-event times for the first Np_rec processes
#ifndef NOOUT
  for (i=0 ; i<Np_rec ; i++) {
    // create a rank plot, which is equivalent to the survival function
    qsort(tau[i],Nevent[i],sizeof(double),compare_gillespie);
    for (j=0 ; j<Nevent[i] ; j++)
      if ( (j<0.99*Nevent[i] && j%(Nevent[i]/1000)==0) || (j>=0.99*Nevent[i] && j<0.999*Nevent[i] && j%(Nevent[i]/10000)==0) || (j>=0.999*Nevent[i]) )
	// avoid plotting all points to reduce the figure size
	cout << tau[i][j] << " " << 1.0 - (double)(j+0.5)/Nevent[i] << endl;
    cout << endl;
  }
  for (i=0 ; i<Np_rec ; i++)
    delete[] tau[i];
  delete[] tau;
#endif // ifndef NOOUT

  return 0;
}

// for qsort
int compare_gillespie(const void *a, const void *b) {
  if (*(double*)a<*(double*)b)
    return (-1);
  else if (*(double*)a>*(double*)b)
    return (1);
  return (0);
}

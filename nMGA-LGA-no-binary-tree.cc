/* non-Markovian Gillespie algorithm (nMGA) and Laplacian Gillespie algorithm without using binary tree data structure

   Copyright (C) 2016, Naoki Masuda, All rights reserved.                          Please cite the following paper when you use this code.
   Naoki Masuda, Luis E. C. Rocha. A Gillespie algorithm for non-Markovian stochastic processes. arXiv (2016).

   The process i procuding each event is searched linearly along array lambda[]
   A power-law distribution of inter-event times assumed.
   p(lambda) is distributed according to a gamma distribution.

   Usage: a.out algorithm Np

     algorithm = 0 for nMGA and = 1 for the Laplace Gillespie algorithm
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

#include "gamma.c"
/* gamma.c is necessary for using function rgamma(,) to generate gamma variates
   Save the code avilable at
   http://www.know-all.net/articles/view/56
   as gamma.c and delete the main() block.
   Note: gamma.c requires mt19937ar.c */

int compare_gillespie(const void *a, const void *b); // for qsort

int main(int argc, char **argv) {

  if (argc != 3) {
    cerr << "nMGA-LGA-no-binary-tree.out algorithm Np" << endl;
    cerr << "algorithm = 0: nMGA, 1: Laplace Gillespie" << endl;
    cerr << "Np: number of processes running in parallel" << endl;
    exit(8);
  }

  init_genrand(time(NULL)); // initialize the random number generator
  
  int algorithm = atoi(argv[1]);
  if (algorithm <= -1 || algorithm >= 2) {
    cerr << "algorithm must be 0 or 1" << endl;
    cerr << "algorithm = 0: nMGA, 1: Laplace Gillespie" << endl;
    exit(8);
  }
  
  int i,j; // counters
  double alpha[3]; // shape parameter of the gamma distribution, psi(tau) = alpha/(1+tau)^{alpha+1}
  alpha[0] = 1.0;
  alpha[1] = 1.5;
  alpha[2] = 2.0;
  int Np = atoi(argv[2]); // # processes

  double lambda[Np]; // rate of Poisson processes
  double sum_lambda = 0.0;

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
    if (algorithm==1) { // Laplace Gillespie
      lambda[i] = rgamma(alpha[i%3],1);
      sum_lambda += lambda[i];
    }
    t_prev[i] = 0.0;
    Nevent[i] = 0;
  }
  
  double t = 0.0;
  double dt;
  double ra;
  i = 0; // initialization (dummy)

  // renewal process dynamics start here
  while (Nevent[i] < Ntau) { // till any process generates Ntau events

    if (algorithm==0) { // initialization for nMGA
      sum_lambda = 0.0;
      for (j=0 ; j<Np ; j++) {
	lambda[j] = alpha[j%3]/(1+t-t_prev[j]);
	sum_lambda += lambda[j];
      }
    }

    // determine the process that has fired an event
    ra = (genrand_int32()+0.5)/4294967296.0 * sum_lambda; // 0 < ra < 1
    // ra is uniformly distributed between 0 and sum_lambda
    i = 0;
    while (ra > lambda[i]) { // linear search
      ra -= lambda[i];
      i++;
    }
    if (i >= Np) {
      cerr << "selection of the process failed" << endl;
      exit(8);
    }
    // the process that fired the event (i.e., i) determined

    // spend another uniformly generated random variate 
    dt = -log((genrand_int32()+0.5)/4294967296.0)/sum_lambda;

    // avoid spending another random variate (Appendix A)
    // dt = -log(ra/lambda[i])/sum_lambda;

    t += dt; // advance time

#ifndef NOOUT // record the inter-event time
    if (i < Np_rec)
      tau[i][Nevent[i]] = t - t_prev[i];
#endif // ifndef NOOUT
    t_prev[i] = t;
    Nevent[i]++;
      
    // draw a new rate and update sum_lambda
    if (algorithm==1) {
      sum_lambda -= lambda[i];
      lambda[i] = rgamma(alpha[i%3],1);
      sum_lambda += lambda[i];
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

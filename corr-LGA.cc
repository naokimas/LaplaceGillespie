/* correlated Laplace Gillespie algorithm to generate renewal processes with correlated inter-event times 

   Copyright (C) 2016, Naoki Masuda, All rights reserved.                          Please cite the following paper when you use this code.
   Naoki Masuda, Luis E. C. Rocha. A Gillespie algorithm for non-Markovian stochastic processes. arXiv (2016).

   A power-law distribution of inter-event times assumed.
   p(lambda) is distributed according to a gamma distribution.
   Run a single renewal process.

   Usage: a.out alpha p
     alpha: shape parameter of the gamma distribution
     p: prob that the rate of a Poisson process is not redrawn. p=0 corresponds to the uncorrelated Laplace Gillespie algorithm

   -DSURV: output (passed to cout) the survival function for a single run of a single process (Fig. 2(c), (d))
   If SURV is undefined, the memory coefficient (avg +- std) is calculated (Fig. 2(a), (b)) and passed to cout.
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

  if (argc != 3) {
    cerr << "corr-LGA.out alpha p" << endl;
    cerr << "alpha: shape parameter of the gamma distribution of lambda" << endl;
    cerr << "p: probability that lambda is not redrawn upon an event" << endl;
    exit(8);
  }

  init_genrand(time(NULL)); // initialize the random number generator

  int i,j; // counters
  double alpha = atof(argv[1]); // shape parameter of the gamma distribution, psi(tau) = alpha/(1+tau)^{alpha+1}
  int trials = 1000; // # samples from which we calculate the mean and std of the memory coefficient
  // Different from typical usage of Gillespie algorithms, renewal processes run one after another, not in parallel.

  double lambda; // rate of the Poisson process

  double t_prev; // time of the previous event
  int Nevent; // # events
  int Ntau = 100000; // # inter-event times generated for one event sequence
#ifdef SURV // calculate the survival function Psi(tau)
  Ntau = 1000000;
  trials = 1;
#endif // SURV
  double tau[Ntau]; // inter-event time
  double ave_tau[2],std_tau[2];
  
  double t;
  double dt;

  double p = atof(argv[2]); // probability to use the same exp dist
  if (p > 1 || p < 0) {
    cerr << "0 <= p <= 1 violated" << endl;
    exit(8);
  }

  double MC; // memory coefficient
#ifndef SURV // if SURV is not defined
  double MC_ave = 0.0;
  double MC_std = 0.0;
#endif // ifndef SURV
  
  for (i=0 ; i<trials ; i++) {

    // initialization
    lambda = rgamma(alpha,1);
    t = t_prev = 0.0;
    Nevent = 0;

    // renewal process dynamics start here
    while (Nevent < Ntau) {

      dt = -log((genrand_int32()+0.5)/4294967296.0)/lambda;
      t += dt;
      tau[Nevent] = t - t_prev; // inter-event time
      t_prev = t;
      Nevent++;

      if ((genrand_int32()+0.5)/4294967296.0 < 1-p) // select a new rate lambda with probability 1-p
	lambda = rgamma(alpha,1);
      
    } // renewal process dynamics completed

    // calculate the memory coefficient (Goh & Barabasi, EPL, 2008)
    ave_tau[0] = ave_tau[1] = std_tau[0] = std_tau[1] = MC = 0.0;
    for (j=0 ; j<Nevent-1 ; j++) {
      ave_tau[0] += tau[j];
      std_tau[0] += tau[j]*tau[j];
      ave_tau[1] += tau[j+1];
      std_tau[1] += tau[j+1]*tau[j+1];
      MC += tau[j]*tau[j+1];
    }
    for (j=0 ; j<2 ; j++) {
      ave_tau[j] /= Nevent-1;
      std_tau[j] = sqrt(std_tau[j]/(Nevent-1) - ave_tau[j]*ave_tau[j]);
    }
    MC = (MC/(Nevent-1) - ave_tau[0]*ave_tau[1])/std_tau[0]/std_tau[1];

// output the distribution of inter-event times
#ifdef SURV
    // create a rank plot, which is equivalent to the survival function
    qsort(tau,Nevent,sizeof(double),compare_gillespie);
    for (j=0 ; j<Nevent ; j++)
      if ( (j<0.99*Nevent && j%(Nevent/1000)==0) || (j>=0.99*Nevent && j<0.999*Nevent && j%(Nevent/10000)==0) || (j>=0.999*Nevent) )
	// avoid plotting all points to reduce the figure size
	cout << tau[j] << " " << 1.0 - (double)(j+0.5)/Nevent << endl;
    cout << endl;
    cerr << "memory coef = " << MC << endl;
#else // if SURV is not defined
    MC_ave += MC;
    MC_std += MC*MC;
#endif // SURV

  } // all trials done
  
// average and standard deviation of the memory coefficient
#ifndef SURV // if SURV is not defined
  MC_ave /= trials;
  MC_std = sqrt(MC_std/trials - MC_ave*MC_ave);
  cout << p << " " << MC_ave << " " << MC_std << endl;
  cerr << "memory coef = " << MC_ave << " +- " << MC_std << endl;
#endif // ifndef SURV
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

/* Numerical simulations of SIR on networks using the Laplace Gillespie algorithm 

   Copyright (C) 2016, Naoki Masuda, All rights reserved.                          Please cite the following paper when you use this code.
   Naoki Masuda, Luis E. C. Rocha. A Gillespie algorithm for non-Markovian stochastic processes. arXiv (2016).

   The initial infected node is always the first node.

   Implemented using binary tree data structure.

   A power-law distribution of inter-event times assumed.
   p(lambda) is distributed according to a gamma distribution.
   However, if alpha<0, the exponential distribution (i.e., Poisson process)
   is assumed.

   Usage: a.out nV alpha p (if MIX is defined)
        : a.out infilename alpha p (if MIX is undefined)

     nV: number of nodes in the complete graph
     infilename: name of the file having a list of links
     alpha: shape parameter of the gamma distribution
     p: prob that the rate of a Poisson process is not redrawn. p=0 corresponds to the uncorrelated Laplace Gillespie algorithm

   Output (passed to cout): averaged final size for a range of effective infectious rate values

  -DMIX: well-mixed population (= complete graph).
  If this option is off, read an external file having a list of links such that each row of the file has two endpoint nodes in the first and second columns and the link weight (=1) in the third column. In the input file, the node index must range from 1 to nV, not from 0 to nV-1.

  -DTIME: calculate the time. If interested only in the final size, this option should be turned off. */

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cstring>

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

int infrom_mat(int *nE, int *E, char *filename);
// read network data in which each row represents a link

int main(int argc, char **argv) {

  if (argc != 4) {
#ifdef MIX
    cerr << "Usage: sir-node-centric-LGA.out nV alpha p" << endl;
    cerr << "Well-mixed population assumed" << endl;
#else // MIX not defined
    cerr << "Usage: sir-node-centric-LGA.out infilename alpha p" << endl;
#endif // MIX
    cerr << "exponential distribution of inter-event times if alpha<0" << endl;
    cerr << "p: prob lambda does not change between consecutive inter-event times" << endl;
    exit(8);
  }

  init_genrand(time(NULL));
  int i,j,ii; // counters
  int vs; // source node of infection
  int ve; // destination node of infection
  int nV; // # nodes

#ifdef MIX // well-mixed population
  nV = atoi(argv[1]);
  if (nV <= 0 || nV > 1e+7) {
    cerr << "nV must be integer not too big" << endl;
    exit(8);
  }
#else // MIX not defined
  int nE; // # undirected links
  // preread nE
  char dummy[10];
  ifstream fin(argv[1]);
  if (!fin) {
    cerr << "cannot open infile.mat" << endl;
    exit(8);
  }
  fin >> dummy >> nV >> nE;
  fin.close();
  int *E = new int[2*nE]; // link list
  infrom_mat(&nE,E,argv[1]);
  fin.close();
  // nV (number of nodes), nE (number of links), and E[] (link list) determined
  
  int k[nV]; // degree
  for (i=0 ; i<nV ; i++) k[i] = 0; // initialization
  for (i=0 ; i<2*nE ; i++) k[E[i]]++;
  int kmax = 0; // max degree
  for (i=0 ; i<nV ; i++)
    if (k[i] > kmax) kmax = k[i];
  cerr << "kmax = " << kmax << endl;
  int adj_list[nV][kmax]; // list of adjacent nodes for each node
  for (i=0 ; i<nV ; i++) k[i] = 0; // to reuse k[i] as a counter
  for (i=0 ; i<nE ; i++) {
    vs = E[2*i];
    ve = E[2*i+1];
    adj_list[vs][k[vs]] = ve;
    adj_list[ve][k[ve]] = vs;
    k[vs]++;
    k[ve]++;
  }
#endif // MIX
  
  double alpha = atof(argv[2]); // shape parameter of the gamma distribution, psi(tau) = alpha/(1+tau)^{alpha+1}

  // preparation to construct a binary tree
  int log2_Np_upper = (int)(log(2*nV-1e-6)/log(2)) + 1;
  // 2^{log2_Np_upper-1} < nP <= 2^{log2_Np_upper}
  int pow2[log2_Np_upper+1]; // power[i] = 2^i.
  pow2[0]=1;
  for (i=0 ; i<log2_Np_upper ; i++) {
    pow2[i+1] = 2 * pow2[i];
  }

  // event rates and the binary tree
  double lambda[2*pow2[log2_Np_upper]];
  /* rate of Poisson processes and associated binary tree
     lambda[0], ..., lambda[nV-1]: recovery rate of nV nodes
       = mu for an infected node and =0 otherwise 
     lambda[nV], ..., lambda[2*nV-1]: infection rate of nV nodes
       As determined by p(lambda)
   */
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
    
  double p = atof(argv[3]); // probability to use the same exp dist
  if (p > 1 || p < 0) {
    cerr << "0 <= p <= 1 violated" << endl;
    exit(8);
  }

#ifdef TIME
  double t, dt; // time
#endif // TIME  

  double final_size; // final % R nodes 
  int st[nV]; // node's state. 0: S, 1: I, 2: R
  int nI,nR; // nI: # infected nodes, nR: # recovered nodes
  int if_infect; // working var
  double eff_inf_rate; // effective infection rate = <tau_rec> / <tau_inf>, as defined in Boguna et al. PRE 2014
  double eff_inf_rate_min = 0.1;
  double eff_inf_rate_max = 5;
  double mu; // recovery rate
  int Neff_inf_rate = 100; // 120; // # effective infection rate values
  int ind,tr; // counter
  int trials=1000;
  double ra; // random variate
  
  for (ind=0 ; ind<Neff_inf_rate ; ind++) {

    eff_inf_rate =(Neff_inf_rate>=2)?
	  eff_inf_rate_min + (eff_inf_rate_max-eff_inf_rate_min)*(double)ind/(Neff_inf_rate-1) : eff_inf_rate_min;
    mu = (alpha<0)? 1.0/eff_inf_rate : (alpha-1)/eff_inf_rate;
    /* adjust the recovery rate to realize the specified eff_inf_rate
       Note that alpha<0 corresponds to the Poisson case */
    
    final_size = 0.0; // initialization

    for (tr=0 ; tr<trials ; tr++) { // "trials" realizations of SIR dynamics in total, starting from the same initial condition

      // initialization
      nI = 1; // # I
      nR = 0; // # R
      st[0] = 1; // node 0 is initially infected
      lambda[0] = mu; // recovery rate
      for (i=1 ; i<nV ; i++) {
	st[i] = 0; // all the other nodes are initially susceptible
	lambda[i] = 0.0; // recovery rate
      }
      for (i=nV ; i<2*nV ; i++) // infection rate
	lambda[i] = (alpha<0)? 1.0 : rgamma(alpha,1);
        // alpha < 0 corresponds to Poisson processes
      for (i=2*nV ; i<2*pow2[log2_Np_upper] ; i++)
	lambda[i] = 0.0; // dummy
      for (i = 0 ; i < log2_Np_upper ; i++) { // construct binary tree
	for (j = 0 ; j < pow2[log2_Np_upper-i] ; j++)
	  lambda[start[i+1]+j/2] += lambda[start[i]+j];
      }

      // SIR dynamics
#ifdef TIME
      t = 0.0;
#endif // TIME      
      while (nI>0) { // there are infectious nodes

	// advance time and determine the process that has yielded the event
	i = 0;
#ifdef TIME
	dt = -log((genrand_int32()+0.5)/4294967296.0)/lambda[2*pow2[log2_Np_upper]-2];
#endif // TIME
	ra = (genrand_int32()+0.5)/4294967296.0 * lambda[2*pow2[log2_Np_upper]-2];
	// ra is uniformly distributed between 0 and lambda[2*pow2[log2_Np_upper]-2]
	for (j = log2_Np_upper-1 ; j >= 0 ; j--) { // fast search using binary tree
	  i *= 2;
	  if (ra > lambda[start[j]+i]) {
	    ra -= lambda[start[j]+i]; 
	    i++;
	  }
	}
	if (i >= 2*nV) {
	  cerr << "selection of the process failed" << endl;
	  exit(8);
	}
	// the node or link to fire (i.e., i) determined

#ifdef TIME
	t += dt;
#endif // TIME	

	if (i<nV) { // I -> R
	  st[i] = 2;
	  nI--;
	  nR++;

	  Dlambda = -mu;
	  lambda[i] = 0.0;
	  ii = i;
	  for (j = 0 ; j < log2_Np_upper ; j++) {
	    ii /= 2;
	    lambda[start[j+1]+ii] += Dlambda; // update the binary tree
	  }

	  Dlambda = -lambda[nV+i];
	  lambda[nV+i] = 0.0;
	  ii = nV+i;
	  for (j = 0 ; j < log2_Np_upper ; j++) {
	    ii /= 2;
	    lambda[start[j+1]+ii] += Dlambda; // update the binary tree
	  }
	} else { // Potentially S -> I
	  // update lambda[i], where nV <= i < 2*nV, corresponding to the (i-nV)th node, no matter whether infection occurs or not
	  if ((genrand_int32()+0.5)/4294967296.0 < 1-p) { // select a new rate lambda[i] with probability 1-p
	    Dlambda = lambda[i];
	    lambda[i] = (alpha<0)? 1.0 : rgamma(alpha,1);
	    Dlambda = lambda[i] - Dlambda; // increment in lambda[i]
	    ii = i;
	    for (j = 0 ; j < log2_Np_upper ; j++) {
	      ii /= 2;
	      lambda[start[j+1]+ii] += Dlambda; // update the binary tree
	    }
	  }
	  vs = i-nV; // source node of potential infection
#ifdef MIX // well-mixed population
	  ve = (vs + 1 + genrand_int32()%(nV-1)) % nV;
#else // MIX not defined
	  ve = adj_list[vs][genrand_int32()%k[vs]];
#endif // MIX
	  // infect ve if ve is not infected
	  if_infect = 0; // initialization
	  if (st[vs]==1 && st[ve]==0)
	    if_infect = 1; // infection will occur
	  else if (st[vs]==0 && st[ve]==1) { // swap vs and ve
	    j = vs;
	    vs = ve;
	    ve = j;
	    if_infect = 1;
	  }
	  if (if_infect==1) { // infection occurs
	    st[ve] = 1; // ve is now infected
	    nI++;
	    lambda[ve] = mu; // recovery rate of ve

	    Dlambda = mu;
	    ii = ve;
	    for (j = 0 ; j < log2_Np_upper ; j++) {
	      ii /= 2;
	      lambda[start[j+1]+ii] += Dlambda; // update the binary tree
	    }
	  } // if (if_infect==1) done
	} // a single potential state transition completed
      } // one realization of the SIR dynamics completed

      final_size += nR; // to calculate the prevalence

    } // all the realizations (= one value of mu) completed

    final_size /= nV*trials;
    cout << eff_inf_rate << " " << final_size << " " << mu << endl;
  } // all values of mu completed

  return 0;
}

// nE: # links, E[]: link list, filename: input file name
int infrom_mat(int *nE, int *E, char *filename) {
  int i;
  int nV; // # nodes
  ifstream fin(filename); // matlab .mat format
  if (!fin) {
      cerr << "infrom_mat: cannot open input file" << endl;
    exit(8);
  }
  double wdummy; // dummy var for the link weight (3rd column)
  char dummy[20];
  fin >> dummy >> nV >> *nE; // 1st row contains (# nodes) and (# links)
  for (i=0; i<(*nE) ; i++) {
    fin >> E[2*i] >> E[2*i+1] >> wdummy;
  }
  for (i=0 ; i<2*(*nE) ; i++) E[i]--;
  /* In the input file, the node index runs 1, ..., nV.
     It is transformed to 0, ..., nV-1. */
  fin.close();
  return nV;
}

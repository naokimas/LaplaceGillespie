/* Generate gamma-distributed random variables using the algorithm given by
Marsaglia, G. & Tsang, W. W. (2000).
A simple method for generating gamma variables.
ACM Transactions on Mathematical Software, 26, 363-372.
doi:10.1145/358407.358414 */

#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib> // atoi
#include <cstring>
#include <cmath>

#ifdef MAIN
#include "mt19937ar.c"
#endif // MAIN
/* Mersenne Twister to generate random variates on [0,1].
   mt19937ar.c is available at 
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c

   If one prefers to use the in-built random number generator,
   (genrand_int32()+0.5)/4294967296.0 should be replaced by (double)rand()/RAND_MAX
   and
   init_genrand(time(NULL)) should be replaced by, e.g., srand(time(NULL)) */

// Gaussian distribution with mean 0 and std sigma
double gauss(double sigma);

/* Marsaglia-Tsang algorithm
alpha: shape parameter
kappa: scale parameter */
double rgamma(double alpha, double kappa) {
    double d, c, x, v, ra1;
    if(alpha >= 1.0){
        d = alpha - 1.0/3.0; 
        c = 1.0/sqrt(9.0*d);
        while(1) {
            do {
                x=gauss(1.0); 
                v=1.0 + c*x;
            } while(v<1e-12);
            v = v * v * v; 
            ra1 = (genrand_int32()+0.5)/4294967296.0; // uniform on [0, 1]
            if (ra1 < 1.0 - 0.0331*(x*x)*(x*x)) {
                return d*v*kappa;
            }
            if (log(ra1) < 0.5*x*x + d*(1.0-v+log(v))) {
                return d*v*kappa;
            }
        }
    } else {
        x = rgamma(alpha+1, kappa);
        x = x * pow((genrand_int32()+0.5)/4294967296.0, 1.0/alpha); 
        return x;
    }
}

// Gaussian distribution with mean 0 and std sigma
double gauss(double sigma) {
    double x1, x2, z;
    do {
        x1 = 2.0 * (genrand_int32()+0.5)/4294967296.0 - 1.0; // uniform on [-1, 1]
        x2 = 2.0 * (genrand_int32()+0.5)/4294967296.0 - 1.0;
        z = x1 * x1 + x2 * x2;
    } while (z >= 1.0);

    return sigma * x1 * sqrt((-2.0 * log(z)) / z);
    // x1/sqrt(z) = cos(theta), where theta is uniformly distributed on [0, 2*pi]
}

#ifdef MAIN
// generate samples and compare the mean and std between the samples and theory
int main() {
    init_genrand(time(NULL));
    int i; // counter
    int n = 10000; // number of samples
    double x[n]; // generated samples
    double alpha = 4.0; // shape parameter
    double kappa = 1.0; // scale parameter
    double ave = 0.0; // sample average
    double std = 0.0; // sample standard deviation

    for(i=0; i<n; i++){
        x[i] = rgamma(alpha, kappa);
        cout << x[i] << endl;
        ave += x[i];
        std += x[i] * x[i];
    }
    ave /= n;
    std = sqrt(std/n - ave*ave);
    cerr << "Empirical: ave = " << ave << ", std = " << std << endl;
    cerr << "Theoretical: ave = " << alpha*kappa << ", std = " << sqrt(alpha*kappa*kappa) << endl; // theoretical mean and std of the gamma dist
    return 0;
}
#endif // MAIN
/*

HE-KDE estimator, wrapper for C++ library HEAAN-1.0.
To be called by julia code after compilation.

8-core processor is assumed; adjust `numThread` in `main` according to actual cores.

Compile by:
```bash
g++ -I./HEAAN-1.0/HEAAN/src/ HEKDE.cpp -o HEKDE.out -L./HEAAN-1.0/HEAAN/lib -lHEAAN -L.ntl-11.4.3/src -lntl -L/usr/lib/x86_64-linux-gnu/ -lgmp -lpthread
```
run by
```bash
./HEKDE.out 11 128 0.1 14 3 3000000000000000000000.2
```
where the example numbers used are
- 11 = length of query point vector
- 128 = number of observations
- 0.1 = bandwidth h
- 14 = degree d
- 3 = shifth tau
- 3000000000000000000000.2 = normalising constant c_{d, tau}
(Note: this value is not the actual normalising constant for d=14, tau=3)

*/


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "HEAAN.h"

using namespace NTL;
using namespace std;

double* readX(int n_x){
  double *x = new double[n_x];
  ifstream file_in;
  file_in.open("./temp/x_data.txt");
  if (file_in.is_open()){
    double element;
    int i = 0;
    while (i < n_x){
      file_in >> element;
      x[i] = element;
      // cout << "index " << i << " element: " << element << endl;
      ++i;
    }
  }
  else throw "file does not exist";
  return x;
}

double* readXQuery(int n_query){
  double *x_query = new double[n_query];
  ifstream file_in;
  file_in.open("./temp/x_query_data.txt");
  if (file_in.is_open()){
    double element;
    int i = 0;
    while (i < n_query){
      file_in >> element;
      x_query[i] = element;
      // cout << "index " << i << " element: " << element << endl;
      ++i;
    }
  }
  else cerr << "file does not exist" << endl;
  return x_query;
}


void writeResult(double *result, int n_query){
  ofstream file_out;
  file_out.open("./temp/result.txt");
  if (file_out.is_open()){
    for(int i = 0; i<n_query; i++){
      // cout << "result[" << i << "] = " << result[i] << endl;
      file_out << result[i] << "\n";
    }
    file_out.close();
  }
}


double myRiemannZeta(double x){
  double zeta = 0.0;
  int reps = 1000;
  for (int rep=1; rep < reps; ++rep){
    zeta += 1 / pow(rep, x);
  }
  return zeta;
}


double kernelWeight(int k){
  double w;
  w = pow(-1, k + 1) / k + pow(-1, k) * myRiemannZeta(k + 1.0);
  return w;
}

/*
Compute Ktilde(u)/(n*h) where Ktilde(u) is the normalised HE-kernel.
*/
Ciphertext normalisedKernel(Ciphertext c_u, Scheme &scheme,
                 long logn, long logp,
                 double h, int d, int tau, double alpha){
   // numbers
   const double gamma =  0.5772156649015328606065120900824024310421;
   const double pi =     3.1415926535897932384626433832795028841971;
   long n = (1 << logn);
   long n_x = (1 << (logn - 1));
   double norm_factor = 1 / (alpha * n_x * h * pi);

   // allocate ciphertexts
   Ciphertext c_kernel; // 1-tau+u^2
   Ciphertext c_power;  // (1-tau+u^2)^j for j=1,...,k
   Ciphertext c_result; // sum of terms
   Ciphertext c_term;   // actual term to add to result = w_j*(1-tau+u^2)^j for j=1,...,k

   // prepare terms for powers
   c_kernel = scheme.square(c_u);
   scheme.reScaleByAndEqual(c_kernel, logp); // u^2
   scheme.addConstAndEqual(c_kernel, 1 - tau ,logp); // 1-tau+u^2
   // ciphertext to compute powers iteratively
   c_power = c_kernel;
   // first term
   double w_1 = kernelWeight(1) * norm_factor;
   c_result = scheme.multByConst(c_power, w_1, logp); // w_1*(1-tau+u^2)/(alpha*n_x*h*pi)
   scheme.reScaleByAndEqual(c_result, logp);
   // other terms
   for (int k = 2; k <= d; ++k){
     // power
     c_power = scheme.mult(c_power, c_kernel);  // (1-tau+u^2)^k
     scheme.reScaleByAndEqual(c_power, logp);
     // coefficient
     double w_k = kernelWeight(k) * norm_factor;
     c_term = scheme.multByConst(c_power, w_k, logp); // normal_factor*(w_k*(1-tau+u^2)^k)
     scheme.reScaleByAndEqual(c_term, logp);
     // add to result
     scheme.modDownToAndEqual(c_result, c_term.logq);  // sum_{j=1}^{k} w_j * (1-tau+u^2)^j
     c_result = scheme.add(c_result, c_term);
   }
   // add normalised gamma
   scheme.addConstAndEqual(c_result, gamma * norm_factor, logp);  // normal_factor * (gamma + sum_{j=1}^{k} w_j * (1-tau+u^2)^j)

   return c_result;
}

/*
Compute the encryption `ftilde_{d,tau}(x0)` where
- `c_query` is the encryption of the scalar query point x0 (broadcasted to an
            n_x-long vector)
- `c_x` is the encryption of the iid sample from the estimand density
- `alpha` is the normalising factor as defined in the julia code
*/
Ciphertext hekde(Ciphertext c_query, Ciphertext c_x, Scheme &scheme,
                 long logn, long logp,
                 double h, int d, int tau, double alpha){
  // allocate ciphertexts
  Ciphertext c_u;
  Ciphertext c_result;

  // input to kernel
  c_u = scheme.sub(c_query, c_x);  // x-x_i
  scheme.multByConstAndEqual(c_u, 1/h, logp);  // (x-x_i)/h =: u
  scheme.reScaleByAndEqual(c_u, logp);
  // normalised kernel
  c_result = normalisedKernel(c_u,scheme, logn,logp, h,d,tau,alpha);
  return c_result;
}



int main(int argc, char *argv[]){
  // parse arguments
  if (argc<7){
    cerr << "not enough arguments" << endl;
    return 1;
  }
  int n_query = atoi(argv[1]);
  long n = atol(argv[2]);
  double h = atof(argv[3]);
  int d = atoi(argv[4]);
  int tau = atoi(argv[5]);
  double alpha = atof(argv[6]);

  // constants
  const double pi = 3.141592653589793;

  // he parameters
  long logn = log2(n);  // N
  long logq = 1380;  // q_L
  long logp = 80;  // log2(Delta)
  long n_x = n / 2;  // N/2 length of vector to encrypt
  long slots = n_x;
  long numThread = 8;  // number of threads for parallel computing
  int zoH = 64; // parameter of ternary distribution HWT(zoH)

  // reading data
  double *x_query = readXQuery(n_query);
  double *x = readX(n_x);

  // ------------------------------- HE-KDE:
  // setup and keygen
  srand(time(NULL));
  SetNumThreads(numThread);
  SecretKey secretKey(logn, 64);
  Context context(logn, logq);
  Scheme scheme(secretKey, context);
  // encrypt sample data
  Ciphertext c_x = scheme.encrypt(x,slots,logp,logq);
  // he-kde for each point in query points
  double *result = new double[n_query];
  for (int i = 0; i < n_query; ++i){
    // broadcast x_query[i] to an n_x-long array
    double *x_query_broadcast = new double[n_x];
    for (int j = 0; j < n_x; ++j){
      x_query_broadcast[j] = x_query[i];
    }
    // encrypt broadcasted array
    Ciphertext c_query = scheme.encrypt(x_query_broadcast,slots,logp,logq);
    // compute he-kde
    Ciphertext c_result = hekde(c_query,c_x, scheme,logn,logp, h,d,tau,alpha);
    // decrypt result
    complex<double>* hekde_out = scheme.decrypt(secretKey, c_result);
    // summation
    result[i] = 0.0; //real(hekde_out[0]);
    for (int j = 0; j < n_x; ++j){
      // cout << "hekde out index " << j << ": " << hekde_out[j] << endl;
      result[i] += real(hekde_out[j]);
      // cout << "result[" << i << "] =" << result[i] << endl;
    }
  }

  // writing result
  writeResult(result, n_query);

  return 0;
}

/*
double *x_query = new double[n_query];
ifstream file_in;
file_in.open("./temp/x_querydata.txt");
if (file_in.is_open()){
  double element;
  int i = 0;
  while (file_in >> element){
    x_query[i++] = element;
  }
}
*/

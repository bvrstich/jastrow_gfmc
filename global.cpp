#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

namespace global{

   int Lx;
   int Ly;

   int omp_num_threads;

   double f;

   Random RN;

   /**
    * @param Lx_in x dimension of the square lattice
    * @param Ly_in y dimension of the square lattice
    * @param f_in jastrow wavefunction parameter
    */
   void init(int Lx_in,int Ly_in,double f_in){

      Lx = Lx_in;
      Ly = Ly_in;

      f = f_in;

#ifdef _OPENMP
      omp_num_threads = omp_get_max_threads();
#else
      omp_num_threads = 1;
#endif

   }

   //!function which generates random complex numbers uniformly on a square of side 2 [(-1,1):(-1,1)]
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

   //!function which generates uniform random numbers between [0:1]
   template<>
      double rgen_pos(){ 

         return RN();

      }

}

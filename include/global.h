#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

#include "Random.h"

class Walker;

namespace global {

   extern Random RN;

   //!x dimension of the lattice, nr of cols
   extern int Lx;

   //!y dimension of the lattice, nr of rows
   extern int Ly;

   //! parameter of the jastrow wavefunction
   extern double f;

   //!number of omp threads
   extern int omp_num_threads;

   void init(int,int,double);

   template<typename T>
      T rgen();

   template<typename T>
      T rgen_pos();
   
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;
using std::vector;

#include "include.h"

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);
   int Nw = atoi(argv[2]);
   double f = atof(argv[3]);

   double dtau = 0.01;

   //initialize the dimensions of the problem, set the trial
   global::init(L,L,f,dtau);

   GFMC gfmc(Nw);
   gfmc.walk(1000000);

}

#ifndef GFMC_H
#define GFMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

class Walker;
class Distribution;

class GFMC {

   public:
   
      //constructor with input trialwavefunction
      GFMC(int);
      
      //Destructor
      virtual ~GFMC();
      
      //Let the walkers propagate for n_steps steps
      void walk(int);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double propagate();

      void PopulationControl(double);
      
      //Write the projected energy, target energy
      void write(int);

      //Setup the walkers
      void SetupWalkers();

      void dump(const char *);

   private:
      
      //!The total desired number of walkers
      int Nw;

      //!projected energy at current timestep
      double EP;

      //!The walkers
      std::vector<Walker> walker;

      //!Distribution of possible final states given a walker state
      std::vector< Distribution > dist;
      
};

#endif

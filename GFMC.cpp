#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;
using std::complex;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

using namespace global;

/**
 * constructor of the GFMC object, takes input parameters that define the QMC walk.
 * @param Nw_in number of Walker states
 */
GFMC::GFMC(int Nw_in){

   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
GFMC::~GFMC(){ }

/**
 * initialize the walkers
 */
void GFMC::SetupWalkers(){

   walker.resize(Nw);

   char filename[200];
   sprintf(filename,"/home/bright/bestanden/results/jastrow/vmc/%dx%d/data/f=%f.walk",Lx,Ly,f);

   ifstream in(filename);

   bool tmp;

   for(int i = 0;i < Nw;++i)
      for(unsigned int j = 0;j < walker[i].size();++j){

         in >> tmp;

         walker[i][j] = tmp;

      }

}

void GFMC::walk(const int n_steps){

   char filename[200];
   sprintf(filename,"output/%dx%d/f=%f.txt",Lx,Ly,global::f);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t" << endl;
   output.close();

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      double wsum = propagate();

      double scaling = Nw / wsum;

      ET = 1.0/dtau * ( 1 - 1.0/scaling);

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP / (double)(Lx*Ly) << endl;
      cout << "         E_T = " << ET / (double)(Lx*Ly) << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

      PopulationControl(scaling);

      double max_ov = 0.0;
      double min_ov = 1.0;

      for(unsigned int i = 0;i < walker.size();++i){

         if(max_ov < std::abs(walker[i].gOverlap()))
            max_ov = std::abs(walker[i].gOverlap());

         if(min_ov > std::abs(walker[i].gOverlap()))
            min_ov = std::abs(walker[i].gOverlap());

      }

#ifdef _DEBUG
      cout << "Minimal Overlap:\t" << min_ov << endl;
      cout << "Maximal Overlap:\t" << max_ov << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double GFMC::propagate(){

   double total_weight = 0.0;
   double sum = 0.0;

#pragma omp parallel for reduction(+:sum,total_weight)
   for(unsigned int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(walker[i],dtau,0.0);
      double nrm = dist[myID].normalize();

      //draw new walker
      int pick = dist[myID].draw();

      //copy the correct walker
      walker[i] = dist[myID].gwalker(pick);

      //set the energy
      sum += walker[i].gWeight() * dist[myID].energy();
      
      //multiply weight
      walker[i].multWeight(nrm);

      total_weight += walker[i].gWeight();

   }

   EP = sum / (double) Nw;

   return total_weight;

}

/**
 * write output to file
 */
void GFMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/%dx%d/f=%f.txt",Lx,Ly,global::f);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP/(double)(Lx*Ly) << "\t" << ET /(double)(Lx*Ly) << endl;
   output.close();

}

/**
 * dump the walkers to a single file
 */
void GFMC::dump(const char *filename){

   ofstream out(filename);
   out.precision(16);

   for(unsigned int i = 0;i < walker.size();++i){

      for(int row = 0;row < Ly;++row)
         for(int col = 0;col < Lx;++col)
            out << walker[i][row*Lx + col] << " ";
      out << endl;

   }

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void GFMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   for(unsigned int i = 0;i < walker.size();i++)
      walker[i].multWeight(scaling);

   for(unsigned int i = 0;i < walker.size();i++){

      double weight = walker[i].gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.1){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + RN());

         if(nCopies == 0){

#ifdef _DEBUG
            cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

            walker.erase(walker.begin() + i);

         }
         else
            walker[i].sWeight(1.0);

      }

      if(weight > 2.0){ //statically energy doesn't change

         int nCopies =(int) ( weight + RN());
         double new_weight = weight / (double) nCopies;

         walker[i].sWeight(new_weight);

#ifdef _DEBUG
         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

         for(int n = 1;n < nCopies;++n){

            Walker nw = walker[i];
            walker.push_back(nw);

         }

      }

   }

   double sum = 0.0;

   for(unsigned int i = 0;i < walker.size();++i)
      sum += walker[i].gWeight();
   
   //rescale the weights to unity for correct ET estimate in next iteration
   for(unsigned int i = 0;i < walker.size();++i)
      walker[i].multWeight((double)Nw/sum);


#ifdef _DEBUG
   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight is " << minw << " ." << endl;
   cout << "The max. encountered weight is " << maxw << " ." << endl;
   cout << endl;
#endif

}

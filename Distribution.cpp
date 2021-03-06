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

using namespace global;

/** 
 * empty constructor
 */
Distribution::Distribution() : vector<double> () { }

/** 
 * copy constructor
 */
Distribution::Distribution(const Distribution &dist_copy) : vector<double>(dist_copy){

   list = dist_copy.glist();

}

/** 
 * empty destructor
 */
Distribution::~Distribution(){ }

/** 
 * @return number of grid points
 */
int Distribution::gn() const {

   return this->size();

}

ostream &operator<<(ostream &output,const Distribution &dist_p){

   for(int i = 0;i < dist_p.gn();++i)
      output << i << "\t" << dist_p[i] << endl;

   return output;

}

/** 
 * @return the list of final Walker states
 */
const vector< Walker > &Distribution::glist() const {

   return list;

}

/** 
 * @return the final state Walker with index i
 * const version
 */
const Walker &Distribution::gwalker(int index) const {

   return list[index];

}

/** 
 * @return the final state Walker with index i
 */
Walker &Distribution::gwalker(int index) {

   return list[index];

}

/**
 * construct and fill the distribution by calculating the matrix elements <0|1-dtau * H|i> for all i = 0,...,n
 * @param walker_i input walker, the list and distribution is constructed from this
 * @param tau time step value, needs to be small enough
 * @param ET target energy, used to rescale the population, should be 'exact' energy
 */
void Distribution::construct(const Walker &walker_i,double dtau,double ET){

   //first reset the lists
   list.clear();
   this->clear();

   //f == i first 
   list.push_back(walker_i);

   //first horizontal 'final' states
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx;++c){

         if(walker_i[r*Lx + c] != walker_i[r*Lx + (c + 1)%Lx]){

            Walker walker_f(walker_i);

            walker_f[r*Lx + c] = !(walker_i[r*Lx + c]);
            walker_f[r*Lx + (c + 1)%Lx] = !(walker_i[r*Lx + (c + 1)%Lx]);

            list.push_back(walker_f);

         }

      }

   }

   //then vertical 'final' states
   for(int c = 0;c < Lx;++c){

      for(int r = 0;r < Ly;++r){

         if(walker_i[r*Lx + c] != walker_i[( (r + 1)%Ly )*Lx + c]){

            Walker walker_f(walker_i);

            walker_f[r*Lx + c] = !(walker_i[r*Lx + c]);
            walker_f[( (r + 1)%Ly )*Lx + c] = !(walker_i[( (r + 1)%Ly )*Lx + c]);

            list.push_back(walker_f);

         }

      }

   }

   for(unsigned int i = 0;i < list.size();++i)
      list[i].calc_overlap();

   this->resize(list.size());

   (*this)[0] = 1.0 - dtau * (walker_i.pot_en() - ET);

   for(unsigned int i = 1;i < list.size();++i)
      (*this)[i] = 0.5 * dtau * ( list[i].gOverlap() / list[0].gOverlap() );

}

/** 
 * get the local energy for the 'seed' walker, first the distribution has to be filled!
 */
double Distribution::energy() const {

   double energy = list[0].pot_en();

   for(unsigned int i = 1;i < list.size();++i)
      energy -= 0.5 * list[i].gOverlap() / list[0].gOverlap();

   return energy;

}

/**
 * draw a number from the distribution
 */
int Distribution::draw() const {

   //Get what you should do
   int trial = RN()*this->size();
   double x = RN();

   while((*this)[trial] < x){

      trial = RN()*this->size();
      x = RN();

   }

   return trial;

}

/** 
 * normalize the distribution
 * @return the total weight/norm of the distrubition
 */
double Distribution::normalize(){

   double nrm = 0.0;

   for(unsigned int i = 0;i < this->size();++i)
      nrm += (*this)[i];

   for(unsigned int i = 0;i < this->size();++i)
      (*this)[i] /= nrm;

   return nrm;

}

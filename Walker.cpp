#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

#include "include.h"

using namespace global;

/**
 * construct a Walker object: initialize on AF state
 * @param n_trot_in number of trotter terms
 * @param option start with spin up or spin down
 */
Walker::Walker(int option) : std::vector< bool >( Lx * Ly ){

   weight = 1.0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         if( (r + c + option)%2 == 0)
            (*this)[ r*Lx + c ] = true;
         else
            (*this)[ r*Lx + c ] = false;

      }

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< bool >(walker) {

   this->weight = walker.gWeight();
   this->overlap = walker.gOverlap();

}

/**
 * destructor
 */
Walker::~Walker(){ }

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}

/** 
 * @return the overlap of the walker with the Trial
 */
double Walker::gOverlap() const{

   return overlap; 

}

ostream &operator<<(ostream &output,const Walker &walker_p){

   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx;++c)
         output << walker_p[r*Lx + c] << " ";

      output << endl;

   }

   return output;

}

/**
 * get the 'potential' energy of the walker: < \sum_i Sz_i Sz_{i+1} >
 */
double Walker::pot_en() const {

   double tmp = 0.0;

   //first horizontal
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx;++c){

         //Sz Sz
         if( (*this)[r*Lx + c] == (*this)[r*Lx + (c + 1)%Lx] )//up up or down down
            tmp += 0.25;
         else //up down or down up
            tmp -= 0.25;

      }

   }

   //then vertical
   for(int c = 0;c < Lx;++c){

      for(int r = 0;r < Ly;++r){

         //Sz Sz
         if( (*this)[r*Lx + c] == (*this)[( (r + 1)%Ly ) *Lx + c] )//up up or down down
            tmp += 0.25;
         else //up down or down up
            tmp -= 0.25;

      }

   }

   return tmp;

}

/**
 * save the walker to file
 * @param filename name of the file
 */
void Walker::save(const char *filename){

   ofstream fout(filename);

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col)
         fout << row << "\t" << col << "\t" << (*this)[row*Lx +col] << endl;

}

/**
 * load the walker
 */
void Walker::load(const char *filename){

   ifstream fin(filename);

   bool tmp;

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         fin >> row >> col >> tmp; 
         (*this)[row*Lx + col] = tmp;

      }

}

/**
 * calculate the overlap of the Walke with the jastrow trial
 */
void Walker::calc_overlap(){

   overlap = 1.0;

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         int ind = row*Lx + col;

         if( (*this)[ind] == true ){

            //neighbours
            int nr = ( (row + 1)%Ly ) * Lx + col;
            int nu = row * Lx + (col + 1)%Lx;

            if((*this)[ind] == (*this)[nr])
               overlap *= global::f;

            if((*this)[ind] == (*this)[nu])
               overlap *= global::f;

         }

      }

}

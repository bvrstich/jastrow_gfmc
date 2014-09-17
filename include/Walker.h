#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

/**
 * class definition of Walker, made to describe the product state walkers. An array of L*L booleans vector. Each site represents a spin up (0) or spin down (1)
 */
class Walker : public vector< bool > {

   friend std::ostream &operator<<(std::ostream &output,const Walker &walker_p);

   public:

      //empty contstructor
      Walker(int option = 0);
   
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();

      double gWeight() const;

      void sWeight(double);

      void multWeight(double);

      void calc_overlap();

      double gOverlap() const;

      double gEL() const;

      double pot_en() const;

      void save(const char *);

      void load(const char *);

  private:

      //!The walker weight
      double weight;

      //!local energy
      double EL;

      //! overlap
      double overlap;


};

#endif

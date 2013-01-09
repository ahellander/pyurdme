

#include <stdlib.h>
#include "propensities.h"
#include "mex.h"
#include <math.h>

#define MEM                     1
#define PATCH                   2
#define PATCH2                  2

#define NR 1

#define A 0

#define MALLOC malloc
#define FREE free

const double ka = 1.0e100;
//const double ka_micro = 

const double rpatch = 0.05;

double rFun1(const int *x, double t, const double vol, const double *data, int sd)
/* A+B -> EmptySet*/
{
   /* Distande from south-pole */
  // double dist = sqrt(data[0]*data[0]+data[1]*data[1]);
   //return ka*x[A]*x[B];
    
   //if (dist <= rpatch && data[2]>0)
   if (sd==PATCH||sd==PATCH2)
    return ka*x[A];
   else 
     return 0.0;
}


PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)MALLOC(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
    
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  FREE(ptr);
}

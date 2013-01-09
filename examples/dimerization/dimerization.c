/* Model file for simple dimerization example. */ 

#include <stdlib.h>
#include <math.h>
#include "propensities.h"
#include <stdio.h>

#define A 0
#define B 1
#define C 2

#define NR		 2
#define NS       3

# define pi 3.141592653589793
/* Rate constants. */

//#define ELF
//#define HELL
const double NA   = 6.0221415e23;

/* parameters[0] is the microscopic constant. */
const double rho = 2e-9;
const double D   = 2e-14;
const double L   = 500e-9; 

/*inline double F(x)
{
  return log(1/x)/(1-pow(pow(x,2),2));  
}*/

inline double F(double x)
{
  return (4*log(1.0/x)-(1.0-x*x)*(3-x*x))/(4*(1-x*x)*(1-x*x));   
}
inline double mesoreac2D(double rho,double vol,double gamma,double kr)
{
   double R;
   R = sqrt(vol/pi);
   double lambda;
   lambda = rho/R;
   double alpha;
   alpha = kr/(2.0*pi*gamma);
   
   return pi*R*R/kr*(1.0+alpha*F(lambda));
    
}


/* Reaction propensities. */
double rFun1(const int *x, double t, const double vol, const double *data, int sd)
/* A + B -> C */
{
    
    #if defined(ELF)
      /* Discretization independent propensities of Fange et. al. (2D) */
      double h = sqrt(vol/pi)-rho;
      double alpha =parameters[0]/(2.0*pi*D);
      double beta = rho/(rho+h);
      double ka_meso = parameters[0]/(1.0+alpha*log(1.0+0.544*(1.0-beta)/beta));
      return ka_meso*x[A]*x[B]/vol; 
    #elif defined(HELL)
      /* Hellander correction */ 
      double h = sqrt(vol);
      double lsq;
      lsq= pow(L/h,2);
      double mhs;
      mhs = h*h/(4.0*D);
      double ka_meso = (1.0+lsq)/(mesoreac2D(rho,L*L,D,parameters[0])-mhs*(lsq*log(lsq)/pi+0.1951*lsq));
      if (ka_meso<0.0) // h < h^*
        ka_meso=1e80;
      return ka_meso*x[A]*x[B];  
    #else
      return parameters[0]*x[A]*x[B]/vol;
    #endif
}

double rFun2(const int *x, double t, const double vol, const double *data, int sd)
/*  C -> A+B */
{   
    
    double ka_meso;
    
    #if defined(ELF)
      /* Discretization independent propensities of Fange et. al. (2D) */
      double h = sqrt(vol/pi)-rho;
      double alpha =parameters[0]/(2.0*pi*D);
      double beta = rho/(rho+h);
      ka_meso = parameters[0]/(1.0+alpha*log(1.0+0.544*(1.0-beta)/beta));
    #elif defined(HELL)
      /* Hellander correction */ 
      double h = sqrt(vol);
      double lsq;
      lsq= pow(L/h,2);
      double mhs;
      mhs = h*h/(4.0*D);
      ka_meso = (1.0+lsq)/(mesoreac2D(rho,L*L,D,parameters[0])-mhs*(lsq*log(lsq)/pi+0.1951*lsq));
      if (ka_meso<0.0)
         ka_meso=1e80;
      ka_meso=ka_meso*vol;
    #else
      ka_meso=parameters[0];
    #endif
    
    double kd_meso = parameters[1]*ka_meso/parameters[0];
    return kd_meso*x[C];      
}


PropensityFun *ALLOC_propensities(void)
{
  PropensityFun *ptr = (PropensityFun *)malloc(sizeof(PropensityFun)*NR);
  
  ptr[0]=rFun1;
  ptr[1]=rFun2;
 	
  return ptr;
}

void FREE_propensities(PropensityFun* ptr)
{
  free(ptr);
}


#include "reactions.h"
#include "species.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>


double transform_rate_constant(reaction *r,species **spec, double vol)
{
	
	double k;
	species *s1,*s2;
	s1 = spec[r->reactants[0]];
	s2 = spec[r->reactants[1]];
	double D = s1->gamma+s2->gamma;
	double sigma = s1->sigma+s2->sigma;
	
	/* Conventional mesoscopic reaction constant in 3D. */
	//k = 4.0*pi*sigma*D*r->k/(4.0*pi*sigma*D+r->k);
    /* Ad hoc consant in 2D */
	//k = 4.0*pi*D;
	k = r->k;
	return k;
    
}

/* Correction to the activity coefficients due to finite sized molecules. */

double epsilon(int l, int *x, int sv, double vol, species **spec,int Mspecies)
{
    double pi = 3.14159265358979323846;
    /* Dimension, 2D or 3D */
    double D = 2.0;

    double epsilon = 0.0;
    
    int j;
    for (j=0;j<Mspecies;j++){
        epsilon += (x[sv*Mspecies+j]/vol)*pow(2.0*spec[j]->sigma,l);
    }
    epsilon *= (pi/(2.0*D));
    return epsilon;
}

double activity_hard_disk(int *x, int sv, double vol, species **spec,int Mspecies,int s){
/* Compute the activity coefficient gamma for a hard disk (2D)
 
   Asserts: gamma >= 0.0

 */
    
    double gamma = 0.0;
    double e0 = epsilon(0,x,sv,vol,spec,Mspecies);
    double e1 = epsilon(1,x,sv,vol,spec,Mspecies);
    double e2 = epsilon(2,x,sv,vol,spec,Mspecies);
        
    double a = 1.0-e2;
    
    gamma -= log(a);
    gamma += 4.0*e1/a*spec[s]->sigma;
    gamma += ((4.0*e0/a + 4.0*e1*e1/(a*a))*(spec[s]->sigma*spec[s]->sigma));
    
   // assert(exp(gamma) >= 1.0);
    if (exp(gamma) < 1.0){
        printf("gamma: %f\n",exp(gamma));
        exit(-1);
    }
        
    

    return exp(gamma);

}

/* Evaluate mass action elementary function for reaction r with reactant 1 from sv 1 and reactant 2 from sv2.
 
    Asserts: 
        Scaling factor (non-ideal contribution) due to crowding, corr_factor >= 0.0
*/

double evaluate_elementary_propensity(reaction *r,int *x,species **spec, int sv1, double vol1, int sv2, double vol2,int Mspecies)
{
	
	double a=0.0;
	int order = r->order;  
	double k = r->k;
    double corr_factor = 1.0;

  	if (order == 0){
		a = k;
	}
	else if (order == 1){
#ifdef HARD_SPHERES
        corr_factor = activity_hard_disk(x,sv1,vol1,spec,Mspecies,r->reactants[0])/activity_hard_disk(x,sv1,vol1,spec,Mspecies,r->products[0]);
#endif
      //  assert(corr_factor>=1.0);
#ifdef __DEBUG
        if (corr_factor <0.0){
            printf("corr_factor: %f\n",corr_factor);
            exit(-1);
            
        }
#endif
	    a = k*corr_factor*x[sv1*Mspecies+r->reactants[0]];
	}
	else{
        /* Compute the available volume */
#ifdef HARD_SPHERES
        corr_factor = activity_hard_disk(x,sv1,vol1,spec,Mspecies,r->reactants[0])*activity_hard_disk(x,sv1,vol1,spec,Mspecies,r->reactants[1])/activity_hard_disk(x,sv1,vol1,spec,Mspecies,r->products[0]);
#endif

#ifdef __DEBUG
        if (corr_factor <0.0){
            printf("corr_factor: %f\n",corr_factor);
            exit(-1);
 
        }
#endif
		a = k*corr_factor*x[sv1*Mspecies+r->reactants[0]]*x[sv2*Mspecies+r->reactants[1]]/vol1;
		
	}
    return a;
	
}

/* 
 
 Evaluate propensity function.
    
 Asserts:
    Propensity is positive, a >= 0.0

*/
double evaluate_propensity(reaction *r,int *x,species **spec, int sv1, double vol1, int sv2, double vol2,int Mspecies)
{
    double a = evaluate_elementary_propensity(r,x,spec,sv1,vol1,sv2,vol2,Mspecies);
#ifdef __DEBUG
    if (a < 0.0){
        printf("a: %f\n",a);
        exit(-1);

    }
#endif
    
    return a;
    
}

/* Execute a reaction in subvolume sv. */
void execute_reaction(reaction *r, int *x, int sv, int Mspecies)
{
	execute_elementary_reaction(r,x,sv,Mspecies);

}

/* 
 
 Execute mass action reaction in subvolume sv.
 
 Asserts:
    Reaction, r, is not NULL.
    No copy number becomes negative.
 
*/
void execute_elementary_reaction(reaction *r, int *x, int sv, int Mspecies)
{
	
//assert(r != NULL);
    
	int s;
	for (s=0; s<Mspecies; s++) {
		x[sv*Mspecies+s]+=r->nr[s];
        if (x[sv*Mspecies+s]<0) {
            print_reaction(r);
            exit(-1);
        }
	}
}

void print_reaction(reaction *r)
{
	
  int i;	
  printf("\n");
  printf("Order: %i\n",r->order);
  printf("Reactants: ");
  for (i=0; i<r->nr_reactants; i++) {
	  printf("%i ",r->reactants[i]);
  }	
  printf("\n");	
  printf("Products: ");
  for (i=0; i<r->nr_products; i++) {
	 printf("%i ",r->products[i]);
  }	
  printf("\n");		
}

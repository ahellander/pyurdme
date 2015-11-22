/* reactions.h */

/* A. Hellander 2012-01-27. */  

#ifndef REACTIONS__H
#define REACTIONS__H

#include "species.h"

/* Definition of a chemical reaction. */
typedef struct{
	
	/* Stoichiometry of reaction. 0 = zeroth order, 1 = first order, 2 = bimolecular. 
	   Only relevant for elementary mass action kinetics. */
    
	int order;
	
    /* List of reactants.*/
	int nr_reactants;
	int *reactants;
	
    /* List of products. */
	int nr_products;
	int *products;
	
	/* Stoichiometry of the reaction.*/
	int *nr;
	
    /* rate constant. */
	double k;
	
	/* Other reactions affected by this reaction (for optimization) 
	   Not supported yet, will replace the dependency graph. */
	int *gr;
	
} reaction;

reaction **initialize_reactions(void);


int number_of_reactions(void);


double evaluate_propensity(reaction *r,int *x, species **list_of_spec,int sv1, double vol1, int sv2, double vol2,int Mspecies);
double evaluate_elementary_propensity(reaction *r,int *x, species **list_of_spec, int sv1, double vol1, int sv2, double vol2,int Mspecies);

void execute_reaction(reaction *r, int *x, int sv, int Mspecies);
void execute_elementary_reaction(reaction *r, int *x, int sv, int Mspecies);


void print_reaction(reaction *r);

double transform_rate_constant(reaction *r,species **spec, double vol);
double epsilon(int l, int *x, int sv, double vol, species **spec,int Mspecies);
double activity_hard_disk(int *x, int sv, double vol, species **spec,int Mspecies,int s);


double *parameters;

#endif /* REACTIONS__H */

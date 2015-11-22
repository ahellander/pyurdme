/* A. Hellander 2012-01-29 */

#ifndef SPECIES__H
#define SPECIES__H

/* Defines a chemical species. */
typedef struct{

	/* Can be used to name the species */
	char *name;
	
	/* Diffusion constant. Should maybe be a double * of size Ncells to
	   allow for spatially varying rate constants. */
	double gamma;
	//double *gamma;
	
	/* Reaction radius. */
	double sigma;
	
	
} species;

species **initialize_species(void);
int number_of_species(void);

#endif
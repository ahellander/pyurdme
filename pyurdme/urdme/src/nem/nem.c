/* nlnsm solver. */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "reactions.h"
#include "species.h"
#include "nem.h"
#include "urdmemodel.h"
#include "outputwriter.h"
#include "mesh.h"
#include "time.h"

#include "hdf5.h"
#include "hdf5_hl.h"

/* Compute the weight for a reaction between dof i and
 dof j. We will test different startegies here to find out what is
 best and simplest. 0.5 (should) give the method from Fange et. al. */
inline double weight(int sv1,int sv2)
{	
	return 0.25;	
}


int main(int argc, char *argv[])

{
	
	char *infile,*outfile;
	int i, nt=1;

	if (argc < 3){
		printf("To few arguments to nsm.");
	    exit(-1);	
	}
	
	/* Input file. */
	infile  = argv[1];
	/* Output file (or directory). */
	outfile = argv[2];
	
	/* Read model specification */
	urdme_model *model;
	model = read_model(infile);
	model->infile = infile;
	
	
	if (model == NULL){
		printf("Fatal error. Couldn't load model file or currupt model file.");
		return(-1);
	}
	
	/* Check model file for optional report level and seed. */ 
	MATFile *input_file;
	input_file = matOpen(infile,"r"); 
	mxArray *temp;
	
	if (input_file == NULL){
		printf("Failed to open mat-file.\n");
		//return NULL;	
		exit(-1);
	}
	/*temp = matGetVariable(input_file, "report");
	
	if (temp != NULL) 
		report_level = (int) mxGetScalar(temp);
	else
		if (nt > 1)
			report_level=0;
	    else 
			report_level=1;*/
	
	model->num_extra_args=1;
	model->extra_args=(void **)malloc(model->num_extra_args*sizeof(void *));
	
	/* Look for seed */
	temp = matGetVariable(input_file, "seed");
	
	long int seed;
	if (temp != NULL) {
		seed = (long int) mxGetScalar(temp);
		srand48(seed);
	}
	
	/* If seed is provided as a parameter, it takes precedence. */
	if (argc > 3) {
		srand48((long int)atoi(argv[3]));  
	}
	else {
	      /* Not a foolproof solution */
	      srand48((long int)time(NULL)+(long int)(1e9*clock()));
	}
	

	/* Look for an optional parameter matrix. */
	const double *matfile_parameters; 
	int mpar = 0;
	int npar = 0;
	int param_case=0;
	
	temp = matGetVariable(input_file, "parameters");
        if (temp != NULL) {
		matfile_parameters = (double *)mxGetPr(temp);
		mpar = mxGetM(temp);
		npar = mxGetN(temp); 
	}
	
	/* Look if a parameter case if supplied as a parameter. */
	if (argc > 4) {
	    param_case = (int)atoi(argv[4]);
	}
	
	if (param_case > npar ) {
		printf("nsmcore: Fatal error, parameter case is larger than n-dimension in parameter matrix.\n");
		exit(-1);
	}
	
	/* Create global parameter variable for this parameter case. */
    parameters = (double *)malloc(mpar*sizeof(double));
	memcpy(parameters,&matfile_parameters[npar*param_case],mpar*sizeof(double));
	
	fem_mesh *mesh;
	mesh = (fem_mesh *)malloc(sizeof(fem_mesh));
	
	temp = matGetVariable(input_file, "p");
	double *mesh_p;
	if (temp != NULL) 
		mesh_p = (double *)mxGetPr(temp);
	else{
		printf("Error: No mesh vertex information in model file. ");
		return NULL;
	}
	
	mesh->p = mesh_p;
	
	/* Adjacency matrix */
	size_t *jcK,*irK;
	double *prK;
	temp = matGetVariable(input_file, "K");
	if (temp != NULL){ 
		irK = (size_t *)mxGetIr(temp);
		jcK = (size_t *)mxGetJc(temp);
		prK = (double *)mxGetPr(temp);
	}
	else{
		printf("Error: No connectivity matrix in model file. ");
		exit(-1);
	}
	
	mesh->jcK=jcK;
	mesh->irK=irK;
	mesh->prK=prK;
	
	model->extra_args[0]=(void *)mesh;
	
	matClose(input_file);
	
	/* Allocate memory to hold nt solutions. */
	init_sol(model,nt);
	
   /* Get a writer to store the output trajectory on a hdf5 file. */
    urdme_output_writer *writer;
    writer = get_urdme_output_writer(model,outfile);
    
    nem(model, writer);
	//}

	/* Print result to file(s)*/
	/*int didprint;
	didprint=dump_results(model, outfile,"single");
	if (!didprint){
		printf("Error: Failed to print results to file.\n");
		exit(-1);
	}*/
    /* Write the timespan vector to the output file */
    write_tspan(writer,model);
	
    destroy_output_writer(writer);
	destroy_model(model);
	free(parameters);
	
	return(0);
	
}

/* Temporary wrapper for nsm-core solver (while we wait for me to rewrite it on the new form). */
void *nem(void *data, urdme_output_writer *writer){
	
	/* Unpack input */
	urdme_model* model;
	model = (urdme_model *)data;
	int Ndofs, *U;
	
	/* nsm_core uses a report function with optional report level. This is
	 passed as the first extra argument. */ 
	//int report_level = *(int *)model->extra_args[0];
	
	/* Output array (to hold a single trajectory) */
	Ndofs = model->Ncells*model->Mspecies;
	U = (int *)malloc(model->tlen*Ndofs*sizeof(int));
	
	/* !!!!!!!!!!!!!!!!!!!!!!!! */
	model->Mreactions = number_of_reactions();
	model->Mspecies   = number_of_species();
	
	/* Call the core simulation routine.  */
	nemcore((fem_mesh *)model->extra_args[0],model->irD, model->jcD, model->prD, model->u0,
			 model->irN, model->jcN, model->prN, model->irG,
			 model->jcG, model->tspan, model->tlen, 
             model->vol, model->data, model->sd, model->Ncells,
			 model->Mspecies, model->Mreactions, model->dsize, writer);
	
	
}


 

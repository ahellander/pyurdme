/* A. Hellander and B. Drawert. */
#include <string.h>
#include "urdmemodel.h"
//#define OUTPUT_MAT
#define OUTPUT_HDF5
#ifdef OUTPUT_HDF5
#include "hdf5.h"
#include "hdf5_hl.h"
#endif

// This is necessary as 'parameters' is extern in "propensities.h" (where it is declared).
// It must be defined (set to a value) in one and only one '.o' file
#include "propensities.h"
double *parameters = NULL;

/*
 
 URDME_LIBMAT is defined at compile time, 
 and dictates wether URMDE will use Mathworks' MATLAB libraries
 for handling .mat files, or our own implmentation.
 
 The Matlab libraries are currently needed for 
 sparse output.  
 
 */

#ifdef URDME_LIBMAT
#include "read_matfile.h"
#else
#include "mat.h" 
#include "mex.h"
#include "matrix.h"
#endif


/*
 
   Read model input file (.mat format) and initialize urdme model struct. 
   If any of the required fields (for the nsm solver) is missing, 
   it returns a NULL-pointer.
   
*/

urdme_model *read_model(char *file)
{
	
	int i,Ndofs;
	
	mxArray *D,*N,*G,*vol,*u0,*tspan,*data,*sd,*K;
	
	urdme_model *model;
	model = (urdme_model *)malloc(sizeof(urdme_model));
	
	/* Open mat-file (read-only) */
	MATFile *input_file;
	input_file = matOpen(file,"r"); 
	
	if (input_file == NULL){
		printf("Failed to open mat-file.\n");
		return NULL;	
	}

	/* Get D-matrix. */
	D = matGetVariable(input_file, "D");
	if (D == NULL){
		printf("The diffusion matrix D is missing in the model file.\n");
		return NULL;
	}

    Ndofs  = mxGetN(D);
    size_t *mxirD;
    size_t *mxjcD;
    double *mxprD;
    mxirD = mxGetIr(D);
    mxjcD = mxGetJc(D);
    mxprD = mxGetPr(D);
    int nnzD = mxGetNzmax(D);
    
    model->jcD = (size_t *)malloc((Ndofs+1)*sizeof(size_t));
    memcpy(model->jcD,mxjcD,(Ndofs+1)*sizeof(size_t));
    model->irD = (size_t*)malloc(nnzD*sizeof(size_t));
    memcpy(model->irD,mxirD,nnzD*sizeof(size_t));
    model->prD = (double *)malloc(nnzD*sizeof(double));
    memcpy(model->prD,mxprD,nnzD*sizeof(double));
    mxDestroyArray(D);

    /* initial condition */
	u0 = matGetVariable(input_file, "u0");
	if (u0 == NULL){
		printf("Initial condition (u0) is missing in the model file.\n");
		return NULL;
	}
	/* Typecast */
    model->Mspecies   = (int) mxGetM(u0);

	int *u0int;
	u0int = (int*)malloc(Ndofs*sizeof(int));
	double *u0temp;
	u0temp = mxGetPr(u0);
	for (i=0;i<Ndofs;i++)
		u0int[i] = (int) u0temp[i];
	model->u0 = u0int;
    mxDestroyArray(u0);
    
	/* Stoichiometry matrix */
	N = matGetVariable(input_file, "N");
    if (N == NULL){
		printf("The stoichiometry matrix is missing in the model file.\n");
		return NULL;
	}
    
    if (mxIsSparse(N)){
        /* This is always true if the matrix is created by the Matlab interface. */
        
        /* Typecast to int */
        double *tempN;
        tempN = mxGetPr(N);
        int nnzN = (int) mxGetNzmax(N);
        int *prN;
        prN = (int*)malloc(nnzN*sizeof(int));
        for (i=0;i<nnzN;i++)
            prN[i] = (int) tempN[i];
        model->prN = prN;
        
        
        model->irN = (size_t *)malloc(nnzN*sizeof(size_t));
        memcpy(model->irN,(size_t *)mxGetIr(N),nnzN*sizeof(size_t));
        model->jcN = (size_t *)malloc((mxGetN(N)+1)*sizeof(size_t));
        memcpy(model->jcN,(size_t *)mxGetJc(N),(mxGetN(N)+1)*sizeof(size_t));
    }
    else{
        if (mxGetN(N)>0||mxGetM(N)>0){
            printf("The stoichiometry matrix must be a sparse CCS matrix.");
            return NULL;
        }
        else{
            /* This little pice of code is needed to cover the case where
               there are no reactions and the matrix is passed as an empty dense
               array, as is necessary if it is generated by pyurdme using scipy. */
            model->irN = (size_t *)malloc(0*sizeof(size_t));
            model->jcN = (size_t *)malloc(0*sizeof(size_t));
            model->prN = (int *)malloc(0*sizeof(int));
        }
    }
       model->Mreactions = (int) mxGetN(N);
       model->Ncells=Ndofs/model->Mspecies;
   
    mxDestroyArray(N);



    /* Connectivity matrix */
    K = matGetVariable(input_file, "K");
    if (K == NULL){
    	printf("The Connectivity matrix K is missing in the model file.\n");
    	return NULL;
    }

    int nK = mxGetN(K);
    size_t *mxirK;
    size_t *mxjcK;
    double *mxprK;
    mxirK = mxGetIr(K);
    mxjcK = mxGetJc(K);
    mxprK = mxGetPr(K);
    int nnzK = mxGetNzmax(K);
    
    model->jcK = (size_t *)malloc((nK+1)*sizeof(size_t));
    memcpy(model->jcK,mxjcK,(nK+1)*sizeof(size_t));
    model->irK = (size_t*)malloc(nnzK*sizeof(size_t));
    memcpy(model->irK,mxirK,nnzK*sizeof(size_t));
    model->prK = (double *)malloc(nnzK*sizeof(double));
    memcpy(model->prK,mxprK,nnzK*sizeof(double));
    mxDestroyArray(K);


    /* Volume vector */
	vol = matGetVariable(input_file,"vol");
	if (vol == NULL){
		printf("The volume vector is missing in the model file.\n");
		return NULL;
	}
    model->vol = (double *)malloc(model->Ncells*sizeof(double));
    memcpy(model->vol,(double *)mxGetPr(vol),model->Ncells*sizeof(double));
    mxDestroyArray(vol);

	/* Dependency graph */
	G = matGetVariable(input_file, "G");
    if (G == NULL){
		printf("The dependency graph (G) is missing in the model file.\n");
		return NULL;
	}
    if (mxIsSparse(G)){
        int nnzG = mxGetNzmax(G);
        model->irG = (size_t *)malloc(nnzG*sizeof(size_t));
        memcpy(model->irG,(size_t *)mxGetIr(G),nnzG*sizeof(size_t));
        model->jcG = (size_t *)malloc((mxGetN(G)+1)*sizeof(size_t));
        memcpy(model->jcG,(size_t *)mxGetJc(G),(mxGetN(G)+1)*sizeof(size_t));
        mxDestroyArray(G);
    }
    else{
        model->irG = (size_t *)malloc(0*sizeof(size_t));
        model->irG = (size_t *)malloc(0*sizeof(size_t));
        model->jcG = (size_t *)malloc(0*sizeof(size_t));
    }
	

    /* time span */
	tspan = matGetVariable(input_file, "tspan");
	if (tspan == NULL){
		printf("Time span (tspan) is missing in the model file.\n");
		return NULL;
	}
	model->tlen = mxGetNumberOfElements(tspan);
    model->tspan = (double *)malloc(model->tlen*sizeof(double));
    memcpy(model->tspan,(double *)mxGetPr(tspan),model->tlen*sizeof(double));
    mxDestroyArray(tspan);
    
	/* Subdomain vector */
	sd = matGetVariable(input_file, "sd");
	if (sd == NULL){
		printf("Subdomain vector (sd) is missing in the model file.\n");
		return NULL;
	}
	
	/* typecast */
	int *sdint = (int*)malloc(model->Ncells*sizeof(int));
	double *sdtemp;
	sdtemp = mxGetPr(sd);
	for (i=0;i<model->Ncells;i++)
		sdint[i]=(int) sdtemp[i];
    
	model->sd = sdint;
    mxDestroyArray(sd);

    
	/* data matrix */
	data = matGetVariable(input_file, "data");
	if (data == NULL){
		printf("Data matrix (data) is missing in the model file.\n");
		return NULL;
	}
	
	model->dsize = mxGetM(data);
    model->data = (double *)malloc(model->dsize*model->Ncells*sizeof(double));
    memcpy(model->data,(double *)mxGetPr(data),model->dsize*model->Ncells*sizeof(double));
    mxDestroyArray(data);
	
    /* Maximum number of solutions defaults to one. */
    model->nsolmax = 1;
    
	matClose(input_file);	
	
	return model;
	
	
}

//--------------------------------------------
void read_solution(urdme_model *model, char*file){
	mxArray *tspan,*U;
	MATFile *input_file;
	/* Open mat-file (read-only) */
	input_file = matOpen(file,"r"); 
	if (input_file == NULL){
		printf("Failed to open mat-file '%s'.\n",file);
		return;	
	}
    printf("read in '%s'\n'",file);
	
	/* Get U-matrix. */
    init_sol(model,1);
    U = matGetVariable(input_file, "U");
	if (U == NULL){
        printf("No U solution variable in mat-file '%s'\n.",file);
        return;
    }
    int Usz = mxGetNumberOfElements(U);
	model->U[0] =(int *)malloc(Usz*sizeof(int));
    double*Utmp = mxGetPr(U);
	model->nsol = 1;
    int i;
    for(i=0;i<Usz;i++){
        model->U[0][i] = (int) Utmp[i];
    }

	/* time span (optional) */
	tspan = matGetVariable(input_file, "tspan");
	if (tspan != NULL){
        model->tspan = mxGetPr(tspan);
        model->tlen = mxGetNumberOfElements(tspan);
	}
}


/* Utility function. Initialize the solution field of the urdme_model struct. */
void init_sol(urdme_model *model, int nsolmax)
{
	model->nsolmax = nsolmax;
	model->nsol = 0;
	
}

/* Free model. Always call this destructor before exiting. */
int destroy_model(urdme_model *model)
{
    
    int i;
    
    /* D-matrix */
    free(model->irD);
    free(model->jcD);
    free(model->prD);
    
    /* Volume vector */
    free(model->vol);
    
    /* Dependency graph */
    free(model->jcG);
    free(model->irG);
    
    /* Stoichiometry matrix */
    free(model->irN);
    free(model->jcN);
    free(model->prN);
	
    /* Connectivity matrix */
    free(model->irK);
    free(model->jcK);
    free(model->prK);
	
    free(model->tspan);
    free(model->u0);
	free(model->sd);
    free(model->data);
    
    for (i=0; i<model->num_extra_args; i++)
      if (model->extra_args[i]!=NULL){
        free(model->extra_args[i]);
	}
    
	free(model);
    
	return 0;
}


/***** DEPRECATED *****/

/*
 Print the result trajectory that is attached to urdme_model to a mat/hdf5 file.
 The output trajectory will be stored in a variable/dataset named "U".
 */
int dump_results(urdme_model* model, char *filename, char *type){
    /* IN PYURDME, THIS IS NEVER CALLED. IT IS ALWAYS USING HDF5. */
    
    int Ndofs,tlen,nsol;
    int *U;
    U = model->U[0];
    
    
    Ndofs = model->Ncells*model->Mspecies;
    tlen = model->tlen;
    nsol = model->nsol;
    
    
#ifdef OUTPUT_MAT
    
    
    MATFile *output;
    mxArray *Uout;
    
    
#ifndef URDME_OUTPUT_SPARSE
    Uout = mxCreateDoubleMatrix(Ndofs,tlen,mxREAL);
    double *data;
    data = mxGetPr(Uout);
    /* Dense output data */
    for (i=0;i<Ndofs*model->tlen;i++){
        data[i]=(double) U[i];
    }
    
#else
    /*
     If sparse output, we write the matrix on dictionary of keys (DOK) format.
     The application using this matrix will then have to convert to whatever format needed.
     the keys are (iU,jU) and values sU. To instantiate a sparse matrix in Matlab, load the output file and do
     >> U = sparse(iU,jU,sU,mU,nU);
     */
    
    /* Count number of non-zeros */
    int nnz = 0,nnz_col;
    for (i=0; i<Ndofs*tlen; i++) {
        if (U[i]>0.0)
            nnz++;
    }
    
    mxArray* iU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* jU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    mxArray* sU = mxCreateDoubleMatrix(nnz,1,mxREAL);
    
    // Dimesions of the matrix, mU (row) and nU(col). In-house matlib does not have mxCreateScalar.
    mxArray* mU = mxCreateDoubleMatrix(1,1,mxREAL);
    double *dim;
    dim = mxGetPr(mU);
    dim[0] = (double)Ndofs;
    mxArray* nU = mxCreateDoubleMatrix(1,1,mxREAL);
    dim = mxGetPr(nU);
    dim[0] = (double)tlen;
    
    double *iUdata = mxGetPr(iU);
    double *jUdata = mxGetPr(jU);
    double *sUdata = mxGetPr(sU);
    
    nnz = 0;
    for (j=0; j<tlen; j++){
        for (i=0;i<Ndofs;i++){
            if (U[j*Ndofs+i]>0.0){
                /* NOTE THE +1 HERE, MATLAB SPECIFIC INDEXING. */
                iUdata[nnz] = (double)(i+1);
                jUdata[nnz] = (double)(j+1);
                sUdata[nnz] = (double)U[j*Ndofs+i];
                nnz++;
            }
        }
    }
    
#endif
    
    
    mxArray *Tspanout;
    Tspanout = mxCreateDoubleMatrix(1, model->tlen,mxREAL);
    double *tdata;
    tdata = mxGetPr(Tspanout);
    
    for(i=0;i<model->tlen;i++)
        tdata[i]=model->tspan[i];
    
    
#ifdef URDME_LIBMAT
    output = matOpen(filename,"w"); // using built in URDME mat read/write
#else
    output = matOpen(filename,"w7.3"); // Possibly large files
#endif
    
    if (output == NULL){
        printf("Error: Could not write to output file: %s\n",filename);
        return(-1);
    }
    
#ifndef URDME_OUTPUT_SPARSE
    matPutVariable(output,"U",Uout);
    matPutVariable(output,"tspan",Tspanout);
#else
    matPutVariable(output,"iU",iU);
    matPutVariable(output,"jU",jU);
    matPutVariable(output,"sU",sU);
    matPutVariable(output,"mU",mU);
    matPutVariable(output,"nU",nU);
    matPutVariable(output,"tspan",Tspanout);
#endif
    
    matClose(output);
    
#ifndef URDME_OUTPUT_SPARSE
    mxDestroyArray(Uout);
    mxDestroyArray(Tspanout);
#else
    mxDestroyArray(iU);
    mxDestroyArray(jU);
    mxDestroyArray(sU);
    mxDestroyArray(Tspanout);
    mxDestroyArray(mU);
    mxDestroyArray(nU);
#endif
    
#elif defined(OUTPUT_HDF5)
    /* Write the result as a HDF5 dataset/file */
    
    hid_t h5_output_file;
    herr_t status;
    h5_output_file = H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
    if (h5_output_file == NULL){
        printf("Failed to write U matrix HDF5 file.");
        exit(-1);
    }
    
    hsize_t dims[2]; /* dataset dimensions */
    dims[0] = Ndofs;
    dims[1] = tlen;
    
    // Write the result matrix to the file
    status = H5LTmake_dataset(h5_output_file,"/U",2,dims,H5T_NATIVE_INT,U);
    if (status != 0){
        printf("Failed to write U matrix HDF5 file.");
        exit(-1);
    }
    
    // Write tspan to the file
    dims[0] = 1;
    dims[1] = tlen;
    status = H5LTmake_dataset(h5_output_file,"/tspan",2,dims,H5T_NATIVE_DOUBLE,model->tspan);
    if (status != 0){
        printf("Failed to write tspan vector HDF5 file.");
        exit(-1);
    }
    
    status = H5Fclose(h5_output_file);
    if (status != 0){
        printf("Failed to close output HDF5 file.");
        exit(-1);
    }
    
#endif
    
    return 1;
    
}

/*
    Print all data from model to stdout, for debugging
*/
void debug_print_model(urdme_model* model){
    int i,j;
	/* Constants */
	printf("Mspecies=%i\n",model->Mspecies);
	printf("Mreactions=%i\n",model->Mreactions);
	printf("Ncells=%i\n",model->Ncells);

    int Ndofs = model->Ncells * model->Mspecies;
    
    for (i=0; i<Ndofs; i++) {
        for (j=model->jcD[i]; j<model->jcD[i+1]; j++){
            printf("D: %i\t%i\t%f\n",j,(int)model->irD[j],model->prD[j]);
        }
    }
    printf("end D\n");
	
	/* Volume vector */
    printf("vol:");
    for(i=0;i<model->Ncells;i++){
        printf("%f ",model->vol[i]);
    }
    printf("\n");

    
    for(i=0;i<model->Mspecies;i++){
        for (j=model->jcN[i]; j<model->jcN[i+1]; j++){
            printf("N: %i\t%lu - %lu\t%lu\t%i\n",j,(long unsigned)model->jcN[i],(long unsigned)model->jcN[i+1],(long unsigned)model->irN[j],model->prN[j]);
        }
    }
    printf("end N\n");

    for(i=0;i<model->Mreactions;i++){
        for (j=model->jcG[i]; j<model->jcG[i+1]; j++){
            printf("G %i\t%lu\n",j,(long unsigned)model->irG[j]);
        }
    }
    printf("end G\n");
	
	/* Initial condition */
    printf("u0:");
    for(i=0;i<Ndofs;i++){
        printf("%i ",model->u0[i]);
    }
    printf("\nend u0\n");
	
	/* Time span 
	int tlen;
	double *tspan;*/
    printf("tspan(%i):",model->tlen);
    for(i=0;i<model->tlen;i++){
        printf("%f ",model->tspan[i]);
    }
    printf("\nend tspan\n");
	
	/* Subdomain vector 
	int *sd;*/
    printf("sd:");
    for(i=0;i<model->Ncells;i++){
        printf("%i ",model->sd[i]);
    }
    printf("\nend sd\n");
	
	/* Data vector 
	int dsize;
	double *data;*/
    printf("data:");
    for(i=0;i<model->Ncells*model->dsize;i++){
        printf("%f ",model->data[i]);
    }
    printf("\nend data\n");
    //------------
}

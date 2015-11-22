#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vertex.h"
#include "species.h"

/* Initialize the vertex list. */
vertex **initialize_vertices(int number_of_vertices,int Mspecies)
{
	int i;
	int maxnredges=30;
	/* FIXME */
	int ndim=2;
	
	vertex *v;
	vertex **vlist = (vertex **)malloc(number_of_vertices*sizeof(vertex *));
	
	for (i=0; i<number_of_vertices; i++) {
		
		vlist[i] = (vertex *)malloc(sizeof(vertex));
		v = vlist[i];
		
		/* Index. */
		v->vi=i;
		
		/* Coordinate */
		v->x =(double *)calloc(ndim,sizeof(double));
		
        /* Volumes */
        v->vol =(double *)calloc(Mspecies,sizeof(double));
		
		v->total_rate=0.0;
		v->total_local_rate=0.0;
		v->total_edge_rate = 0.0;
		v->reaction_weight=0.0;
		
		/* Allocate space to hold edge_lists. Currently wastes a lot of space, FIXME LATER!!! */
		v->number_of_outgoing_edges=0;
		v->outgoing_edges = (int *)malloc(maxnredges*sizeof(int));
		v->number_of_incoming_edges=0;
		v->incoming_edges = (int *)malloc(maxnredges*sizeof(int));
		
	}
		
	return vlist;
}

/* Utility function to compute partial values for activity_coefficient_2D */
void compute_eta_2D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,double *eta)
{
    
    int s;
    /* Free volume in a voxel in 2D */
    double pi = 3.1416;
    double eta0 = 0.0;
    double eta1 = 0.0;
    double eta2 = 0.0;
    double rho=0.0;
    
    double sigma;
    
    for (s=0;s<Mspecies;s++){
        rho = xx[vtx->vi*Mspecies+s]/vtx->vol[s];
        sigma = list_of_species[s]->sigma;
        eta0 += rho;
        eta1 += 2.0*sigma*rho;
        eta2 += 4.0*sigma*sigma*rho;
    }
    
    eta0=eta0*(pi/4.0);
    eta1=eta1*(pi/4.0);
    eta2=eta2*(pi/4.0);
    eta[0]=eta0;
    eta[1]=eta1;
    eta[2]=eta2;
}

/* Utility function to compute partial values for activity_coefficient_2D */
void compute_eta_3D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,double *eta)
{
    
    
    int s;
    /* Free volume in a voxel in 2D */
    double pi = 3.1416;
    
    double eta0 = 0.0;
    double eta1 = 0.0;
    double eta2 = 0.0;
    double eta3 = 0.0;
    double a;
    double rho;
    
    double v0,sigma;
    v0=vtx->vol[0];
    
    for (s=0;s<Mspecies;s++){
        rho = xx[vtx->vi*Mspecies+s]/v0;
        sigma = list_of_species[s]->sigma;
        eta0 += rho;
        a = sigma*rho;
        eta1 += 2.0*a;
        eta2 += 4.0*sigma*a;
        eta3 += 8.0*sigma*sigma*a;
    }
    
    eta0=eta0*(pi/6.0);
    eta1=eta1*(pi/6.0);
    eta2=eta2*(pi/6.0);
    eta3=eta3*(pi/6.0);
    
    eta[0]=eta0;
    eta[1]=eta1;
    eta[2]=eta2;
    eta[3]=eta3;
}

/* Calculate the activity coefficient in the 3D case based on the formulas in the paper . */
double activity_coefficient_3D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,int spec,double *eta)
{
    
    
    double gamma;
    /* Free volume in a voxel in 2D */
    
    double sigma =list_of_species[spec]->sigma;
    double x = 1.0-eta[3];
    if (x <=0.0){
        gamma = -INFINITY;
        return gamma;
    }
    gamma = -log(x) + 6.0*eta[2]/x*sigma + (12.0*eta[1]/(x*x)+18.0*eta[2]*eta[2]/(x*x))*sigma*sigma + (8.0*eta[0]/x + 24.0*eta[1]*eta[2]/(x*x)+24.0*eta[2]*eta[2]*eta[2]/(x*x*x))*sigma*sigma*sigma;
    gamma = exp(gamma);
    return gamma;
    
}



/* Calculate the activity coefficient in the 2D case based on the formulas in the paper . */
double activity_coefficient_2D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,int spec,double *eta)
{
    
    
    double gamma;
    /* Free volume in a voxel in 2D */
    
    double sigma =list_of_species[spec]->sigma;
    /*int sums=0;
    for (int s=0;s<Mspecies;s++)
        sums += xx[vtx->vi*Mspecies+s];*/
    double x = 1.0-eta[2];
    if (x <=0.0){
        gamma = -INFINITY;
        return gamma;
    }
    gamma = -log(x) + 4.0*eta[1]/x*sigma + (4.0*eta[0]/x+4.0*eta[1]*eta[1]/(x*x))*sigma*sigma;
    //gamma_test=gamma;
    gamma = exp(gamma);
    /*if (sums==0&&gamma>1.0) {
        exit(-1);
    }*/
    //return 1.0;
    /*if (isnan(gamma)){
        printf("gamma is nan\n,xx: %i, eta: %f %f %f x,gamma test: %f %f\n",xx[vtx->vi*Mspecies+spec],eta[0],eta[1],eta[2],x,gamma_test);
        exit(-1);
    }*/
    return gamma;
    
}

/* Calculate the naive free volume of the voxel as the voxel volume - the volume of 
   all the spherical particles. */
double available_volume(vertex *vtx, int *xx,int Mspecies,species **list_of_species)
{
    
    int s;
    
    /* Free volume in a voxel in 2D */
    
    /* packing fraction? */
    double eta = 1.45;
    /* Compute the free volumes of the vertices */
    double v0,sigma;
    v0=vtx->vol[0];
    
    for (s=0;s<Mspecies;s++){
        sigma = list_of_species[s]->sigma;
        v0-=eta*4.0*sigma*sigma*3.1416*xx[vtx->vi*Mspecies+s];
    }
  
    return fmax(v0,0.0);
}

void free_vertex(vertex *v)
{
	
	free(v->x);
	free(v->outgoing_edges);
	free(v->incoming_edges);
	free(v);
	
}

void free_list_of_vertices(vertex **vlist,int number_of_vertices)
{
	int i;
	for (i=0; i<number_of_vertices; i++) {
		free_vertex(vlist[i]);
	}
	free(vlist);
}



void print_vertex(vertex *v)
{
 
   printf("Index: %i\n",v->vi);
   printf("Number of outgoing edges: %i\n",v->number_of_outgoing_edges);	 
   printf("Number of incoming edges: %i\n",v->number_of_incoming_edges);	
   printf("\n");	
}

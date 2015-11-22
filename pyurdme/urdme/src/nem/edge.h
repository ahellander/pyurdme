/* A. Hellander, 2012-02-01. */

#ifndef EDGE_H	
#define EDGE_H		

#include "vertex.h"
#include "reactions.h"
#include "species.h"

/* Extended edge data structure for use in the next edge method. */
typedef struct
{
	
	/* Self index */
	int index;
	
	/* Vertices */
	vertex *v1,*v2; 	
    
	/* Adjacency map?? */
	//int number_of_adjacent_edges;
	//int *list_of_adjacent_edges;
	
	/* Length */
	double h;
    
	/* Weight (stiffness) of the edge. */
	double weight;
	
	/* Diffusion weights (we will get this from the stiffness matrix). 
	   d will be a Mspecies*2 matrix. */
	double *d;	
 	
	double reaction_weight;
	
	/* Reaction weights (how we will get those is yet to be determinied).
	   Will also be a Mspecies*2 matrix.*/
	double *w;	
 	
    /* Total intensity for this edge. */
	double total_reaction_rate;
	double total_diffusion_rate;
	double total_rate;

	
} edge;

edge **initialize_edges(int number_of_edges, int number_of_vertices, vertex **list_of_vertices,size_t *jcK,size_t *irK, double *prK,int Mspecies);
void edge_initialize_h(edge **list_of_edges,int number_of_edges,double *p, int dim);


double edge_compute_nonlocal_reaction_rates(edge *e, int *xx, reaction **list_of_reactions,species **list_of_species,double *vol,int Mspecies,int Mreactions);

double edge_compute_diffusion_rates(edge *e, int *xx,int Mspecies,species **list_of_species,int dim);

int edge_execute_next_event(edge *e);
double edge_update_rates(edge *e);

double node_compute_reaction_rate(vertex *vtx,int *xx,reaction *r,species **list_of_species,double *vol,int Mspecies);
//double node_compute_reaction_rates(int node,int *xx,reaction **r);

int edge_execute_next_reaction_event(edge *e,int *xx,int Mreactions, reaction **list_of_reactions,int Mspecies,species **list_of_species,double *vol);
int edge_execute_next_diffusion_event(edge *e,int *xx, int Mspecies, species **list_of_species,int dim);

double edge_compute_chemical_potential(edge *e,int *xx,int Mspecies,species **list_of_species);
void edge_diffusion_excluded_volume_correction(edge *e,int *xx,int Mspecies,species **list_of_species,double *delta_mu,int dim);


void free_edge(edge *e);
void free_list_of_edges(edge **elist,int number_of_edges);

double sample_next_edge_event(edge *e);

#endif
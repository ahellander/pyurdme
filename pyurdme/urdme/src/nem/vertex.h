#ifndef VERTEX_H
#define VERTEX_H	

#include "species.h"

typedef struct
{
  /* Vertex index */
  int vi;
    
  /* Coordinate */
  double *x;
  
  /* Volumes (1xMspecies) */
  double *vol;
  	
  /* List of adjacent edges (edges) */
  int number_of_incoming_edges;
  int number_of_outgoing_edges;	
  int *outgoing_edges;
  int *incoming_edges;
	
  /* List of adjacent nodes */
  int number_of_adjacent_vertices;	
  int *adjacent_vertices;
	
  /* Reaction weight for bimolecular reactions. */
  double reaction_weight;	
	
  /* Total event rate for this vertex. */
  double total_rate;	
  /* Total local reaction rate */
  double total_local_rate;	
  /* Sum of total rates for all adjacent edges. */	
  double total_edge_rate;
	
} vertex;


vertex **initialize_vertices(int number_of_vertices,int Mspecies);

void print_vertex(vertex *v);
void free_vertex(vertex *v);
void free_list_of_vertices(vertex **vlist,int number_of_vertices);
double available_volume(vertex *vtx, int *xx,int Mspecies,species **list_of_species);
double activity_coefficient_2D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,int spec,double *eta);
void compute_eta_2D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,double *eta);
double activity_coefficient_3D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,int spec,double *eta);
void compute_eta_3D(vertex *vtx, int *xx,int Mspecies,species **list_of_species,double *eta);



#endif
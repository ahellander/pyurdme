/* A. Hellander, 2012-02-01. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vertex.h"
#include "edge.h"
#include "reactions.h"
#include "species.h"

//inline double weight(int dof1,int dof2){
//	return 0.0;
//}

/* Preprocess the model and create the edge-list with all the proper weights (from the
   weighted adjacency matrix). */
edge **initialize_edges(int number_of_edges,int number_of_vertices,vertex** list_of_vertices,size_t *jcK,size_t *irK,double *prK,int Mspecies)
{

	edge **list_of_edges;
	// Number of unique edges?? From (nnz-Ndofs)/Mspecies. 
	list_of_edges = (edge **)malloc(number_of_edges*sizeof(edge *));
	
	edge *et;
	int i,j,e;
	vertex *v1,*v2;
    
	/* Loop over all elements in the weighted adjacency matrix and initialize
	   the edges. This is not a robust solution; may fail if all species have diffusion
	   constant zero at a node. */
	
	e=0;
	for (i=0; i<number_of_vertices; i++) {
		
		v1 = list_of_vertices[i];			
		
		for (j=jcK[i]; j<jcK[i+1]; j++) {
			
			if (irK[j]!=i) {
				
				/* Create edge */
				list_of_edges[e]=(edge *)malloc(sizeof(edge)); 
			    et=list_of_edges[e];
				et->index=e;
				
				/* Vertices for this edge */
			    v2=list_of_vertices[irK[j]];
				et->v1 = v1;
				et->v2 = v2;
				
				/* Weigth (stiffness) for the edge. */
				et->weight = prK[j]/v1->vol[0];
				
				/* Add edge to vertex adjacency lists. This edge is an outgoing edge for vertex v1 ... */
				v1->outgoing_edges[v1->number_of_outgoing_edges]=e;
				v1->number_of_outgoing_edges++;
				
				/* and an incoming edge for v2. */
				v2->incoming_edges[v2->number_of_incoming_edges]=e;
				v2->number_of_incoming_edges++;
				
				et->reaction_weight=0.0;
				et->total_rate=0.0;
				et->total_reaction_rate=0.0;
				et->total_diffusion_rate=0.0;
				
				/* Length of edge */
				//x1 = v1->x;
				//x2 = v2->x;
				//et->h=sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]));
				
				e++;
			    
				
			}
			
		}
	
	}

	return list_of_edges;
}

void edge_initialize_h(edge **list_of_edges,int number_of_edges,double *p, int dim)
/* Compute the length of all edges */
{

	int j,e;
	vertex *v1,*v2;
	double *x,*y;
	
	for (e=0; e<number_of_edges; e++) {
	
		list_of_edges[e]->h=0.0;
		
		v1 = list_of_edges[e]->v1;
		v2 = list_of_edges[e]->v2;
		
		x = v1->x;
		y = v2->x;
		
		for (j=0; j<dim; j++) {
			list_of_edges[e]->h += (x[j]-y[j])*(x[j]-y[j]);
		} 
		
		list_of_edges[e]->h = sqrt(list_of_edges[e]->h);
		
	}
	
}

void free_edge(edge *e)
{
	free(e);	
}

/* Deallocate a list of edges. */
void free_list_of_edges(edge **edges,int number_of_edges)
{
  	
  	int e;
	for (e=0; e<number_of_edges; e++) {
		free_edge(edges[e]);
	}
	free(edges);

}

double node_compute_reaction_rate(vertex *vtx,int *xx,reaction *r,species **list_of_species,double *vol,int Mspecies){
    /* Compute the total reaction rate at a vertex */
	double a;
	int node = vtx->vi;
    a = evaluate_propensity(r,xx,list_of_species,node,vol[node],node,vol[node],Mspecies);
	if (r->order==2) {
		a*=vtx->reaction_weight;
	}
	return a;
}

double edge_compute_nonlocal_reaction_rates(edge *e, int *xx, reaction **list_of_reactions,species **list_of_species,double *vol,int Mspecies, int Mreactions){
/* Compute the non-local reaction rates for an edge. */

	double a0=0.0;
	double a1;
	reaction *reac;
	int r;
	vertex *v1,*v2;
	
	v1 = e->v1;
	v2 = e->v2;
	
	for (r=0; r<Mreactions; r++) {
		
		reac = list_of_reactions[r];
		/* If a bimolecular reaction, compute all possible combinations for reactions along this edge.  */
		if (reac->order == 2) {
			
			/* A_1 + B_2 -> C_1 */
		   	a1  = evaluate_propensity(reac,xx,list_of_species,v1->vi,vol[v1->vi],v2->vi,vol[v2->vi],Mspecies);	
			a0 += e->reaction_weight*a1;
			/* A_1 + B_2 -> C_2 */
			a0 += e->reaction_weight*a1;
			
		}
	}
    
	return a0;
}


/* Find the next non-local reaction event along this edge and update the state accordingly. */
int edge_execute_next_reaction_event(edge *e, int *xx, int Mreactions, reaction **list_of_reactions,int Mspecies, species **list_of_species,double *vol){
	
	double a1;
	reaction *reac;
	int r;
	vertex *v1,*v2;
	
	v1 = e->v1;
	v2 = e->v2;
	
	double ur;
	ur = (1.0-drand48())*e->total_reaction_rate;
	double cumsum = 0.0;
	
	for (r=0; r<Mreactions; r++) {
		
		reac = list_of_reactions[r];
		/* If a bimolecular reaction, compute all possible combinations for reactions along this edge.  */
		if (reac->order == 2) {
			
			/* A_1 + B_2 -> C_1 */
		   	a1  = evaluate_propensity(reac,xx,list_of_species,v1->vi,vol[v1->vi],v2->vi,vol[v2->vi],Mspecies);	
			cumsum += e->reaction_weight*a1;
			
			if (ur <= cumsum) {
				
				xx[v1->vi*Mspecies+reac->reactants[0]]+=reac->nr[reac->reactants[0]];
				xx[v2->vi*Mspecies+reac->reactants[1]]+=reac->nr[reac->reactants[1]];;
				xx[v1->vi*Mspecies+reac->products[0]]+=reac->nr[reac->products[0]];;
				break;
			}
			
			/* A_1 + B_2 -> C_2 */
		   	//a2  = evaluate_propensity(reac,xx,list_of_species,v1->vi,vol[v1->vi],v2->vi,vol[v2->vi],Mspecies);	
			cumsum+= e->reaction_weight*a1;
        
			if (ur <= cumsum) {
			
				xx[v1->vi*Mspecies+reac->reactants[0]]+=reac->nr[reac->reactants[0]];
				xx[v2->vi*Mspecies+reac->reactants[1]]+=reac->nr[reac->reactants[1]];
				xx[v2->vi*Mspecies+reac->products[0]]+=reac->nr[reac->products[0]];
				break;
			}
		}
	}
	
	if (r>=Mreactions) {
		printf("Error sampling non-local reaction\n.");
		r--;
	}

	return r;
	
	
}

int edge_execute_next_diffusion_event(edge *e,int *xx,int Mspecies,species **list_of_species,int dim){
/* Find the next diffusion event along an edge and update the state. */
	
	int s;
	double ur;
	int v1,v2;
	v1 = (e->v1)->vi;
	v2 = (e->v2)->vi;

	ur = (1.0-drand48())*e->total_diffusion_rate;
    
	double dij = e->weight;
    
    double cumsum=0.0;
   // double *delta_mu;
  //  delta_mu = calloc(Mspecies,sizeof(double));
    double delta_mu[Mspecies];
    edge_diffusion_excluded_volume_correction(e,xx,Mspecies,list_of_species,delta_mu,dim);
    
	for (s=0; s<Mspecies; s++) {
        cumsum+=list_of_species[s]->gamma*dij*xx[v1*Mspecies+s]*delta_mu[s];
		if (ur <= cumsum) {
			break;
		}
	}
	
	if (s>=Mspecies) {
		printf("Error sampling species to diffuse along the edge, total_diffusion_rate: %e ur: %e cumsum: %e.\n",e->total_diffusion_rate,ur,cumsum);
		s--;
	}
   	
	/* Species s diffuses from v1 to v2 */
	xx[v1*Mspecies+s]--;
	xx[v2*Mspecies+s]++;
	//free(delta_mu);
	return s;

}


void edge_diffusion_excluded_volume_correction(edge *e,int *xx,int Mspecies,species **list_of_species,double *delta_mu,int dim)
/* Compute the diffusion correction factor based on the chemical potential energy barrier to climb if a molecule of species spec is to jump from vertex 1 to vertex 2.  */
{
    vertex *v2;
    v2 = e->v2;
    int s;
    double gamma;

#ifdef HARD_SPHERES
    if (dim==2){
        double eta[3];
        compute_eta_2D(v2,xx,Mspecies,list_of_species,eta);
        for (s=0;s<Mspecies;s++){
            gamma = activity_coefficient_2D(v2,xx,Mspecies,list_of_species,s,eta);
            delta_mu[s] = 1.0/gamma;
        }
    }else{
        double eta[4];
        compute_eta_3D(v2,xx,Mspecies,list_of_species,eta);
        for (s=0;s<Mspecies;s++){
            gamma = activity_coefficient_3D(v2,xx,Mspecies,list_of_species,s,eta);
            delta_mu[s] = 1.0/gamma;
        }
    }
#else
    //delta_mu = 1.0;
#endif
    
}

double edge_compute_diffusion_rates(edge *e,int *xx,int Mspecies,species **list_of_species,int dim)
/* Compute the total diffusion rate along an edge. */
{
	
	double total_diffusion_rate=0.0;

	int node,s;
	double dij=e->weight;
	node = (e->v1)->vi;
    
   // double *delta_mu;
   // delta_mu = calloc(Mspecies,sizeof(double));
    double delta_mu[Mspecies];

    edge_diffusion_excluded_volume_correction(e,xx,Mspecies,list_of_species,delta_mu,dim);

	for (s=0; s<Mspecies; s++) {
        total_diffusion_rate+=dij*list_of_species[s]->gamma*xx[node*Mspecies+s]*delta_mu[s];
	}

    //free(delta_mu);
	return total_diffusion_rate;
	
}

int edge_execute_next_event(edge *e){
	
	int islocal = 0;
	return islocal;
}

double edge_update_rates(edge *e){
	return 0.0;
}

double sample_next_edge_event(edge *e){
	return 0.0;
}






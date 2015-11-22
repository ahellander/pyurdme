/* 
 
 Next Edge Method (NEM). Essentially the Next Subvolume Method (NSM) applied to edges since this is 
 more natural for interactions between adjacent neighbours. Bimolecular reactions occurs 
 along edges. 
 
*/

/* A. Hellander  2012-01-27  */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "species.h"
#include "reactions.h"
#include "vertex.h"
#include "edge.h"
#include "mesh.h"
#include "nem.h"
#include "binheap.h"


#include "hdf5.h"
#include "hdf5_hl.h"

//#define __DEBUG
#define pi 3.1416



void nemcore(fem_mesh *mesh,size_t *irD,size_t *jcD, double *prD,
			 const int *u0, const size_t *irN, const size_t *jcN, const int *prN,
			 const size_t *irG, const size_t *jcG, const double *tspan, const size_t tlen,
			 double *vol, const double *data, const int *sd, const size_t Ncells,
			 const size_t Mspecies, const size_t Mreactions, const size_t dsize,urdme_output_writer *writer)

{
	
	int node;
	int r,e,s;
	int Ndofs = Ncells*Mspecies;
	int number_of_vertices = Ncells;
	
	/*!!!!!!!!!!!!!!!!!*/
	int dim = 2;
	
	edge *next_edge,*adjacent_edge;
    vertex *vtx,*adjacent_node;
	
	/* Allocate the list of reactions and the list of species. */
	reaction  **list_of_reactions;
	list_of_reactions = initialize_reactions();	
	species **list_of_species;	
	list_of_species = initialize_species();	
	
	/* Initialize vertex data-structure. */
	vertex **list_of_vertices;
	list_of_vertices = initialize_vertices(number_of_vertices,Mspecies);
    for (node=0; node<number_of_vertices; node++) {
        for (s=0; s<Mspecies; s++) {
            list_of_vertices[node]->vol[s] = vol[node];
        }
    }
    
	/* Initialize the edge data-structure. ***/	
	double *p;
	p=mesh->p;
	
	size_t *jcK,*irK;
	double *prK;
	jcK=mesh->jcK;
	irK=mesh->irK;
	prK=mesh->prK;

	int number_of_edges;
	edge **list_of_edges;
	number_of_edges=(jcK[number_of_vertices] - number_of_vertices);
	
	list_of_edges=initialize_edges(number_of_edges,number_of_vertices,list_of_vertices,jcK,irK,prK,Mspecies);
   	
	/* Initialize the reaction weights. We will do this from an externally assembled weight-matrix later. */
	double vertex_weight=1.0;
    // To get the conventional RDME use neighbour_weight=0.0
	double neighbour_weight=0.0;
    
	/* There are four possible reaction combinations for each edge. */
	//double neighbour_weight = 0.5*(1.0-vertex_weight);
	//double neighbour_weight = (1.0-vertex_weight);

	
	edge_initialize_h(list_of_edges,number_of_edges,p,dim);
	double h;
	for (node=0; node<number_of_vertices; node++) {
        
       vtx=list_of_vertices[node];
	   h = sqrt(vol[0]);
       vtx->reaction_weight=vertex_weight;

	   for (e=0; e<vtx->number_of_outgoing_edges; e++) {
         list_of_edges[vtx->outgoing_edges[e]]->reaction_weight=neighbour_weight;
       }
	}
	

	/* Create a copy of the initial state matrix. */
	int *xx;
	xx=(int *)malloc(Ndofs*sizeof(int));
	memcpy(xx,u0,Ndofs*sizeof(int));
	
	double *node_rrates;
	node_rrates = (double *)calloc(number_of_vertices*Mreactions,sizeof(double));
		
	/* Compute all the local reaction rates (vertex rates). An edge will have these rates associated with its vertices. 
	   Several edges will have common node rates. */	
	double a0;
	for (node=0; node<number_of_vertices; node++)
	{
		a0=0.0;
		for (r=0; r<Mreactions; r++)
		{
		    node_rrates[node*Mreactions+r]=node_compute_reaction_rate(list_of_vertices[node],xx,list_of_reactions[r],list_of_species,vol,Mspecies);
			a0+=node_rrates[node*Mreactions+r];
		}
		list_of_vertices[node]->total_local_rate+=a0;
		list_of_vertices[node]->total_rate+=list_of_vertices[node]->total_local_rate;
	}
	
	/* Compute all the non-local reaction rates for all edges. */
	for (e=0; e<number_of_edges; e++) {
		a0=edge_compute_nonlocal_reaction_rates(list_of_edges[e],xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
		list_of_edges[e]->total_reaction_rate=a0;
		list_of_edges[e]->total_rate+=a0;
	}
	
	/* Compute the rates for the diffusion events for all edges. */
	for (e=0; e<number_of_edges; e++) {
		list_of_edges[e]->total_diffusion_rate=edge_compute_diffusion_rates(list_of_edges[e],xx,Mspecies,list_of_species,dim);
		list_of_edges[e]->total_rate+=list_of_edges[e]->total_diffusion_rate;
        //printf("Initilaization, total diffusion rate: %f\n",list_of_edges[e]->total_diffusion_rate);
	}
	
	/* Add total diffusion rates to vertices */
	for (node=0; node<number_of_vertices; node++) {
		vtx = list_of_vertices[node];
		for (e=0; e<vtx->number_of_outgoing_edges; e++) {
			adjacent_edge = list_of_edges[vtx->outgoing_edges[e]];
			vtx->total_edge_rate+=adjacent_edge->total_diffusion_rate;
		}
		vtx->total_rate+=vtx->total_edge_rate;
	}

	
	/* Initialize binary (min)heap. */
	int *index,*heap;
	index=(int *)malloc(number_of_vertices*sizeof(int));
	heap=(int *)malloc(number_of_vertices*sizeof(int));
	
	/* Sample the next time until an event for all edges. This time is the 
	   minimum of the times for the next local reaction, the next diffusion
	   event or the next non-local reaction. */
	double *taunode;
	taunode = (double *)malloc(number_of_vertices*sizeof(double));
	for (node=0; node<number_of_vertices; node++){
		taunode[node] = -log(1.0-drand48())/list_of_vertices[node]->total_rate +tspan[0];
		index[node]=node;
		heap[node]=node;
	}
	
	initialize_heap(taunode, index, heap, number_of_vertices);
    int num_diffusion_events=0;
	int it=0;
	double global_time = tspan[0];
	vertex *next_node;
	double ur,cumsum;
	double total_edge_rate;
	double total_local_rate;
	double total_rate;
	int rr;
	reaction *reac;
	double told;
    
    
	/* Sample events until done. */
	for (;; ) {
		
		/* Get the time until the next edge fires. */
		told = global_time;
		global_time = taunode[0];
	    //printf("global_time: %.4e\n",global_time);
		next_node = list_of_vertices[index[0]];
		//printf("Next node: %i\n",next_node->vi);
		
		/* Output the state at times given in tspan. */
		if (global_time>=tspan[it]||isinf(global_time)) {
			for (; it<tlen && (global_time>=tspan[it]||isinf(global_time)); it++){
                write_state(writer,xx);
            }
			if (it>=tlen){
                flush_buffer(writer);
                break;
            }
		}
        
		/* We know that the event occurs on an edge for which 
		   next_node is a member. Now find the edge that fires next. */
		
		/* The total node reaction rate (local contribution) */
		total_rate       = next_node->total_rate;
        total_local_rate = next_node->total_local_rate; 
		total_edge_rate  = next_node->total_edge_rate;
		
        if ((1.0-drand48())*total_rate<=total_local_rate) {
            /* A local reaction. */
            
			/* Find the local vertex reaction event that occurs next. */
			ur = drand48()*total_local_rate;
			for (r=0, cumsum=node_rrates[next_node->vi*Mreactions]; r < Mreactions && ur>cumsum;
				 r++, cumsum+=node_rrates[next_node->vi*Mreactions+r]);
            
			rr = r;
			/* Execute the reaction. */	
			execute_reaction(list_of_reactions[rr],xx,next_node->vi,Mspecies);

			/* Update all the affected local rates. */
			list_of_vertices[next_node->vi]->total_rate-=list_of_vertices[next_node->vi]->total_local_rate;	
			list_of_vertices[next_node->vi]->total_local_rate=0.0;
			a0=0.0;
			for (r=0; r<Mreactions; r++)
			{
				node_rrates[next_node->vi*Mreactions+r]=node_compute_reaction_rate(next_node,xx,list_of_reactions[r],list_of_species,vol,Mspecies);
				a0+=node_rrates[next_node->vi*Mreactions+r];
			}
			list_of_vertices[next_node->vi]->total_local_rate=a0;
			list_of_vertices[next_node->vi]->total_rate+=list_of_vertices[next_node->vi]->total_local_rate;	
			
			/* Update the total diffusion rates for all outgoing edges */
            reac=list_of_reactions[rr];
			
			next_node->total_rate-=next_node->total_edge_rate;
			for (e=0; e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge=list_of_edges[next_node->outgoing_edges[e]];
				
	            next_node->total_edge_rate-=adjacent_edge->total_rate;
				
				adjacent_edge->total_rate-=adjacent_edge->total_diffusion_rate;
				adjacent_edge->total_diffusion_rate=edge_compute_diffusion_rates(adjacent_edge,xx,Mspecies,list_of_species,dim);
				adjacent_edge->total_rate+=adjacent_edge->total_diffusion_rate;
				
				next_node->total_edge_rate+=adjacent_edge->total_rate;	
			}
			next_node->total_rate+=next_node->total_edge_rate;
            
			/* Update non-local rates for incoming and outgoing reactions */
			
			/* Update all outgoing affected edges from the current node. */
			next_node->total_rate-=next_node->total_edge_rate;
			for (e=0;e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->outgoing_edges[e]];
				next_node->total_edge_rate-=adjacent_edge->total_rate;
				
				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0=edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;
				
				next_node->total_edge_rate+=adjacent_edge->total_rate;
			}
			next_node->total_rate+=next_node->total_edge_rate;
			
			/* Update all incoming edges from neighbouring nodes to next_node. */
			for (e=0; e<next_node->number_of_incoming_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->incoming_edges[e]];
				adjacent_node = adjacent_edge->v1;
				
				adjacent_node->total_rate-=adjacent_node->total_edge_rate;
				adjacent_node->total_edge_rate-=adjacent_edge->total_rate;
				
				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0 = edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;
                
                /* Need to update the diffusion rates of the incoming edges if we have excluded volume effects.*/
                adjacent_edge->total_rate-=adjacent_edge->total_diffusion_rate;
				adjacent_edge->total_diffusion_rate=edge_compute_diffusion_rates(adjacent_edge,xx,Mspecies,list_of_species,dim);
				adjacent_edge->total_rate+=adjacent_edge->total_diffusion_rate;

				
				adjacent_node->total_edge_rate+=adjacent_edge->total_rate;
				adjacent_node->total_rate+=adjacent_node->total_edge_rate;
				
			}
			
			/* Update the heap for next_node affected nodes */
			if (next_node->total_rate>0.0) {	
				taunode[heap[next_node->vi]] = -log(1.0-drand48())/next_node->total_rate + global_time;
			}else{
				taunode[heap[next_node->vi]] = INFINITY;
				
			}
			update(heap[next_node->vi], taunode, index, heap, number_of_vertices);
			
			/* Update heap */
			for (e=0;e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->outgoing_edges[e]];
				adjacent_node = adjacent_edge->v2;
				
				if (adjacent_node->total_rate>0.0) {	
					taunode[heap[adjacent_node->vi]] = -log(1.0-drand48())/adjacent_node->total_rate + global_time;
				}else{
					taunode[heap[adjacent_node->vi]] = INFINITY;
					
				}
				update(heap[adjacent_node->vi], taunode, index, heap, number_of_vertices);
				
			}
			
					
			
		} 
		else {
			

			/* Find the egde on which the next event occurs. */
			ur = (1.0-drand48())*next_node->total_edge_rate;
			cumsum =0.0;
			for (e=0; e<next_node->number_of_outgoing_edges; e++) {
				adjacent_edge = list_of_edges[next_node->outgoing_edges[e]];
				cumsum += adjacent_edge->total_rate;
				if (ur <= cumsum) {
					next_edge = adjacent_edge;
					break;
				}
				
			}
			
			if (e>=next_node->number_of_outgoing_edges){
				printf("Error sampling next edge.\n");
				e--;
				next_edge=list_of_edges[next_node->outgoing_edges[e]];
			}
			
			/* Now find the type of event that occured on that edge (reaction or diffusion). */
			if ((1.0-drand48())*next_edge->total_rate <= next_edge->total_diffusion_rate) {
                num_diffusion_events++;
//printf("A diffusion event on edge %i %.4e %.4e %.4e.\n",next_edge->index,next_edge->total_rate,next_edge->total_diffusion_rate,next_edge->total_reaction_rate);
				/* Find the species that diffuses along the next edge, and update the state. */
				s=edge_execute_next_diffusion_event(next_edge,xx,Mspecies,list_of_species,dim);
								
			}
			else {
				
				/* Find the non-local reaction event that occured along the edge. */
				//printf("A non-local reaction event at edge %i. %.4e %.4e %.4e\n",next_edge->index,next_edge->total_rate,next_edge->total_diffusion_rate,next_edge->total_reaction_rate);
				//exit(-1);
				r = edge_execute_next_reaction_event(next_edge,xx,Mreactions,list_of_reactions,Mspecies,list_of_species,vol);
			
			}
			
			/* Update all affected rates (for vertices and edges) */
			
			/* Adjust all affected local reaction rates. */
			list_of_vertices[next_node->vi]->total_rate-=list_of_vertices[next_node->vi]->total_local_rate;	
			list_of_vertices[next_node->vi]->total_local_rate=0.0;
			a0=0.0;
			for (r=0; r<Mreactions; r++)
			{
				node_rrates[next_node->vi*Mreactions+r]=node_compute_reaction_rate(next_node,xx,list_of_reactions[r],list_of_species,vol,Mspecies);
				a0+=node_rrates[next_node->vi*Mreactions+r];
			}
			list_of_vertices[next_node->vi]->total_local_rate+=a0;
			list_of_vertices[next_node->vi]->total_rate+=list_of_vertices[next_node->vi]->total_local_rate;	
			
			vtx = next_edge->v2;
			list_of_vertices[vtx->vi]->total_rate-=list_of_vertices[vtx->vi]->total_local_rate;	
			list_of_vertices[vtx->vi]->total_local_rate=0.0;
			a0=0.0;
			for (r=0; r<Mreactions; r++)
			{
				node_rrates[vtx->vi*Mreactions+r]=node_compute_reaction_rate(vtx,xx,list_of_reactions[r],list_of_species,vol,Mspecies);
				a0+=node_rrates[vtx->vi*Mreactions+r];
			}
			list_of_vertices[vtx->vi]->total_local_rate+=a0;
			list_of_vertices[vtx->vi]->total_rate+=list_of_vertices[vtx->vi]->total_local_rate;	
			
			
			/* Adjust all the affected diffusion edge-rates and the total vertex rate. */
			next_node->total_rate-=next_node->total_edge_rate;
			for (e=0; e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge=list_of_edges[next_node->outgoing_edges[e]];
				
	            next_node->total_edge_rate-=adjacent_edge->total_rate;
				
				adjacent_edge->total_rate-=adjacent_edge->total_diffusion_rate;
				adjacent_edge->total_diffusion_rate=edge_compute_diffusion_rates(adjacent_edge,xx,Mspecies,list_of_species,dim);
				adjacent_edge->total_rate+=adjacent_edge->total_diffusion_rate;
		
				next_node->total_edge_rate+=adjacent_edge->total_rate;	
			}
			next_node->total_rate+=next_node->total_edge_rate;
            
			/* Adjust all the affected diffusion edge-rates and the total vertex rate. */
			vtx=next_edge->v2;
			vtx->total_rate-=vtx->total_edge_rate;
			for (e=0; e<vtx->number_of_outgoing_edges; e++) {
				
				adjacent_edge=list_of_edges[vtx->outgoing_edges[e]];
			
	            vtx->total_edge_rate-=adjacent_edge->total_rate;
				adjacent_edge->total_rate-=adjacent_edge->total_diffusion_rate;
				adjacent_edge->total_diffusion_rate=edge_compute_diffusion_rates(adjacent_edge,xx,Mspecies,list_of_species,dim);
				adjacent_edge->total_rate+=adjacent_edge->total_diffusion_rate;
				
				vtx->total_edge_rate+=adjacent_edge->total_rate;	
			}
			vtx->total_rate+=vtx->total_edge_rate;
            
			/* Update all outgoing affected edges from the current node. */
			next_node->total_rate-=next_node->total_edge_rate;
			for (e=0;e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->outgoing_edges[e]];
				next_node->total_edge_rate-=adjacent_edge->total_rate;

				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0=edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;
				
				next_node->total_edge_rate+=adjacent_edge->total_rate;
			}
			next_node->total_rate+=next_node->total_edge_rate;
			
			/* Update all incoming edges from neighbouring nodes to next_node. */
			for (e=0; e<next_node->number_of_incoming_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->incoming_edges[e]];
				adjacent_node = adjacent_edge->v1;
				
				adjacent_node->total_rate-=adjacent_node->total_edge_rate;
				adjacent_node->total_edge_rate-=adjacent_edge->total_rate;
			
				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0 = edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;
				
				adjacent_node->total_edge_rate+=adjacent_edge->total_rate;
				adjacent_node->total_rate+=adjacent_node->total_edge_rate;
				
			}
			
			/* Update all outgoing affected edges from the second node. */
			vtx = next_edge->v2;
			vtx->total_rate-=vtx->total_edge_rate;
			for (e=0;e<vtx->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[vtx->outgoing_edges[e]];
				vtx->total_edge_rate-=adjacent_edge->total_rate;
				
				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0=edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;

				vtx->total_edge_rate+=adjacent_edge->total_rate;
			}
			vtx->total_rate+=vtx->total_edge_rate;
			
			/* Update all incoming edges from neighbouring nodes to the second node. */
			for (e=0; e<vtx->number_of_incoming_edges; e++) {
				
				adjacent_edge = list_of_edges[vtx->incoming_edges[e]];
				adjacent_node = adjacent_edge->v1;
				
				adjacent_node->total_rate-=adjacent_node->total_edge_rate;
				adjacent_node->total_edge_rate-=adjacent_edge->total_rate;

				adjacent_edge->total_rate-=adjacent_edge->total_reaction_rate;
				a0 = edge_compute_nonlocal_reaction_rates(adjacent_edge,xx,list_of_reactions,list_of_species,vol,Mspecies,Mreactions);
				adjacent_edge->total_reaction_rate=a0;
				adjacent_edge->total_rate+=adjacent_edge->total_reaction_rate;
                
                /* Need to update diffusion rates if we have excluded volume effects */
                adjacent_edge->total_rate-=adjacent_edge->total_diffusion_rate;
				adjacent_edge->total_diffusion_rate=edge_compute_diffusion_rates(adjacent_edge,xx,Mspecies,list_of_species,dim);
				adjacent_edge->total_rate+=adjacent_edge->total_diffusion_rate;


				adjacent_node->total_edge_rate+=adjacent_edge->total_rate; 
				adjacent_node->total_rate+=adjacent_node->total_edge_rate;
				
			}
			
			
			/* Update heap */
			for (e=0;e<next_node->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[next_node->outgoing_edges[e]];
				adjacent_node = adjacent_edge->v2;
				
				if (adjacent_node->total_rate>0.0) {	
					taunode[heap[adjacent_node->vi]] = -log(1.0-drand48())/adjacent_node->total_rate + global_time;
				}else{
					taunode[heap[adjacent_node->vi]] = INFINITY;
					
				}
				update(heap[adjacent_node->vi], taunode, index, heap, number_of_vertices);
				
			}
			
			for (e=0;e<vtx->number_of_outgoing_edges; e++) {
				
				adjacent_edge = list_of_edges[vtx->outgoing_edges[e]];
				adjacent_node = adjacent_edge->v2;
				
				if (adjacent_node->total_rate>0.0) {	
					taunode[heap[adjacent_node->vi]] = -log(1.0-drand48())/adjacent_node->total_rate + global_time;
				}else{
					taunode[heap[adjacent_node->vi]] = INFINITY;
					
				}
				update(heap[adjacent_node->vi], taunode, index, heap, number_of_vertices);
				
			}
			
		}
		
		/* Some error checking. */
#if defined(DEBUG)
		for (int j=0; j<number_of_vertices; j++) {
			vtx = list_of_vertices[j];
			if (fabs(vtx->total_rate - (vtx->total_local_rate+vtx->total_edge_rate))>1e-10) {
				printf("Error vtx: %.4e %.4e %.4e %.4e\n",vtx->total_rate,vtx->total_local_rate,vtx->total_edge_rate,vtx->total_local_rate+vtx->total_edge_rate);
				exit(-1);
			}
		}
		
		for (int j=0; j<number_of_edges; j++) {
			next_edge = list_of_edges[j];
			if (fabs(next_edge->total_rate - (next_edge->total_diffusion_rate+next_edge->total_reaction_rate))>1e-10) {
				printf("Error edge: %.4e %.4e %.4e %.4e\n",next_edge->total_rate,next_edge->total_diffusion_rate,next_edge->total_reaction_rate,next_edge->total_diffusion_rate+next_edge->total_reaction_rate);
				exit(-1);
			}
			
		}
		
		for (int j=0; j<Mspecies*number_of_vertices; j++) {
			if (xx[j]<0) {
				printf("Negative copy number.\n");
				exit(-1);
			}
		}
#endif
		
			
	}
  
	
	free_list_of_vertices(list_of_vertices,number_of_vertices);
	free_list_of_edges(list_of_edges,number_of_edges);
	free(list_of_species);
	free(list_of_reactions);
	free(node_rrates);
	free(heap);
	free(index);
	free(taunode);
  	
    printf("Num_diffusion_events: %i\\",num_diffusion_events);
    
}


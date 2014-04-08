/* nsmcore.c - Core NSM solver. Generates one trajectory of the process.  */

/* A. Hellander 2012-06-15 (Revision) */
/* P. Bauer and S. Engblom 2012-04-10 (Revision) */
/* A. Hellander 2009-11-24 (Revision) */
/* J. Cullhed 2008-06-18 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "propensities.h"
#include "nsm.h"
#include "binheap.h"
#include "report.h"

#include "hdf5.h"
#include "hdf5_hl.h"


int flush_solution_to_file(hid_t trajectory_dataset,int *buffer,int column_offset, int num_columns, int col_size)
{

    /* This  is the column offset in the hdf5 datafile. */
    hsize_t start[2];
    hsize_t count[2];
    hsize_t block[2];
    
    /* Some parameters for the hyperslabs we need to select. */
    start[0] = column_offset;
    start[1] = 0;
    
    count[0] = 1;
    count[1] = 1;
    
    block[0] = num_columns;
    block[1] = col_size;
    
    herr_t status;
    
    /* A memory space to indicate the size of the buffer. */
    hsize_t mem_dims[2];
    mem_dims[0] = num_columns;
    mem_dims[1] = col_size;
    hid_t mem_space = H5Screate_simple(2, mem_dims, NULL);
    
    hid_t file_dataspace = H5Dget_space(trajectory_dataset);
    //  if (file_dataspace == NULL){
    //      printf("Failed to get the dataspace of the dataset.");
    //      return(-1);
    //  }
    status =H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, start, NULL,count,block);
    if (status){
        printf("Failed to select hyperslab.");
        return(-1);
    }
    status = H5Dwrite(trajectory_dataset, H5T_NATIVE_INT, mem_space,file_dataspace,H5P_DEFAULT,buffer);
    if (status){
        printf("Failed to write to the dataset.");
        return(-1);
    }
    
    H5Sclose(mem_space);
    H5Sclose(file_dataspace);
    return(0);
}

void nsm_core(const size_t *irD,const size_t *jcD,const double *prD,
              const int *u0,
              const size_t *irN,const size_t *jcN,const int *prN,
              const size_t *irG,const size_t *jcG,
              const double *tspan,const size_t tlen,
              const double *vol,const double *data,const int *sd,
              const size_t Ncells,
              const size_t Mspecies,const size_t Mreactions,
              const size_t dsize,int report_level, hid_t output_file,
              const size_t *irK,const size_t *jcK,const double *prK)

/* Specification of the inputs:
 
 Ncells
 Number of subvolumes.
 
 Mspecies
 Number of species.
 
 Hence Ndofs = Ncells*Mspecies.
 
 Mreactions
 Total number of reactions.
 
 dsize
 Size of data vector sent to propensities.
 
 tlen
 Number of sampling points in time.
 
 report_level
 The desired degree of feedback during simulations. 0, 1, and 2 are
 currently supported options.
 
 Diffusion matrix D. Double sparse (Ndofs X Ndofs).
 Macroscopic diffusion matrix. D(i,j) is the diffusion rate from dof #j to
 dof #i. This matrix uses the CSR-format and not CSC because fast access to
 rows is needed.
 
 Initial state vector u0. Integer (Mspecies X Ncells).
 Gives the initial copy number of the species in each subvolume.
 
 Stochiometric matrix N. Integer sparse (Mspecies X Nreactions).
 N(:,j) describes how reaction j changes the number of species.
 
 Dependency graph G. Integer sparse (Mreactions X Mspecies+Mreactions).
 G(i,Mspecies+j) is non-zero if executing reaction j means that reaction i
 needs to be re-evaluated. The first Mspecies columns of G similarily cover
 diffusion events.
 
 tspan. Double vector.
 Output times. tspan[0] is the start time and tspan[length(tspan)-1] is the
 stop time.
 
 vol. Double vector (length Ncells).
 vol[i] gives the volume of cell #i.
 
 data. Double matrix (dsize X Ncells).
 Generalized data matrix, data(:,j) gives a data vector for cell #j.
 
 sd. Integer vector (length Ncells).
 Subdomain number. sd[i] is the subdomain of cell #i. The vector sd can also
 be used to separate boundaries, line segments and points.
 
 Format of sparse matrices:
 G, N and S are sparse matrices in compressed column format (CCS). D is sparse
 but in compressed row format (CRS), or equivalently, a transposed matrix in
 CCS format.
 jcD, irD, prD (double *)
 jcN, irN, prN (int *)
 jcG, irG (int *)
 
 Propensities:
 a vector of function pointers (length Mreactions) is input by
 linking with the prototypes in propensities.h and function
 definitions in a user-specified .c-file. The type of this vector is
 PropensityFun which defines the input to a property function. See
 propensities.h for more details.
 
 Ordering of the dofs:
 Dof #i is located in cell #(i/Mspecies), and the dofs located in
 cell #j is u0(:,j). Thus, u0 is understood as a matrix of size
 Mspecies X Ncells.
 
 The output is a matrix U (Ndofs X length(tspan)).
 U(:,j) contains the state of the system at tspan(j).
 */
{
    double told,tt = tspan[0];
    double rdelta,rrdelta;
    double rand,cum,old;
    double *srrate,*rrate;
    double *sdrate,*Ddiag;
    double *rtimes;
    double old_rrate = 0.0,old_drate = 0.0;
    
    double totrate;
    
    int *node,*heap,*xx;
    long int total_reactions = 0;
    long int total_diffusion = 0;
    int dof,col,s;
    
    int subvol,event,re,spec,errcode = 0;
    size_t i,j,it = 0;
    size_t to_node,to_vol = 0;
    const size_t Ndofs = Ncells*Mspecies;
    
    
    PropensityFun *rfun;
    rfun = ALLOC_propensities();
    
    ReportFun report;
    if (report_level)
    report = &reportFun1;
    else
    report = NULL;
    
    
    /* Add some metadata to the output file and create datasets */
    // Write tspan to the file
    herr_t status;
    hsize_t dataset_dims[2]; /* dataset dimensions */
    hsize_t chunk_dims[2];
    
    dataset_dims[0] = 1;
    dataset_dims[1] = tlen;
    status = H5LTmake_dataset(output_file,"/tspan",2,dataset_dims,H5T_NATIVE_DOUBLE,tspan);
    if (status != 0){
        printf("Failed to write tspan vector HDF5 file.");
        exit(-1);
    }
    
    hid_t trajectory_dataset, datatype,trajectory_dataspace;
    datatype = H5Tcopy(H5T_NATIVE_INT);
    
    /* This is the maximal buffer size we use to store the solution before writing to file. */
    size_t max_buffer_size = 1048576*8;
    
    /* How many timepoints do we log before the buffer is full? */
    size_t column_size = Ndofs*sizeof(int);
    int num_columns = max_buffer_size / column_size;
    if (num_columns > tlen){
        num_columns = tlen;
    }
    size_t buffer_size = num_columns*Ndofs;
    int *buffer = (int *)calloc(buffer_size,sizeof(int));
    int num_columns_since_flush = 0;
    int chunk_indx=0;
    
    dataset_dims[0] = tlen;
    dataset_dims[1] = Ndofs;
    
    chunk_dims[0] = num_columns;
    chunk_dims[1] = Ndofs;
    
    trajectory_dataspace = H5Screate_simple(2, dataset_dims, NULL);
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    //H5Pset_chunk(plist,2,dataset_dims);
    H5Pset_chunk(plist,2,chunk_dims);
    
    //status = H5Pset_deflate (plist, 4);
    
    trajectory_dataset = H5Dcreate2(output_file, "/U", datatype, trajectory_dataspace, H5P_DEFAULT,plist,H5P_DEFAULT);
    
    /* Set xx to the initial state. xx will always hold the current solution. */
    xx = (int *)malloc(Ndofs*sizeof(int));
    memcpy(xx,u0,Ndofs*sizeof(int));
    
    /* Create reaction rate matrix (Mreactions X Ncells) and total rate
     vector. In rrate we store all propensities for chemical rections,
     and in srrate the sum of propensities in every subvolume. */
    rrate = (double *)malloc(Mreactions*Ncells*sizeof(double));
    srrate = (double *)malloc(Ncells*sizeof(double));
    
    
    /* Calculate the propensity for every reaction and every
     subvolume. Store the sum of the reaction intensities in each
     subvolume in srrate. */
    for (i = 0; i < Ncells; i++) {
        srrate[i] = 0.0;
        for (j = 0; j < Mreactions; j++) {
            //rrate[i*Mreactions+j] =
            //(*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i],i,xx,irK,jcK,prK);
            //srrate[i] += rrate[i*Mreactions+j];
            rrate[i*Mreactions+j] =
            (*rfun[j])(&xx[i*Mspecies],tt,vol[i],&data[i*dsize],sd[i]);
            srrate[i] += rrate[i*Mreactions+j];
        }
    }
    
    /* Total diffusion rate vector (length Mcells). It will hold
     the total diffusion rates in each subvolume. */
    sdrate = (double *)malloc(Ncells*sizeof(double));
    
    /* The diagonal value of the D-matrix is used frequently. For
     efficiency, we store the negative of D's diagonal in Ddiag. */
    Ddiag = (double *)malloc(Ndofs*sizeof(double));
    for (i = 0; i < Ndofs; i++) {
        Ddiag[i] = 0.0;
        for (j = jcD[i]; j < jcD[i+1]; j++)
        if (irD[j] == i) Ddiag[i] = -prD[j];
    }
    
    /* Calculate the total diffusion rate for each subvolume. */
    for(i = 0; i < Ncells; i++) {
        sdrate[i] = 0.0;
        for(j = 0; j < Mspecies; j++)
        sdrate[i] += Ddiag[i*Mspecies+j]*xx[i*Mspecies+j];
    }
    
    /* Create binary (min)heap. */
    rtimes = (double *)malloc(Ncells*sizeof(double));
    node = (int *)malloc(Ncells*sizeof(int));
    heap = (int *)malloc(Ncells*sizeof(int));
    
    /* Calculate times to next event (reaction or diffusion)
     in each subvolume and initialize heap. */
    for (i = 0; i < Ncells; i++) {
        rtimes[i] = -log(1.0-drand48())/(srrate[i]+sdrate[i])+tspan[0];
        heap[i] = node[i] = i;
    }
    initialize_heap(rtimes,node,heap,Ncells);
    int total_columns_written = 0;
    /* Main loop. */
    for (;;) {
        /* Get the subvolume in which the next event occurred.
         This subvolume is on top of the heap. */
        told = tt;
        tt   = rtimes[0];
        subvol = node[0];
        
        /* Store solution if the global time counter tt has passed the
         next time is tspan. */
        if (tt >= tspan[it] || isinf(tt)) {
	  
            for (; it < tlen && (tt >= tspan[it] || isinf(tt)); it++) {
                
                if (report)
                    report(tspan[it],tspan[0],tspan[tlen-1],total_diffusion,total_reactions,0,report_level);
                
                
                /* Write to the buffer */
                for (i=0; i<Ncells;i++){
                    for (s=0; s<Mspecies; s++){
                        buffer[Ndofs*num_columns_since_flush+s*Ncells+i] = xx[i*Mspecies+s];
                    }
                }
                
                
                num_columns_since_flush++;

                if (num_columns_since_flush == num_columns){

                    status = flush_solution_to_file(trajectory_dataset,buffer,chunk_indx*num_columns,num_columns_since_flush,Ndofs);
                    if (status){
                        printf("Failed to flush buffer to file.");
                        exit(-1);
                    }
                    
                    num_columns_since_flush = 0;
                    chunk_indx++;
                    total_columns_written += num_columns;
                }
                
            }
            
            /* If the simulation has reached the final time, exit. */
            if (it >= tlen){
                /* Flush the buffer */
		
		if (num_columns_since_flush > 0){
		  status = flush_solution_to_file(trajectory_dataset,buffer,chunk_indx*num_columns,num_columns_since_flush,Ndofs);
		  if (status){
                    printf("Failed to flush buffer to file.\n");
                    exit(-1);
		  }
		}
                total_columns_written += num_columns_since_flush;
                break;
            }
            
        }
        
        /* First check if it is a reaction or a diffusion event. */
        totrate = srrate[subvol]+sdrate[subvol];
        rand = drand48();
        
        
        if (rand*totrate <= srrate[subvol]) {
            /* Reaction event. */
            event = 0;
            
            /* a) Determine the reaction re that did occur (direct SSA). */
            rand *= totrate;
            for (re = 0, cum = rrate[subvol*Mreactions]; re < Mreactions && rand > cum; re++, cum += rrate[subvol*Mreactions+re])
            ;
            
            if (re >= Mreactions){
                re--;
            }
            
            /* b) Update the state of the subvolume subvol and sdrate[subvol]. */
            for (i = jcN[re]; i < jcN[re+1]; i++) {
                xx[subvol*Mspecies+irN[i]] += prN[i];
                if (xx[subvol*Mspecies+irN[i]] < 0) errcode = 1;
                sdrate[subvol] += Ddiag[subvol*Mspecies+irN[i]]*prN[i];
            }
            
            /* c) Recalculate srrate[subvol] using dependency graph. */
            for (i = jcG[Mspecies+re], rdelta = 0.0; i < jcG[Mspecies+re+1]; i++) {
                old = rrate[subvol*Mreactions+irG[i]];
                j = irG[i];
                //rdelta +=
                //(rrate[subvol*Mreactions+j] =
                // (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol],subvol,xx,irK,jcK,prK)
                // )-old;
                rdelta +=
                (rrate[subvol*Mreactions+j] =
                 (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol])
                 )-old;
            }
            srrate[subvol] += rdelta;
            
            total_reactions++; /* counter */
        }
        else {
            /* Diffusion event. */
            event = 1;
            
            /* a) Determine which species... */
            rand *= totrate;
            rand -= srrate[subvol];
            
            for (spec = 0, dof = subvol*Mspecies, cum = Ddiag[dof]*xx[dof];
                 spec < Mspecies && rand > cum;
                 spec++, cum += Ddiag[dof+spec]*xx[dof+spec]);
            
            
            /* b) and then the direction of diffusion. */
            col = dof+spec;
            rand = drand48()*Ddiag[col];
            
            /* Search for diffusion direction. */
            for (i = jcD[col], cum = 0.0; i < jcD[col+1]; i++)
            if (irD[i] != col && (cum += prD[i]) > rand)
            break;
            
            /* paranoia fix: */
            if (i >= jcD[col+1]){
               i--;
            }
            
            to_node = irD[i];
            to_vol = to_node/Mspecies;
            
            /* c) Execute the diffusion event (check for negative elements). */
            xx[subvol*Mspecies+spec]--;
            if (xx[subvol*Mspecies+spec] < 0){
                    errcode = 1;
            }
            xx[to_node]++;
            
            /* Save reaction and diffusion rates. */
            old_rrate = srrate[to_vol];
            old_drate = sdrate[to_vol];
            
            /* Recalculate the reaction rates using dependency graph G. */
            if (Mreactions > 0){
                for (i = jcG[spec], rdelta = 0.0, rrdelta = 0.0; i < jcG[spec+1]; i++) {
                    
                    j = irG[i];
                    old = rrate[subvol*Mreactions+j];
                    
                    //rdelta +=
                    //  (rrate[subvol*Mreactions+j] =
                    //    (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol],subvol,xx,irK,jcK,prK)
                    //  )-old;
                    rdelta +=
                      (rrate[subvol*Mreactions+j] =
                        (*rfun[j])(&xx[subvol*Mspecies],tt,vol[subvol],&data[subvol*dsize],sd[subvol])
                      )-old;
                    old = rrate[to_vol*Mreactions+j];
                    
                    //rrdelta += 
					//  (rrate[to_vol*Mreactions+j] = 
                    //    (*rfun[j])(&xx[to_vol*Mspecies],tt,vol[to_vol],&data[to_vol*dsize],sd[to_vol],to_vol,xx,irK,jcK,prK)
                    //  )-old;
                    rrdelta += 
					  (rrate[to_vol*Mreactions+j] = 
                        (*rfun[j])(&xx[to_vol*Mspecies],tt,vol[to_vol],&data[to_vol*dsize],sd[to_vol])
                      )-old;
                }
                
                srrate[subvol] += rdelta;
                srrate[to_vol] += rrdelta;
            }
            
            /* Adjust diffusion rates. */
            sdrate[subvol] -= Ddiag[subvol*Mspecies+spec];
            sdrate[to_vol] += Ddiag[to_vol*Mspecies+spec];
            
            total_diffusion++; /* counter */
            
        }
        
        /* Compute time to new event for this subvolume. */
        totrate = srrate[subvol]+sdrate[subvol];  
        if (totrate > 0.0)
        rtimes[0] = -log(1.0-drand48())/totrate+tt;
        else
        rtimes[0] = INFINITY;
        
        /* Update the heap. */
        update(0,rtimes,node,heap,Ncells);
        
        /* If it was a diffusion event, also update the other affected
         node. */
        if (event) {
            totrate = srrate[to_vol]+sdrate[to_vol];      
            if (totrate > 0.0) {
                if (!isinf(rtimes[heap[to_vol]]))
                rtimes[heap[to_vol]] = 
                (old_rrate+old_drate)/totrate*(rtimes[heap[to_vol]]-tt)+tt;
                else
                /* generate a new waiting time */
                rtimes[heap[to_vol]] = -log(1.0-drand48())/totrate+tt;
            } 
            else
            rtimes[heap[to_vol]] = INFINITY;
            
            update(heap[to_vol],rtimes,node,heap,Ncells);
        } 
        
        /* Check for error codes. */
        if (errcode) {
            /* Report the error that occurred. */
            if (report)
            report(tt,tspan[0],tspan[tlen-1],total_diffusion,total_reactions,errcode,report_level);
            /* Cannot continue. Clear this solution and exit. */
            break;
        }
    }
    
    /* Close the datasets */
    H5Sclose(trajectory_dataspace);
    H5Dclose(trajectory_dataset);
    
    FREE_propensities(rfun);
    free(buffer);
    free(heap);
    free(node);
    free(rtimes);
    free(Ddiag);
    free(sdrate);
    free(srrate);
    free(rrate);
    free(xx); 
    
}

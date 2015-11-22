/* nsm.h */

/* J. Cullhed    2008-06-18. 
   A. Hellander  2010-01-16 */

#ifndef RDME__H
#define RDME__H

#include "mesh.h"
#include "outputwriter.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void *nem(void *data, urdme_output_writer *);
void nemcore(fem_mesh *mesh,size_t *irD,size_t *jcD,double *prD,
		 const int *u0, const size_t *irN, const size_t *jcN, const int *prN,
		 const size_t *irG, const size_t *jcG, const double *tspan, const size_t tlen,
		 double *vol, const double *data, const int *sd, const size_t Ncells,
		 const size_t Mspecies, const size_t Mreactions, const size_t dsize, urdme_output_writer *writer);




#endif /* RDME__H */


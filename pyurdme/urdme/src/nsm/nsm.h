/* nsm.h - Header file for NSM-solver. */

/* P. Bauer and S. Engblom 2012-05-04 (Revision) */
/* A. Hellander 2009-11-24 (Revision) */
/* J. Cullhed 2008-06-18 */

#ifndef __nsm_h
#define __nsm_h

#include "hdf5.h"
#include "hdf5_hl.h"

void nsm(void *data, hid_t output_file);
void nsm_core(const size_t *irD,const size_t *jcD,const double *prD,
              const int *u0,
              const size_t *irN,const size_t *jcN,const int *prN,
              const size_t *irG,const size_t *jcG,
              const double *tspan,const size_t tlen,
              const double *vol,const double *data,
              const int *sd,const size_t Ncells,
              const size_t Mspecies,const size_t Mreactions,
              const size_t dsize,
              int report_level, hid_t output_file);

#endif /* __nsm_h */

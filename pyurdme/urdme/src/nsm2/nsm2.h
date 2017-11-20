/* nsm2.h - Header file for NSM-solver. */

/* P. Bauer and S. Engblom 2012-05-04 (Revision) */
/* A. Hellander 2009-11-24 (Revision) */
/* J. Cullhed 2008-06-18 */

#ifndef __nsm2_h
#define __nsm2_h

#include "hdf5.h"
#include "hdf5_hl.h"
#include "urdmemodel.h"
#include "outputwriter.h"

void nsm2(void *data, urdme_output_writer *writer);
double rfun(const int *x, double t, const double vol, const double *data, int sd,
            const int ireaction,
            const double *R, const size_t Mreactions,
            const int *I, const int *S, const size_t Msubdomains);
void nsm2_core(const size_t *irD,const size_t *jcD,const double *prD,
              const int *u0,
              const size_t *irN,const size_t *jcN,const int *prN,
              const size_t *irG,const size_t *jcG,
              const double *tspan,const size_t tlen,
              const double *vol,const double *data,const int *sd,
              const size_t Ncells,
              const size_t Mspecies,const size_t Mreactions,
              const size_t Msubdomains,
              const size_t dsize,int report_level,
			  const size_t *irK,const size_t *jcK,const double *prK, 
              const double *R, const int *I, const int *S,
              urdme_output_writer *writer);

#endif /* __nsm2_h */

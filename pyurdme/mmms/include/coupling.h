/* Utility routines for micro/meso coupling. */

#ifndef COUPLING__H
#define COUPLING__H
//#include "micro.h"
#include <vector>
//#include "micro.h"
//#include "matmodel.h"
#include "mesh.h"
#include "urdmemodel.h"
#include "structs.h"

#define micro2meso micro2meso3
#define meso2micro meso2micro2

int micro2meso1(particle *, fem_mesh *);
int micro2meso2(particle *, fem_mesh *, int Mspecies);
int micro2meso3(particle *particle, fem_mesh *mesh, int Mspecies);
int micro2meso_fp(particle *particle, fem_mesh *mesh);


void meso2micro(particle *,int voxel, int species, double rr, double d, int dim,fem_mesh* mesh);
void meso2micro2(particle *out,int voxel, int species, double rr, double d,int dim,fem_mesh* mesh);
void meso2micro1(particle *out,int voxel, int species, double rr, double d,int dim,fem_mesh* mesh);

vector <plane> voxel_boundaries(urdme_model *model, fem_mesh *mesh);
inline void trinormal(triangle *tri,double *n);
static inline void vec_in_plane(double *normal, double *v1, double *v2);



#endif COUPLING__H

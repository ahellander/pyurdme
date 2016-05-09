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

//define micro2meso micro2meso3
//define meso2micro meso2micro2

//int micro2meso1(particle *, fem_mesh *);
//int micro2meso2(particle *, fem_mesh *, int Mspecies);
//int micro2meso3(particle *particle, fem_mesh *mesh, int Mspecies);
//int micro2meso_fp(particle *particle, fem_mesh *mesh);
int micro2meso(particle *particle,vector <species>& specs, fem_mesh *mesh);


void meso2micro(particle *,int voxel, int species, double rr, double d, int dim,fem_mesh* mesh);
void meso2micro2(particle *out,fem_mesh* mesh);
void meso2micro1(particle *out,int voxel, int species, double rr, double d,int dim,fem_mesh* mesh);

void set_meso_mesh(urdme_model *model,fem_mesh *mesh,vector <voxel>& voxels,vector <plane>& boundary,simulation *sys,vector <species>& specs);
vector <plane> voxel_boundaries(urdme_model *model, fem_mesh *mesh);
inline void trinormal(triangle *tri,double *n);
inline void vec_in_plane(double *normal, double *v1, double *v2);
//static inline int nn(size_t *jcK,size_t *irK, double *p, double *pos, int pv,int spec,int Mspecies);
void check_inside_surrounding_cells(particle *p, vertex *vtx, fem_mesh *mesh, int *out);



#endif

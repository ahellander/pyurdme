/* Utility routines related to the primal and dual mesh. */

#ifndef MESH__H
#define MESH__H

#include <stdlib.h>

/* The following structures represent (extented) triangles and tetrahedra. 
   It contains information about the primal mesh as well as the dual
   cells. */

typedef struct{
	
	/* Points in xx order. */
	int v1p;
	double v1[3];
	
	int v2p;
	double v2[3];
	
	int v3p;
	double v3[3];
	
	/* Area of primal triangle*/
	double area;
	
	/* Centroid */
	double c[3];
	
	/* The volume of the contribution to the dual element 
	 associated with the vertex. */
	double a1;
	double a2;
	double a3;
	
	/* The normals and length of the edges in the dual element. The normals point
	 from point 1 towards 2 etc (dictated by the name of the normal). */
	
	double n12[3];
	double l12;
	
	double n13[3];
	double l13;
	
	double n23[3];
	double l23;
	
}triangle;


typedef struct{
	
	/* Points in xx order (primal mesh). */
	int v1p;
	double v1[3];
	int v2p;
	double v2[3];
	int v3p;
	double v3[3];
	int v4p;
	double v4[3];

	/* Volume of primal tetrahedron */
	double vol;
	
	/* Centroid of tetrahedron. This is a point in all the planes 
	   defining the dual element. */
	double c[3];
	
	/* The volume of the contribution to the dual element 
	   associated with the vertex. */
	
	double vol1;
	double vol2;
	double vol3;
	double vol4;
	
	/* The normals and area of the facets in the dual element. The normals point
	   from point 1 towards 2 etc (dictated by the name of the normal). */
	double n12[3];
	double a12;
	
	double n13[3];
	double a13;
	
	double n14[3];
	double a14;
	
	double n23[3];
	double a23;
	
	double n24[3];
	double a24;
	
	double n34[3];
	double a34;
	
	/* Frequently, we will need to calulate the barycentric coordinates of 
	   a point inside the tetrahedron. Thus, we precompute and store
	   the inverse of the matrix. */
	double T[3][3];
	
	
}tetrahedron;

typedef struct{
	
	/* Vertex coordinates (primal mesh) (3xNcells). */
	int Ncells;
	double *p;
	/* Boundary elements (primal mesh) (3xntri) */
	int ntri;
	int *e;
	/* Tetrahedra (index into p in primal mesh) (4xntet) */
	int ntet;
	int *t;
	
	/* Tetrahedra (internal format) */
	tetrahedron **tets;
	triangle **bnd;
	
	/* Connectivity matrix (needed by micro2meso,nem). */
	size_t *jcK;
	size_t *irK;
	double *prK;
	
	/* Which tetrahedra does a vertex belong to? */
	double *prv2t;
	size_t *irv2t;
	size_t *jcv2t;
	
	/* Which bundary elements does a vetrex belong to?*/
	double *prv2e;
	size_t *irv2e;
	size_t *jcv2e;
	
	/* Local mesh size */
	double *h;
	
	
}fem_mesh;

void print_tet(tetrahedron *tet);

/* Simple geometry */
double dot(double *x,double *y);
double norm(double *x);
void normalize(double *y);

/* Initialize primal/dual mesh format */
void mesh_primal2dual(fem_mesh* mesh);
void tet_index_to_internal(fem_mesh *mesh);
void tri_index_to_internal(fem_mesh *mesh);
void primal2dualtet(tetrahedron *tet);
void primal2dualtri(triangle *tri);
void tet_initialize_T(tetrahedron *tet);

/* Utility routines related to the geometry of tetrahedra/triangles */
void bary2tetpoint(double *a, tetrahedron *tet, double *p);
void bary2tripoint(double *a, triangle *tri, double *p);
int  tet_inside(tetrahedron *tet, double *x,double *minbary);
int  tri_inside(triangle *tri,double *x,double *minbary);


void tetpoint2bary(tetrahedron *tet,double *x,double *xb);
void tripoint2bary(triangle *tri,double *p,double *pb);
void triangle_randunif(triangle *tri, double *ru);
int tri_which_macro_element(triangle *tri,double *x);


double triangle_area(tetrahedron *tet);
double tetrahedron_vol(tetrahedron *tet);

void cross(double *v,double *x, double *y);
double quadarea(double *v1,double *v2,double *v3, double *v4,double *n);
double triarea(double *v1,double *v2,double *v3,double *n);


/* Random sampling */
int tet_which_macro_element(tetrahedron *tet,double *x);
void tetrahedron_randunif(tetrahedron *tet, double *ru);


#endif
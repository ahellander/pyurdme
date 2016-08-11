/* Utility routines related to the primal and dual mesh. */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"
#include "structs.h"
#include "utils.h"
#define tr 0.33333333333333333333

/* Initialize the variables and data structures relates to the
   dual mesh. This will also compute normals to the facets in
   the dual (so this need not be done again upon assembly) */

void mesh_primal2dual(fem_mesh* mesh)
{
	
	int i;
	
	tet_index_to_internal(mesh);
	tri_index_to_internal(mesh);
	
	for (i=0; i<mesh->ntet; i++) {
		primal2dualtet(mesh->tets[i]);
	}
	
	for (i=0; i<mesh->ntri; i++) {
		primal2dualtri(mesh->bnd[i]);
	}
		
}

void tet_initialize_T(tetrahedron *tet)
{
	
	double *v1,*v2,*v3,*v4;
	double b[3][3];
	double a;
	
	int i,j,k,l,maxi;
	
	int perm[] = {0,1,2};
	
	for(i=0;i<3;i++) {
		for(j=0;j<3;j++) {
			b[i][j] = 0.0;
		}
	}
	b[0][0]=1.0;
	b[1][1]=1.0;
	b[2][2]=1.0;
	
	v1 = tet->v1;
	v2 = tet->v2;
	v3 = tet->v3;
	v4 = tet->v4;
	
	/* Assemble matrix. */
	double T[3][3];
	T[0][0]=v1[0]-v4[0];
	T[1][0]=v1[1]-v4[1];
	T[2][0]=v1[2]-v4[2];
	
	T[0][1]=v2[0]-v4[0];
	T[1][1]=v2[1]-v4[1];
	T[2][1]=v2[2]-v4[2];
	
	T[0][2]=v3[0]-v4[0];
	T[1][2]=v3[1]-v4[1];
	T[2][2]=v3[2]-v4[2];
	
	
	i = 0;
	j = 0;
	while(i<3 && j<3) {
		
		/* Find pivot in column j, starting in row i. */
		maxi = i;
		for(k=maxi+1;k<3;k++) {
			if(fabs(T[k][j])>fabs(T[maxi][j])) {
				maxi = k;
			}
		}
		if(T[maxi][j] != 0) {
			/* Swap rows i and maxi */
			if(i!=maxi) {
				int temp = perm[i];
				perm[i] = perm[maxi];
				perm[maxi] = temp;
				double row_temp[] = {T[i][0],T[i][1],T[i][2]};
				T[i][0] = T[maxi][0];
				T[i][1] = T[maxi][1];
				T[i][2] = T[maxi][2];
				T[maxi][0] = row_temp[0];
				T[maxi][1] = row_temp[1];
				T[maxi][2] = row_temp[2];
				
				double row_temp2[] = {b[i][0],b[i][1],b[i][2]};
				b[i][0] = b[maxi][0];
				b[i][1] = b[maxi][1];
				b[i][2] = b[maxi][2];
				b[maxi][0] = row_temp2[0];
				b[maxi][1] = row_temp2[1];
				b[maxi][2] = row_temp2[2];
			}
			
			/* Divide each element in row i by T[i][j]. */
			a = 1/T[i][j];
			for(k=0;k<3;k++) {
				T[i][k] *= a;
				b[i][k] *= a;
			}
			/* Eliminate */
			for(k=i+1;k<3;k++) {
				/* Subtract T[k][j]*row i from row k. */
				a = T[k][j];
				for(l=0;l<3;l++) {
					T[k][l] -= a*T[i][l];
					b[k][l] -= a*b[i][l];
				}
			}
			i+=1;
		}
		else {
			printf("Matrix not invertible.\n");
		}

		j+=1;
	}
	
	/* Backward substitution */
	a=T[1][2]/T[2][2];
    for (j=0; j<3; j++) {
		T[1][j]-=(a*T[2][j]);
		b[1][j]-=(a*b[2][j]);
	}
    
	
	a=T[0][2]/T[2][2];
    for (j=0; j<3; j++) {
		T[0][j]-=(a*T[2][j]);
		b[0][j]-=(a*b[2][j]);
	}
	
	a=T[0][1]/T[1][1];
	for (j=0; j<3; j++) {
		T[0][j]-=(a*T[1][j]);
		b[0][j]-=(a*b[1][j]);
	}

	/* Store inverse */
	for (j=0; j<3; j++) {
		tet->T[0][j]=b[0][j];
		tet->T[1][j]=b[1][j];
		tet->T[2][j]=b[2][j];
	}

}

void print_T(tetrahedron *tet)
{
	printf("[");
	for(int i = 0;i<3;i++)
	{
		for(int j = 0;j<3;j++)
		{
			printf("%e ",tet->T[i][j]);
		}
		printf("\n");
	}
	printf("]\n");
}

void print_vertex(vertex *vtx){
	printf("vertex, id: %i", vtx->id);
	printf("\tcells:\n\t");

	for (int i=0;i<(int)vtx->cells.size();i++)
		cout << vtx->cells[i] << " ";
	cout << "\n";
}

void print_tet(tetrahedron *tet)
{
  printf("v1 = [%.3e %.3e %.3e]; v1p = %i ;\n",tet->v1[0],tet->v1[1],tet->v1[2],tet->v1p);	
  printf("v2 = [%.3e %.3e %.3e]; v2p = %i ;\n",tet->v2[0],tet->v2[1],tet->v2[2],tet->v2p);
  printf("v3 = [%.3e %.3e %.3e]; v3p = %i ;\n",tet->v3[0],tet->v3[1],tet->v3[2],tet->v3p);
  printf("v4 = [%.3e %.3e %.3e]; v4p = %i ;\n",tet->v4[0],tet->v4[1],tet->v4[2],tet->v4p);
  printf("center = [%.3e %.3e %.3e];\n",tet->c[0],tet->c[1],tet->c[2]);
  printf("volume = [%.3e];\n",tet->vol);
  printf("barycenter = [%.3e %.3e %.3e]; \n",tet->c[0],tet->c[1],tet->c[2]);

  printf("n12 = [%.3e %.3e %.3e]; \n",tet->n12[0],tet->n12[1],tet->n12[2]);
  printf("n13 = [%.3e %.3e %.3e]; \n",tet->n13[0],tet->n13[1],tet->n13[2]);
  printf("n14 = [%.3e %.3e %.3e]; \n",tet->n14[0],tet->n14[1],tet->n14[2]);
  printf("n23 = [%.3e %.3e %.3e]; \n",tet->n23[0],tet->n23[1],tet->n23[2]);
  printf("n24 = [%.3e %.3e %.3e]; \n",tet->n24[0],tet->n24[1],tet->n24[2]);
  printf("n34 = [%.3e %.3e %.3e]; \n",tet->n34[0],tet->n34[1],tet->n34[2]);
  //tet_initialize_T(tet);
  //print_T(tet);
}

void print_tri(triangle *tet)
{
  printf("v1 = [%.3e %.3e %.3e]; v1p = %i ;\n",tet->v1[0],tet->v1[1],tet->v1[2],tet->v1p);	
  printf("v2 = [%.3e %.3e %.3e]; v2p = %i ;\n",tet->v2[0],tet->v2[1],tet->v2[2],tet->v2p);
  printf("v3 = [%.3e %.3e %.3e]; v3p = %i ;\n",tet->v3[0],tet->v3[1],tet->v3[2],tet->v3p);
  printf("center = [%.3e %.3e %.3e];\n",tet->c[0],tet->c[1],tet->c[2]);
  printf("n12 = [%.3e %.3e %.3e]; \n",tet->n12[0],tet->n12[1],tet->n12[2]);
  printf("n13 = [%.3e %.3e %.3e]; \n",tet->n13[0],tet->n13[1],tet->n13[2]);
  printf("n23 = [%.3e %.3e %.3e]; \n",tet->n23[0],tet->n23[1],tet->n23[2]);
  
}

/* Scalar product */
/*double dot(double *x,double *y)
{
  return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
}*/

/* Euclidian distance */
double norm(double *x)
{
	return sqrt(dot(x,x));
}

void normalize(double *y)
{
	double a = norm(y);
	y[0]/=a;
	y[1]/=a;
	y[2]/=a;
}

/* Barycentric coordinate to point in tetrahedron/triangle. */
void bary2tetpoint(double *a,tetrahedron *tet,double *p)
{
	p[0]=0.0;
	p[1]=0.0;
	p[2]=0.0;

	double *x;
	
	x = tet->v1;
	
	p[0]+=a[0]*x[0];
	p[1]+=a[0]*x[1];
	p[2]+=a[0]*x[2];
	
	x = tet->v2;
	
	p[0]+=a[1]*x[0];
	p[1]+=a[1]*x[1];
	p[2]+=a[1]*x[2];
	
	x = tet->v3;
	
	p[0]+=a[2]*x[0];
	p[1]+=a[2]*x[1];
	p[2]+=a[2]*x[2];
	
	x = tet->v4;
	
	p[0]+=a[3]*x[0];
	p[1]+=a[3]*x[1];
	p[2]+=a[3]*x[2];
	
	
}

void bary2tripoint(double *a,triangle *tri,double *p)
{
	p[0]=0.0;
	p[1]=0.0;
	p[2]=0.0;
	
	double *x;
	
	x = tri->v1;
	
	p[0]+=a[0]*x[0];
	p[1]+=a[0]*x[1];
	p[2]+=a[0]*x[2];
	
	x = tri->v2;
	
	p[0]+=a[1]*x[0];
	p[1]+=a[1]*x[1];
	p[2]+=a[1]*x[2];
	
	x = tri->v3;
	
	p[0]+=a[2]*x[0];
	p[1]+=a[2]*x[1];
	p[2]+=a[2]*x[2];
	
}



void tetpoint2bary(tetrahedron *tet,double *x, double *xb)
{
	
	double b[3];
/*	printf("in tetpoint2bary 1\n");
	for(int i=0;i<3;i++){
		printf("b_%d=%g\n",i,tet->v4[i]);
	}*/
	b[0]=x[0]-tet->v4[0];
 //printf("in tetpoint2bary 2\n");
	b[1]=x[1]-tet->v4[1];
	b[2]=x[2]-tet->v4[2];
	
	xb[0]=tet->T[0][0]*b[0]+tet->T[0][1]*b[1]+tet->T[0][2]*b[2];
	xb[1]=tet->T[1][0]*b[0]+tet->T[1][1]*b[1]+tet->T[1][2]*b[2];
	xb[2]=tet->T[2][0]*b[0]+tet->T[2][1]*b[1]+tet->T[2][2]*b[2];
	xb[3]=1.0-xb[0]-xb[1]-xb[2];
	
}


/* Area of triangle */
double triangle_area(triangle *tri)
{
   	return 0.0;
}


void tripoint2bary(triangle *tri,double *p,double *pb)
{
	double *v1 = tri->v1;
	double *v2 = tri->v2;
	double *v3 = tri->v3;
	
	double a = (v2[1]-v3[1])*(p[0]-v3[0])+(v3[0]-v2[0])*(p[1]-v3[1]);
	double b = (v2[1]-v3[1])*(v1[0]-v3[0])+(v3[0]-v2[0])*(v1[1]-v3[1]);
	pb[0]=a/b;
	a = (v3[1]-v1[1])*(p[0]-v3[0])+(v1[0]-v3[0])*(p[1]-v3[1]);
	b = (v3[1]-v1[1])*(v2[0]-v3[0])+(v1[0]-v3[0])*(v2[1]-v3[1]);
	pb[1]=a/b;
	pb[2]=1.0-pb[0]-pb[1];
}

int tri_inside(triangle *tri,double *x,double *minbary)
{
	double xb[3];
 	int i;
	tripoint2bary(tri,x,xb);
	*minbary = 0.0;
	for (i=0; i<3; i++) {
		if (xb[i]<0.0) {
			*minbary += xb[i];
		}
	}
	if (*minbary<0.0) 
		return 0;
	
	return 1;
}


/* Volume of tetrahedron */
double tetrahedron_vol(tetrahedron *tet)
{
   	return 0.0;
}

/* Determine if the point x is inside (or on) the tetrahedron tet. */
int tet_inside(tetrahedron *tet, double *x,double *minbary)
{
	double xb[4];
	int i;
	tetpoint2bary(tet,x,xb);
	
	*minbary = 0.0;
	for (i=0; i<4; i++) {
		if (xb[i]<0.0) {
			*minbary += xb[i];
		}
	}
	
	if (*minbary<0.0) 
		return 0;
	
	return 1;
	
}

int tri_which_macro_element(triangle *tri, double *x)
{
	
	double *c = tri->c;
	double *n;
	
	double v[3];
	double a;
	/* A vector from the centroid to the point x.*/ 
	v[0]=x[0]-c[0]; 
	v[1]=x[1]-c[1]; 
	v[2]=x[2]-c[2];
    normalize(v);
	
	/* Normal for edge 1->2 */
	n = tri->n12;
	a = dot(n,v);
	
	if (a < 0.0) {
		n =  tri->n13;
		a = dot(n,v);
		if (a < 0.0) {
			return tri->v1p;
		}
		else {
			return tri->v3p;		
		}
	}
	else {
		n = tri->n23;
		a = dot(n,v);
		if (a < 0.0) {
			return tri->v2p;
		}
		else {
			return tri->v3p;
		}

	}
	
	/* error */
	return -1;
	
	
}

/* Knowing that a point x is inside the tetrahedron tet, 
   find which nodes macroelement it belongs to */
int tet_which_macro_element(tetrahedron *tet,double *x)
{
	/* 1.) Check which side of the plane separating 
	       node 1 and 2 the point is. */
	//printf("OK\n");
	/* The centroid is a point in all the planes. */
	double *c = tet->c;
	double *n;
	
	double v[3];
	double a;

	/* A vector from the centroid to the point.*/ 
	v[0]=x[0]-c[0]; 
	v[1]=x[1]-c[1]; 
	v[2]=x[2]-c[2];
	
	normalize(v);
	
	/* If the scalar product of the vector and the normal to the
	   plane is negative, the point is on the side belonging to 
	   vertex 1. */
	
	/* the normal for plane 1->2 */
	n = tet->n12;
	a = dot(n,v);
	
	if (a < 0.0){
		n = tet->n13;	
	   /* Now we check which side of node 3 it is on */
	    a = dot(n,v);	
		if (a < 0.0){
		 /* check point 4 */
			n = tet->n14;	
			a = dot(n,v);
			if (a < 0.0)
				return tet->v1p;
			else {
				return tet->v4p;
			}

		}
		else {
			/* On node 3's side. Check if in 3 or 4 */
			n = tet->n34;
			a = dot(n,v);
			if (a < 0.0) {
				return tet->v3p;
			}
			else{
			    return tet->v4p;	
			}
		}
	}
	else {
		/* If on nodes 2:s side, check if it is on 
		   node 3's side. */
		n = tet->n23;
		a = dot(n,v);
		if (a < 0.0) {
			/*Check if it on 4's side*/
			n = tet->n24;
			a = dot(n,v);
			if (a < 0.0) {
				return tet->v2p;
			}
			else {
				return tet->v4p;
			}

		}
		else {
			/* Check if on node 4's side */
			n = tet->n34;
			a = dot(n,v);
			if (a < 0.0){
				return tet->v3p;	
			}
			else {
				return tet->v4p;
			}

		}

		
	}
	
	/* error */
	return -1;
	
}

/* Initialize internal tetrahedron format */
void tet_index_to_internal(fem_mesh *mesh)
{
	
    printf("Entering tet_to_internal\n");
	int i,ntet;
	double *p=mesh->p;
	int *t=mesh->t;
	
	double *x;
	ntet = mesh->ntet;
    printf("here\n");
	tetrahedron **tets =(tetrahedron **)malloc(ntet*sizeof(tetrahedron*));
	for (i=0; i<ntet; i++) {
		tets[i]=(tetrahedron *)malloc(sizeof(tetrahedron));
	}
    printf("allocated\n");
    
	tetrahedron *temp;
	
	for (i=0; i<ntet; i++) {
		
	    temp=tets[i];
		
		temp->v1p=t[4*i];
		x=&p[3*temp->v1p];
		temp->v1[0]=x[0];
		temp->v1[1]=x[1];
		temp->v1[2]=x[2];
		
		temp->v2p=t[4*i+1];
		x=&p[3*temp->v2p];
		temp->v2[0]=x[0];
		temp->v2[1]=x[1];
		temp->v2[2]=x[2];
		
		temp->v3p=t[4*i+2];
		x=&p[3*temp->v3p];
		temp->v3[0]=x[0];
		temp->v3[1]=x[1];
		temp->v3[2]=x[2];
		
		temp->v4p=t[4*i+3];
		x=&p[3*temp->v4p];
		temp->v4[0]=x[0];
		temp->v4[1]=x[1];
		temp->v4[2]=x[2];
			
	}
	
	mesh->tets = tets;
	
}

/* Initialize internal triangle format */
void tri_index_to_internal(fem_mesh *mesh)
{
	
	int i,ntri;
	double *p=mesh->p;
	int *e=mesh->e;
	
	double *x;
	ntri = mesh->ntri;
	
	triangle **tri =(triangle **)malloc(ntri*sizeof(triangle *));
	for (i=0; i<ntri; i++) {
		tri[i]=(triangle *)malloc(sizeof(triangle));
	}
	triangle *temp;
	
	for (i=0; i<ntri; i++) {
	    temp=tri[i];
		
		temp->v1p=e[3*i];
		x=&p[3*temp->v1p];
		temp->v1[0]=x[0];
		temp->v1[1]=x[1];
		temp->v1[2]=x[2];
		
		
		temp->v2p = e[3*i+1];
		x=&p[3*temp->v2p];
		temp->v2[0]=x[0];
		temp->v2[1]=x[1];
		temp->v2[2]=x[2];
		
		temp->v3p=e[3*i+2];
		x=&p[3*temp->v3p];
		temp->v3[0]=x[0];
		temp->v3[1]=x[1];
		temp->v3[2]=x[2];
	
	}
	
	mesh->bnd = tri;
	
}

/* Uniform random point in triangle. */ 
void triangle_randunif(triangle *tri,double *ru)
{
	
	ru[0]=0.0;
	ru[1]=0.0;
	ru[2]=0.0;
	
	double bc[3];
	bc[0] = drand48();
	bc[1] = drand48();
	
	if (bc[0]+bc[1]>1.0) {
		bc[0] =1.0-bc[0];
		bc[1] =1.0-bc[1];
	}
	bc[2] = 1.0-bc[0]-bc[1];
	bary2tripoint(bc,tri,ru);
	
}

/* Sample point uniformly in a tetrahedron */
void tetrahedron_randunif(tetrahedron *tet,double *ru)
{
	double *x;
	
	ru[0]=0.0;
	ru[1]=0.0;
	ru[2]=0.0;
	
	double s,t,u,a;
	s = drand48();
	t = drand48();
	u = drand48();
	
	if(s+t>1.0) { // cut'n fold the cube into a prism
		
		s = 1.0 - s;
		t = 1.0 - t;
		
	}
	if(t+u>1.0) { // cut'n fold the prism into a tetrahedron
		
		double tmp = u;
		u = 1.0 - s - t;
		t = 1.0 - tmp;
		
	} else if(s+t+u>1.0) {
		
		double tmp = u;
		u = s + t + u - 1.0;
		s = 1.0 - t - tmp;
		
	}
	a=1.0-s-t-u; // a,s,t,u are the barycentric coordinates of the random point.
	
	x = tet->v1;
	ru[0]+=(a*x[0]);
	ru[1]+=(a*x[1]);
	ru[2]+=(a*x[2]);
	
	
	x = tet->v2;
	ru[0]+=(s*x[0]);
	ru[1]+=(s*x[1]);
	ru[2]+=(s*x[2]);
	
	
	x = tet->v3;
	ru[0]+=(t*x[0]);
	ru[1]+=(t*x[1]);
	ru[2]+=(t*x[2]);
	
	
	x = tet->v4;
	ru[0]+=(u*x[0]);
	ru[1]+=(u*x[1]);
	ru[2]+=(u*x[2]);
	
}

/* Cross product between vectors x,y. Result stored in v. */
void cross_mesh(double *v,double *x, double *y)
{
	
	v[0]=x[1]*y[2]-y[1]*x[2];
	v[1]=-x[0]*y[2]+y[0]*x[2];
	v[2]=x[0]*y[1]-y[0]*x[1];
	
}

/* Returns the surface area of a quadilateral specified by points v1,v2,v3,v4. 
   The surface normal is computed and stored in n */
double quadarea(double *v1,double *v2,double *v3, double *v4,double *n)
{
	
	double d1[3],d2[3];
	double a;
	
	d1[0]=v3[0]-v1[0];
	d1[1]=v3[1]-v1[1];
	d1[2]=v3[2]-v1[2];
   	
	d2[0]=v4[0]-v2[0];
	d2[1]=v4[1]-v2[1];
	d2[2]=v4[2]-v2[2];
	
	cross_mesh(n,d1,d2);
	
	/* Normalize n */
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]/=a;
	n[1]/=a;
	n[2]/=a;
	
	return 0.5*a;
	
	
}

/* From an initialized tetrahedra, calculate and initialize the part related to the dual mesh. */
void primal2dualtri(triangle *tri)
{
	
	/* Barycentric coordinates of all points that we need to compute the dual contributions. */	
	
	/* Centroid of tetrahedron*/
	double cb[] = {tr,tr,tr};	
	
	/* midpoints of edges the primal tetrahedron. */
	double mid12b[] = {0.5,0.5,0.0};	
	double mid13b[] = {0.5,0.0,0.5};
	double mid23b[] = {0.0,0.5,0.5};
	
	double v1[3],v2[3],t[3],e[3];	
	
	/* Centroid of primal triangle*/ 	
	double c[3];
	bary2tripoint(cb,tri,c);
	tri->c[0]=c[0];
	tri->c[1]=c[1];
	tri->c[2]=c[2];
		
	/* A vector in the tangent plane (normal to the primal triangle). */
	v1[0]=tri->v3[0]-tri->v1[0];
	v1[1]=tri->v3[1]-tri->v1[1];
	v1[2]=tri->v3[2]-tri->v1[2];
	
	v2[0]=tri->v2[0]-tri->v1[0];
	v2[1]=tri->v2[1]-tri->v1[1];
	v2[2]=tri->v2[2]-tri->v1[2];
	
	cross_mesh(t,v2,v1);
	
	/* Edge 12 */
	bary2tripoint(mid12b,tri,v1);
	e[0]=v1[0]-c[0];
	e[1]=v1[1]-c[1];
	e[2]=v1[2]-c[2];
	cross_mesh(tri->n12,e,t);
	normalize(tri->n12);
	/* Check if normal is pointing in the right direction. */
	v1[0]=tri->v2[0]-tri->v1[0];
	v1[1]=tri->v2[1]-tri->v1[1];
	v1[2]=tri->v2[2]-tri->v1[2];
	if (dot(v1,tri->n12)<0.0) {
		tri->n12[0]*=-1.0;
		tri->n12[1]*=-1.0;
		tri->n12[2]*=-1.0;
	}
	
	/* Edge 13 */
	bary2tripoint(mid13b,tri,v1);
	e[0]=v1[0]-c[0];
	e[1]=v1[1]-c[1];
	e[2]=v1[2]-c[2];
	cross_mesh(tri->n13,e,t);
	normalize(tri->n13);
	
	/* Check if normal is pointing in the right direction. */
	v1[0]=tri->v3[0]-tri->v1[0];
	v1[1]=tri->v3[1]-tri->v1[1];
	v1[2]=tri->v3[2]-tri->v1[2];
	if (dot(v1,tri->n13)<0.0) {
		tri->n13[0]*=-1.0;
		tri->n13[1]*=-1.0;
		tri->n13[2]*=-1.0;
	}
	
	/* Edge 23 */
	bary2tripoint(mid23b,tri,v1);
	e[0]=v1[0]-c[0];
	e[1]=v1[1]-c[1];
	e[2]=v1[2]-c[2];
	cross_mesh(tri->n23,e,t);
	normalize(tri->n23);
	
	/* Check if normal is pointing in the right direction. */
	v1[0]=tri->v3[0]-tri->v2[0];
	v1[1]=tri->v3[1]-tri->v2[1];
	v1[2]=tri->v3[2]-tri->v2[2];
	if (dot(v1,tri->n23)<0.0) {
		tri->n23[0]*=-1.0;
		tri->n23[1]*=-1.0;
		tri->n23[2]*=-1.0;
	}

}



void primal2dualtet(tetrahedron *tet)
{
	
	/* Barycentric coordinates of all points that we need to compute the dual contributions. */	
	
	/* Centroid of tetrahedron*/
	double cb[] = {0.25,0.25,0.25,0.25};	
	
	/* Centroid of triangle faces in the primal tetrahedron. */
	double c134b[] = {tr, 0.0, tr, tr};	
	double c124b[] = {tr, tr, 0, tr};
	double c123b[] = {tr, tr, tr, 0};
	double c234b[] = {0.0, tr, tr, tr};
	
	/* Midpoint on edges of the primal tetrahedron. */	
	double mid12b[] = {0.5, 0.5, 0.0, 0.0};  	
	double mid13b[] = {0.5, 0.0, 0.5, 0.0};
	double mid23b[] = {0.0, 0.5, 0.5, 0.0};
	double mid14b[] = {0.5, 0.0, 0.0, 0.5};
	double mid24b[] = {0.0, 0.5, 0.0, 0.5};
	double mid34b[] = {0.0, 0.0, 0.5, 0.5};
	
	double v1[3],v2[3],v3[3];	
	/* Centroid */ 	
	double c[3];
	bary2tetpoint(cb,tet,c);
	
	tet->c[0]=c[0];
	tet->c[1]=c[1];
	tet->c[2]=c[2];
	
	double a;	
	
	/* Facet 12 */
	bary2tetpoint(mid12b,tet,v1);
	bary2tetpoint(c123b,tet,v2);
	bary2tetpoint(c124b,tet,v3);
	
	a = quadarea(v1,v2,c,v3,tet->n12);	

   	
	tet->a12 = a;
	
	/* Facet 13 */
	bary2tetpoint(mid13b,tet,v1);
	bary2tetpoint(c134b,tet,v2);
	bary2tetpoint(c123b,tet,v3);
	
	
	a = quadarea(v1,v2,c,v3,tet->n13);	
	tet->a13 = a;
	
	/* Facet 23 */
	bary2tetpoint(mid23b,tet,v1);
	bary2tetpoint(c123b,tet,v2);
	bary2tetpoint(c234b,tet,v3);
	
	a = quadarea(v1,v2,c,v3,tet->n23);	
	tet->a23 = a;
	
	/* Facet 14 */
	bary2tetpoint(mid14b,tet,v1);
	bary2tetpoint(c124b,tet,v2);
	bary2tetpoint(c134b,tet,v3);
	
	a = quadarea(v1,v2,c,v3,tet->n14);	
	tet->a14 = a;
	
	/* Facet 24 */
	bary2tetpoint(mid24b,tet,v1);
	bary2tetpoint(c234b,tet,v2);
	bary2tetpoint(c124b,tet,v3);
	
	a = quadarea(v1,v2,c,v3,tet->n24);	
	tet->a24 = a;
	
	/* Facet 34 */
	bary2tetpoint(mid34b,tet,v1);
	bary2tetpoint(c134b,tet,v2);
	bary2tetpoint(c234b,tet,v3);
	
	a = quadarea(v1,v2,c,v3,tet->n34);	
	tet->a34 = a;
	
	tet_initialize_T(tet);
	
}






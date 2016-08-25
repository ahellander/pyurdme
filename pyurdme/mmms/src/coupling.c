/* Utility routines for micro/meso coupling */

/* Micro particle typedef */
//#include "micro.h"
#include "coupling.h"
#include "mesh.h"
#include "structs.h"
#include "utils.h"
#include <math.h>

#define pi 3.141592653589793

using namespace std;


/* Some utility routines for geometrical properties of mesh elements. */

/* Distance between two points. This routine will be used frequently, 
   and eventually we should optimize this as much as possible. */ 
static inline double distance(double *x,double *y)
{
  return (x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]); 
}

void set_meso_mesh(urdme_model *model,fem_mesh *mesh,vector <voxel>& voxels,vector <plane>& boundary,simulation *sys,vector <species>& specs)
{
    int Ncells = mesh->Ncells;
    int Mspecies = model->Mspecies;
    double *p = mesh->p;
    size_t *jcD = model->jcD;
    size_t *irD = model->irD;
    double *prD = model->prD;
    
    
    voxel temp;
    neighbor tn;
    int to_dof;
    for(int i=0;i<Ncells;i++){
        temp.id = i;
        temp.node[0] = boundary[i].p[0];
        temp.node[1] = boundary[i].p[1];
        temp.node[2] = boundary[i].p[2];
        
        
        int dof = i*Mspecies;

        temp.num_conn=0;
        temp.vol = model->vol[i];
        temp.sd = model->sd[i];
        temp.totalD = 0.0;
        temp.neighbors.resize(0);
        temp.isbnd = boundary[i].isbnd;
        
        
        for (int k=jcD[dof]; k<jcD[dof+1]; k++) {
            
            
            to_dof=irD[k];
            if(to_dof!=dof){
                temp.num_conn++;
                tn.vox = to_dof/Mspecies;
                tn.D = prD[k]/specs[0].D;
                temp.totalD += tn.D;
                temp.neighbors.push_back(tn);
            }
            
        }
        voxels.push_back(temp);
        
//        printf("\n");
//        printf("Ncells=%d,Mspecies=%d\n",Ncells,model->Mspecies);
    }
}

vector <plane> voxel_boundaries(urdme_model *model, fem_mesh *mesh)
{
   		
	int Ncells   = mesh->Ncells;
//	int Mspecies = model->Mspecies;
	 
	/* Boundary triangles */
	double *p = mesh->p;
//	int *sd = model->sd;
	
	vector <plane> boundaries;
	
	int ntri=mesh->ntri; 
	int i;//,j;
	
	double a;

	/* Initialize all normals and vectors */
	plane temp;
	
	temp.n[0]=0.0;
	temp.n[1]=0.0;
	temp.n[2]=0.0;
	temp.v1[0]=0.0;
	temp.v1[1]=0.0;
	temp.v1[2]=0.0;
	temp.v2[0]=0.0;
	temp.v2[1]=0.0;
	temp.v2[2]=0.0;
	temp.p[0]=0.0;
	temp.p[1]=0.0;
	temp.p[2]=0.0;
	temp.isbnd = 0;
	
	for (i=0; i<Ncells; i++) {
		/* Point in plane, the vertex itself. */ 
		temp.id = i;
		boundaries.push_back(temp);
	}
	
	/* Compute the point in the plane as the mean value of the vertex iteself and the 
	 centroid of the surrounding triangles. */
	for (i=0; i<Ncells; i++) {
		boundaries[i].p[0]=p[3*i];
		boundaries[i].p[1]=p[3*i+1];
		boundaries[i].p[2]=p[3*i+2];
	}
	
	
	/* Loop over the boundary triangles and assemble their contribution to the normals */
	triangle *tri;
    double n[3];
	double vec[3];

	int sv;
//   double radius=0;	
	
	size_t *jcD = model->jcD;
	size_t *irD = model->irD;
	
    int dof,to_dof;  
	size_t k;
	double y[3];
	double alpha;

	/* We will need to be able determine whether a dof is on the boundary/surface or not. */
	int *onboundary;
	onboundary = (int *)calloc(Ncells,sizeof(int));
	
	for (i=0; i<ntri; i++) {
	    tri = mesh->bnd[i];
		onboundary[tri->v1p]=1;
		onboundary[tri->v2p]=1;
		onboundary[tri->v3p]=1;
	}
	
	double amax;
	
	for (i=0; i<ntri; i++) {
		
		tri = mesh->bnd[i];
		
		/* Compute the unit normal of that triangle.*/
		trinormal(tri,n);

		/* TODO: Need to check that all normals that we assemble point in a direction that
		   the species may move. */		 
//		radius = norm(tri->v1);
		a = dot(vec,n);
		
		/* Add that normal to accumulated mean normals of the vertices of the triangle.*/
		sv = tri->v1p;
	
			
		alpha = 1.0;
		dof = sv;
			
		/* Find the edge with the largest dot-product with the normal. */
		a=0.0;
		amax =0.0;
		dof = sv;
		/* If the point has neighbours in the interior */
		for (k=jcD[dof]; k<jcD[dof+1]; k++) {
				
			to_dof=irD[k];
			
			y[0]=p[3*to_dof]  -p[3*dof];
			y[1]=p[3*to_dof+1]-p[3*dof+1];
			y[2]=p[3*to_dof+2]-p[3*dof+2];
				
			normalize(y);
			a = dot(n,y);
			if(fabs(a)>fabs(amax)){
				amax = a;
			}
				
		}
		
		/* If the edge in D with the largest component in the normals direction
		  points in the opposite direction, we flip the normal. */
		if (amax < 0.0) {
			alpha = -1.0;
		}

		boundaries[sv].n[0]+=alpha*n[0];
		boundaries[sv].n[1]+=alpha*n[1];
		boundaries[sv].n[2]+=alpha*n[2];
		boundaries[sv].isbnd=1;
		
			
		sv = tri->v2p;
			
		alpha = 1.0;
			
		/* Find the edge with the largest dot-product with the normal. */
		a=0.0;
		amax =0.0;
		dof = sv;

		for (k=jcD[dof]; k<jcD[dof+1]; k++) {
				
			to_dof=irD[k];
			y[0]=p[3*to_dof]  -p[3*dof];
			y[1]=p[3*to_dof+1]-p[3*dof+1];
			y[2]=p[3*to_dof+2]-p[3*dof+2];

			normalize(y);
			a = dot(n,y);
			if(fabs(a)>fabs(amax))
				amax = a;
				
		}
			
		/* If the edge in D with the larges component in the normals direction
		 points in the opposite direction, we flip the normal. */
		if (amax < 0.0) {
			alpha = -1.0;
		}
			
//			printf("n3=[%g %g %g];\n",alpha*n[0],alpha*n[1],alpha*n[2]);

		boundaries[dof].n[0]+=alpha*n[0];
		boundaries[dof].n[1]+=alpha*n[1];
		boundaries[dof].n[2]+=alpha*n[2];
		boundaries[dof].isbnd=1;
		
		
		sv = tri->v3p;

		/* Find the edge with the largest dot-product with the normal. */
		alpha = 1.0;
		dof = sv;
		a=0.0;
		amax =0.0;
			
		for (k=jcD[dof]; k<jcD[dof+1]; k++) {
			
			to_dof=irD[k];
			y[0]=p[3*to_dof]  -p[3*dof];
			y[1]=p[3*to_dof+1]-p[3*dof+1];
			y[2]=p[3*to_dof+2]-p[3*dof+2];
				
			normalize(y);
			a = dot(n,y);
			if(fabs(a)>fabs(amax))
				amax = a;
				
		}
			
		/* If the edge in D with the largest component in the normals direction
			points in the opposite direction, we flip the normal. */
		if (amax < 0.0) {
			alpha = -1.0;
		}
		boundaries[sv].n[0]+=alpha*n[0];
		boundaries[sv].n[1]+=alpha*n[1];
		boundaries[sv].n[2]+=alpha*n[2];
		boundaries[sv].isbnd=1;


	}
	
	/* Normalize the mean normals */
	plane *tmp;
	for (i=0; i<Ncells; i++) {
		
		tmp = &boundaries[i];
		
		a = norm(tmp->n);
		if (a>0.0){
			tmp->n[0]/=a;
			tmp->n[1]/=a;
			tmp->n[2]/=a;
		}
		
	}

	/* Get vectors in the tangent plane */
	for (i=0; i<Ncells;i++) {
		if (boundaries[i].isbnd){ 
			vec_in_plane(boundaries[i].n,boundaries[i].v1,boundaries[i].v2);
		}
		boundaries[i].type.push_back(0);
	}
		
	free(onboundary);

	return boundaries;	

}

/* Normalized surface normal to the plane defined by the triangle tri. */
static inline double *tri_surfnormal(int *tri,double *p)
{
	double a[3];
	double b[3];
	double *n;
	n = (double *)malloc(3*sizeof(double));
	
	double *x1,*x2,*x3;
	x1 = &p[tri[0]];
	x2 = &p[tri[1]];
	x3 = &p[tri[2]];
	
	a[0]=x2[0]-x1[0];
	a[1]=x2[1]-x1[1];
	a[2]=x2[2]-x1[2];

	b[0]=x3[0]-x1[0];
	b[1]=x3[1]-x1[1];
	b[2]=x3[2]-x1[2];
	
	n[0]=a[1]*b[2]-b[1]*a[2];
	n[1]=-(a[0]*b[2]-b[0]*a[2]);
	n[2]=a[0]*b[1]-b[0]*a[1];
	
    double d = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/d;
	n[1]=n[1]/d;
	n[2]=n[2]/d;
	
	return n;
	
}



static inline int cmp_point(double *v1,double *v2)
{
    int i;
   for (i=0; i<3; i++) {
	   if (v1[i]!=v2[i]) {
		   return -1;
	   }
   }
   return 1;	
	
}

void project(double *p,triangle *tri) {
	double *v1 = tri->v1;
	double *v2 = tri->v2;
	double *v3 = tri->v3;
	double *vec1 = (double *)malloc(3*sizeof(double));
	double *vec2 = (double *)malloc(3*sizeof(double));
	vec1[0] = v2[0]-v1[0];
	vec1[1] = v2[1]-v1[1];
	vec1[2] = v2[2]-v1[2];
	vec2[0] = v3[0]-v2[0];
	vec2[1] = v3[1]-v2[1];
	vec2[2] = v3[2]-v2[2];
	double *normal = (double *)malloc(3*sizeof(double));
	cross_mesh(normal,vec1,vec2);
	normalize(normal);
	double dist = normal[0]*(p[0]-v3[0])+normal[1]*(p[1]-v3[1])+normal[2]*(p[2]-v3[2]);
	p[0] -= dist*normal[0];
	p[1] -= dist*normal[1];
	p[2] -= dist*normal[2];
	free(normal);
	free(vec1);
	free(vec2);
}

static inline int nn(size_t *jcK,size_t *irK, double *p, double *pos, int pv,int spec,int Mspecies){
	
	int node,node2,nv;
	size_t j,k;
	double d,dn;
	double *x;
		
	//nv=pv*Mspecies+spec;
	//d = INFINITY;
	
	int dof = pv*Mspecies + spec;
	nv = dof;
     
	/* Distance from the particle to the current guess for the voxel. */
	x = &p[3*pv];
	d = distance(x,pos);
	//printf("d: %.2e\n",d);
	//printf("x: (%.2e %.2e %.2e)\n",x[0],x[1],x[2]);
	//printf("pos: (%.2e %.2e %.2e)\n",pos[0],pos[1],pos[2]);
	/* To be sure to find the nearest voxels in 3D we need to look two layers out from the voxel. */

	for (j=jcK[dof]; j<jcK[dof+1]; j++) {
		
		node = irK[j];
		
		/* Check distance to all neighbours to node */ 
		x  = &p[3*(node-spec)/Mspecies];
		dn = distance(x,pos);
		
		//printf("dn: %.2e\n",dn);
		
		if (dn < d) {
			nv = node;
			d = dn;
		}
		
		/* And to all the neighbours of the neighbour (we are comparing to the same point many times this way...). */
		for (k=jcK[node]; k<jcK[node+1]; k++) {
			
			node2 = irK[k];
			
			// Check distance to all neighbours to node 
			x = &p[3*(node2-spec)/Mspecies];
			dn = distance(x,pos);
			//printf("dn: %.2e\n",dn);

			if (dn < d) {
				nv = node2;
				d = dn;
			}
		}
		
		
	}
	return (nv-spec)/Mspecies;
}



/* Exact sampling of points on dual elements. */
void meso2micro2(particle *out,fem_mesh* mesh)
{
    
//    double *p = mesh->p;
    
    double *prv2t = mesh->prv2t;
    size_t *jcv2t = mesh->jcv2t;
    
    double *prv2e = mesh->prv2e;
    size_t *jcv2e = mesh->jcv2e;
    
    int voxel = out->voxel;
//    printf("voxel=%d\n",voxel);
    int dim = out->dim;
//    out->type   = species;
//    out->voxel  = voxel;
//    out->status = 1;
//    out->dim  = dim;
//    
//    out->sigma = sigma;
//    out->D  = D;
    
    tetrahedron **tets = mesh->tets;
    vertex **vertices = mesh->vertices;
    vertex *vtx;
    triangle    **tris = mesh->bnd;
    
    double ru[3];
    int wm;
    tetrahedron *temp;
    double minbary[1];
    
    /* If the particle is a 3D species */
    if (1==1){//dim == 3){
//        printf("dim=%d\n",dim);
        /* Number of tetrahedra around that vertex. */
//        int ntet = jcv2t[voxel+1]-jcv2t[voxel];
//        int ntet = (int)neighbors.size();
        vtx = vertices[voxel];
        int ntet = (int)(vtx->cells.size());
        
        
//        printf("ntet=%d\n",ntet);
        /* Pick a "random" tetrahedra. -1 to compensate for Matlab indexing starting at 1 instead of 0. */
//        int tet  = (int) prv2t[jcv2t[voxel]+(int)(drand48()*ntet)]-1;
        int tet = vtx->cells[(int)(drand48()*ntet)];
        /* Randomly in that tetrahedron */
        
        tetrahedron_randunif(tets[tet],ru);
        temp = tets[tet];
        wm = tet_which_macro_element(tets[tet],ru);
        if (wm==-1) {
            printf("Problem in whichmacroelement.\n");
        }
        
        
        /* Rejection sampling */
        while (wm!=voxel) {
            //			printf("tet=%d\n",tet);
            tetrahedron_randunif(tets[tet],ru);
            if (!tet_inside(tets[tet],ru,minbary)) {
                printf("Sampled point not inside tetrahedron. This is bad.\n");
            }
            wm = tet_which_macro_element(tets[tet],ru);
            if (wm==-1) {
                printf("Problem in whichmacroelement.\n");
            }
        }
        
        out->pos[0]=ru[0];
        out->pos[1]=ru[1];
        out->pos[2]=ru[2];
        
        
        
    }
    /* If 2D species, place randomly on the surface mesh */
    else if (dim == 2){
        
        int ntri = jcv2e[voxel+1]-jcv2e[voxel];
        int tri  = (int) prv2e[jcv2e[voxel]+(int)(drand48()*ntri)]-1;
        
        triangle_randunif(tris[tri],ru);
        
        if (!tri_inside(tris[tri],ru,minbary)) {
            printf("Sampled point not inside tetrahedron. This is bad.\n");
        }
        wm = tri_which_macro_element(tris[tri],ru);
        if (wm==-1) {
            printf("Problem in whichmacroelement.\n");
        }
        
        while (wm!=voxel) {
            triangle_randunif(tris[tri],ru);
            if (!tri_inside(tris[tri],ru,minbary)) {
                printf("Sampled point not inside tetrahedron. This is bad.\n");
            }
            wm = tri_which_macro_element(tris[tri],ru);
            if (wm==-1) {
                printf("Problem in whichmacroelement.\n");
            }
        }
        
        out->pos[0]=ru[0];
        out->pos[1]=ru[1];
        out->pos[2]=ru[2];
        
        
    }
    else {
        printf("1D species are not yet supported. dim = %d\n",dim);
        exit(-1);
    }
    
}



/* A new, faster micro2meso. */
int micro2meso2(particle *particle, fem_mesh *mesh, int Mspecies)
{
	
	// The connectivity matrix 
	size_t *irK,*jcK;
	irK = mesh->irK;
	jcK = mesh->jcK;
	
	// mesh vetrices 
	double *p;
	p = mesh->p;
	
	// The particles current position 
	double *pos = particle->pos;
	
	// Check which vertex the particle was nearest in the previus step.
	int pv = particle->voxel;
	
	// Find the nearest neighbour among pv's nearest neighbours in the mesh. 
	int nv;
	int spec = particle->type;

	nv = nn(jcK,irK,p,pos,pv,spec,Mspecies);
    
	//printf("pv %i nv %i\n",pv,nv);
	
	// Slide towards the nearest voxel. 
	while (nv != pv) {
		pv = nv;
		nv = nn(jcK,irK,p,pos,pv,spec,Mspecies);
	//	printf("pv %i nv %i\n",pv,nv);

	}
	
	return nv;

}

void check_inside_surrounding_cells(particle *particle, vertex *vtx, fem_mesh *mesh, int *out)
{
		
	    tetrahedron **tets = mesh->tets;
		double minbary[1];
		int wn = -1;
		tetrahedron *tet;
		double mindist=-INFINITY;
		int closest_tet=-1;

//		printf("Number of cells: %d\n",(int)(vtx->cells.size())); 

		for (int k=0; k<(int)(vtx->cells.size()); k++) {
//			printf("numtets=%d, k=%d\n",mesh->ntet,k);			
//			printf("index=%d\n",vtx->cells[k]);
			tet = tets[vtx->cells[k]];
//			printf("tet=%d\n",tet);
			// We need to check if there is a non-zero edge in the connectivity matrix
			// between our vertex and the others in that tetrahedron,
			// otherwise the particle should be restricted from moving there (Not yet implemented).
			
			if (tet_inside(tet,particle->pos,minbary)) {
				wn=tet_which_macro_element(tet,particle->pos);
				break;
			}
			else {
				if (*minbary > mindist) {
					mindist = *minbary;
					closest_tet = vtx->cells[k];
				}
			}

		}

		out[0] = wn; 
		out[1] = closest_tet;

}

/* "Exact" version of micro2meso. */
int micro2meso(particle *particle, vector <species>&specs, fem_mesh *mesh)
{

	
	double *p = mesh->p;
	
	tetrahedron **tets = mesh->tets;
	vertex **vertices = mesh->vertices;
	
	/* If a 3D species, start looking in the tetrahedra that
       covers the macroelement associated with the vertex the particle
	   belonged to in the previous timestep.  */
	
	//int j,k,voxel;
    int voxel;
	int wn = -1;
	tetrahedron *tet;
	//int closest_tet;
	double mindist=INFINITY;
	double dist;

	double minbary[1];
	//int spec=particle->type;
	
	//int i;
	vertex *vtx;
	int out[2];
	//if (particle->dim==3 || particle->dim==1) {
    
    /* TODO: All particles are 3D at the moment. Update this when we introduce membrane-bound particles. */
    
    if (1==1){//specs[particle->type].dim==3){

		voxel = particle->voxel;
		vtx = vertices[voxel];

		// Check if particle is inside one of the tetrahedra surrounding the previous voxel.  
		check_inside_surrounding_cells(particle, vtx, mesh, out);

		if (out[0] > 0){
			particle->voxel=out[0];
			return particle->voxel;
		}

		// If it is not inside one of those, we search to find the closest vertex in the mesh, 
		// TODO: THIS IS JUST TEMPORARY, WE NEED TO RE-INTEGRATE THE STRATEGY WITH SLIDING FORM MICRO2MESO2
		for (int i=0; i<mesh->Ncells; i++){
			dist = dist3(&p[i*3],particle->pos);
			if (dist<mindist){
				mindist = dist;
				wn = i;
			}
		}

		// 
		particle->voxel = wn;
		// Now we check again if the point is contained in one of the tetrahedra surrounding the closest vertex. 
		vtx=vertices[wn];

		check_inside_surrounding_cells(particle, vtx, mesh, out);
		if (out[0] > 0){
			particle->voxel = out[0];
			return out[0];
		}

		// If not, we are in trouble. We now resort to looking in every tetrahedron in the mesh. 
		// This is very expensive, so the timestep of the micro sovler should be chosen so that this only happens rarely.  
		for (int i=0; i<mesh->ntet;i++){
			tet = tets[i];
			if (tet_inside(tet,particle->pos,minbary)) {
				wn=tet_which_macro_element(tet,particle->pos);
				break;
			}
		}
		particle->voxel=wn;
		// At this point, this should not happen. 
		if (wn == -1)
			cout << "Warning, did not find exact voxel for particle.\n";
		
		return particle->voxel;
		
	}

	/*else if (particle->dim==2) {
	
		triangle **tris = mesh->bnd;
		voxel = particle->voxel;
		triangle *tri;
		int closest_tri;
		double *p_temp = (double *)malloc(3*sizeof(double));
		// Check if particle is inside one of the surrounding tetrahedra. 
		for (k=jcv2e[voxel]; k<jcv2e[voxel+1]; k++) {
			
			tri = tris[(int)prv2e[k]-1];
			
			// We need to check if there is a non-zero edge in the connectivity matrix
			// between our vertex and the others in that tetrahedron,
			// otherwise the particle should be restricted from moving there (Not yet implemented).//
			p_temp[0] = particle->p[0];
			p_temp[1] = particle->p[1];
			p_temp[2] = particle->p[2];
			project(p_temp,tri);
			if (tri_inside(tri,p_temp,minbary)) {
				wn=tri_which_macro_element(tri,p_temp);
				free(p_temp);
				return wn;
			}
			else {
				if (*minbary > mindist) {
					mindist = *minbary;
					closest_tri = (int)prv2e[k]-1;
				}
			}
		}
	
		nv = micro2meso2(particle,mesh,Mspecies);
		
		for (k=jcv2e[nv]; k<jcv2e[nv+1]; k++) {
			
			tri = tris[(int)prv2e[k]-1];
			
			p_temp[0] = particle->p[0];
			p_temp[1] = particle->p[1];
			p_temp[2] = particle->p[2];
			project(p_temp,tri);
			if (tri_inside(tri,p_temp,minbary)){
				wn=tri_which_macro_element(tri,p_temp);
				free(p_temp);
				return wn;
			}
			else {
				if (*minbary > mindist) {
					mindist = *minbary;
					closest_tri = (int)prv2e[k]-1;
				}
			}
			
		}
		
		// If this too failed, look in the tetrahedra around the neighbouring vertices.  
		int spec = particle->type;
		int dof = Mspecies*nv+spec;
		for (j=jcK[dof]; j<jcK[dof+1]; j++) {
			
			for (k=jcv2e[(irK[j]-spec)/Mspecies]; k<jcv2e[(irK[j]-spec)/Mspecies+1]; k++) {
				
				tri = tris[(int)prv2e[k]-1];
				p_temp[0] = particle->p[0];
				p_temp[1] = particle->p[1];
				p_temp[2] = particle->p[2];
				project(p_temp,tri);
				if (tri_inside(tri,p_temp,minbary)){
					wn=tri_which_macro_element(tri,p_temp);
					free(p_temp);
					return wn;
				}
				
			}
			
		}
		
		// If this failed, assign the particle to the "closest" tetrahedron, 
		// as determinied by the sum of the negative components of the barycentric
		// coordinates. 
		
		tri = tris[closest_tri];
		p_temp[0] = particle->p[0];
		p_temp[1] = particle->p[1];
		p_temp[2] = particle->p[2];
		project(p_temp,tri);
		wn = tri_which_macro_element(tri,p_temp);
		free(p_temp);
		return wn;
		
	}
 	*/
	else {
		printf("Only 3D species are supported in micro2meso.\n");
        printf("Dimension is: %d. Particle type is: %d\n",specs[particle->type].dim,particle->type);
	}
	
	return -1;

	   
}


static inline double *tet_randunif(int *tet,double *p)
{
	double *ru,*x;
	ru = (double *)malloc(3*sizeof(double));
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

	x = &p[3*tet[0]];
	ru[0]+=(a*x[0]);
	ru[1]+=(a*x[1]);
	ru[2]+=(a*x[2]);


	x = &p[3*tet[1]];
	ru[0]+=(s*x[0]);
	ru[1]+=(s*x[1]);
	ru[2]+=(s*x[2]);


	x = &p[3*tet[2]];
	ru[0]+=(t*x[0]);
	ru[1]+=(t*x[1]);
	ru[2]+=(t*x[2]);


	x = &p[3*tet[3]];
	ru[0]+=(u*x[0]);
	ru[1]+=(u*x[1]);
	ru[2]+=(u*x[2]);

	return ru;

}

/* Make sure that the vectors in the plane make up an orthonormal system. */
inline void orthonormal_plane(plane *pl)
{
	double x[3];
	double a,b;
	/* First vector. */
	
	/* The normal vector is already normalized. Projection of normal on v1 */
	a = dot(pl->n,pl->v1);
	
	x[0]=pl->v1[0]-a*pl->n[0];
	x[1]=pl->v1[1]-a*pl->n[1];
	x[2]=pl->v1[2]-a*pl->n[2];
	
	a = sqrt(dot(x,x));

	pl->v1[0]=x[0]/a;
	pl->v1[1]=x[1]/a;
	pl->v1[2]=x[2]/a;
	
	/* Other vector */
	a = dot(pl->n,pl->v2);
	b = dot(pl->v1,pl->v2);
	
	x[0]=pl->v2[0]-a*pl->n[0]-b*pl->v1[0];
	x[1]=pl->v2[1]-a*pl->n[1]-b*pl->v1[1];
	x[2]=pl->v2[2]-a*pl->n[2]-b*pl->v1[2];
	
	a = sqrt(dot(x,x));
	
	pl->v2[0]=x[0]/a;
	pl->v2[1]=x[1]/a;
	pl->v2[2]=x[2]/a;
	
}


inline void trinormal(triangle *tri,double *n)
{

	double a=0.0;
	double v1[3],v2[3];

	/* Make sure that n is zero. */
	n[0]=0.0;
	n[1]=0.0;
	n[2]=0.0;
	
	/* Vectors in plane */
	v1[0]=tri->v2[0]-tri->v1[0];
	v1[1]=tri->v2[1]-tri->v1[1];
	v1[2]=tri->v2[2]-tri->v1[2];
	
	v2[0]=tri->v3[0]-tri->v2[0];
	v2[1]=tri->v3[1]-tri->v2[1];
	v2[2]=tri->v3[2]-tri->v2[2];
	
	/* Normal: cross product of v1 and v2. */
	cross_mesh(n,v1,v2);
	
	/* normalize */
	a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]/=a;
	n[1]/=a;
	n[2]/=a;
	
}

inline void vec_in_plane(double *normal, double *v1, double *v2)
{
	double r_dir[] = {normal[0]+1,2*normal[1]+2,3*normal[2]+3};
	cross_mesh(v1,normal,r_dir);
	cross_mesh(v2,normal,v1);
	double norm_n = 1.0/norm(normal);
	double norm_v1 = 1.0/norm(v1);
	double norm_v2 = 1.0/norm(v2);
	normal[0] = norm_n*normal[0];
	normal[1] = norm_n*normal[1];
	normal[2] = norm_n*normal[2];
	v1[0] = norm_v1*v1[0];
	v1[1] = norm_v1*v1[1];
	v1[2] = norm_v1*v1[2];
	v2[0] = norm_v2*v2[0];
	v2[1] = norm_v2*v2[1];
	v2[2] = norm_v2*v2[2];
}







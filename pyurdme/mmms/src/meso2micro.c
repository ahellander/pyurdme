*
   Create a microparticle from a mesoscopic dof. Place the particle randomly
   within the mesh voxel.

   A first crude implementation. Places the particle in the centroid of
   a neighbouring tetrahedra or triangle.

   */

void meso2micro1(particle *out,int voxel, int species, double rr, double d,int dim,fem_mesh* mesh)
{

	double *p = mesh->p;
	int *t = mesh->t;
	int *e = mesh->e;

	double *prv2t = mesh->prv2t;
	size_t *irv2t = mesh->irv2t;
	size_t *jcv2t = mesh->jcv2t;

	double *prv2e = mesh->prv2e;
	size_t *irv2e = mesh->irv2e;
	size_t *jcv2e = mesh->jcv2e;


	out->type   = species;
	out->voxel  = voxel;
	out->status = 1;
	out->dim  = dim;

	out->rr = rr;
	out->d  = d;

	int i,j,k;

	double *x;
	double center[3];

	center[0]=0.0;
	center[1]=0.0;
	center[2]=0.0;


	/* If the particle is a 3D species */

	if (dim == 3){

		/* Number of tetrahedra around that vertex. */
		int ntet = jcv2t[voxel+1]-jcv2t[voxel];
		/* Pick a "random" tetrahedra. -1 to compensate for Matlab indexing starting at 1 instead of 0. */
		int tet  = (int) prv2t[jcv2t[voxel]+(int)(drand48()*ntet)]-1;

		/* Randomly in that tetrahedron */
		x = tet_randunif(&t[4*tet],p);
		out->p[0]=x[0];
		out->p[1]=x[1];
		out->p[2]=x[2];

		/* Put in vertex itself - debug... */
		out->p[0]=p[3*voxel];
		out->p[1]=p[3*voxel+1];
		out->p[2]=p[3*voxel+2];

		free(x);

	}
	/* If 2D species, place randomly on the surface mesh */
	else if (dim == 2){
		int ntri = jcv2e[voxel+1]-jcv2e[voxel];
		int tri  = (int) prv2e[jcv2e[voxel]+(int)(drand48()*ntri)]-1;

		for (i=0; i<3; i++) {
			x = &p[3*e[3*tri+i]];
			center[0]+=(x[0]/3.0);
			center[1]+=(x[1]/3.0);
			center[2]+=(x[2]/3.0);
		}

		/* Put in vertex itself - debug... */
		out->p[0]=p[3*voxel];
		out->p[1]=p[3*voxel+1];
		out->p[2]=p[3*voxel+2];


		/*out->p[0]=center[0];
		 out->p[1]=center[1];
		 out->p[2]=center[2];*/


	}
	else {
		printf("1D species are not yet supported. dim = %d\n",dim);
		exit(-1);
	}

}

/* Exact sampling of points on dual elements. */
void meso2micro2(micro_particle *out,int voxel, int species, double rr, double d,int dim,fem_mesh* mesh)
{
	
	double *p = mesh->p;
	
	double *prv2t = mesh->prv2t;
	size_t *jcv2t = mesh->jcv2t; 
	
	double *prv2e = mesh->prv2e;
	size_t *jcv2e = mesh->jcv2e; 
	
	
	out->type   = species;
	out->voxel  = voxel;
	out->status = 1;
	out->dim  = dim;
	
	out->rr = rr;
	out->d  = d;
	
	tetrahedron **tets = mesh->tets;
	triangle    **tris = mesh->bnd;
		
	double ru[3];
	int wm;
	tetrahedron *temp;
	double minbary[1];
	
	/* If the particle is a 3D species */
	if (dim == 3){
		
		/* Number of tetrahedra around that vertex. */	
		int ntet = jcv2t[voxel+1]-jcv2t[voxel]; 
		/* Pick a "random" tetrahedra. -1 to compensate for Matlab indexing starting at 1 instead of 0. */
		int tet  = (int) prv2t[jcv2t[voxel]+(int)(drand48()*ntet)]-1;
		/* Randomly in that tetrahedron */
		
		tetrahedron_randunif(tets[tet],ru);
		temp = tets[tet];
		wm = tet_which_macro_element(tets[tet],ru);
		if (wm==-1) {
			printf("Problem in whichmacroelement.\n");
		}
		
		
		/* Rejection sampling */
		while (wm!=voxel) {
			tetrahedron_randunif(tets[tet],ru);
			if (!tet_inside(tets[tet],ru,minbary)) {
			  printf("Sampled point not inside tetrahedron. This is bad.\n");
			}
			wm = tet_which_macro_element(tets[tet],ru);
            if (wm==-1) {
				printf("Problem in whichmacroelement.\n");
			}
		}
		
		out->p[0]=ru[0];
		out->p[1]=ru[1];
		out->p[2]=ru[2];
	
		
		
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
		
		out->p[0]=ru[0];
		out->p[1]=ru[1];
		out->p[2]=ru[2];
		
		
	}
	else {
		printf("1D species are not yet supported. dim = %d\n",dim);
		exit(-1);
	}
	
}

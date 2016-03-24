/*
   Create a microparticle from a mesoscopic dof. Place the particle randomly
   within the mesh voxel.

   A first crude implementation. Places the particle in the centroid of
   a neighbouring tetrahedra or triangle.

*/
vector <plane> voxel_boundaries2(urdme_model *model, fem_mesh *mesh)
{
    
    int Ncells   = mesh->Ncells;
    int Mspecies = model->Mspecies;
    int Ndofs = Ncells*Mspecies;
    
    /* Boundary triangles */
    double *p = mesh->p;
    int *sd = model->sd;
    
    vector <plane> boundaries;
    
    int ntri=mesh->ntri;
    int i,j;
    
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
    
    for (i=0; i<Ndofs; i++) {
        /* Point in plane, the vertex itself. */
        //for (j=0; j<Mspecies; j++)
        temp.id = i;
        boundaries.push_back(temp);
    }
    
    /* Compute the point in the plane as the mean value of the vertex iteself and the
     centroid of the surrounding triangles. */
    for (i=0; i<Ncells; i++) {
        for (j=0; j<Mspecies; j++) {
            boundaries[i*Mspecies+j].p[0]=p[3*i];
            boundaries[i*Mspecies+j].p[1]=p[3*i+1];
            boundaries[i*Mspecies+j].p[2]=p[3*i+2];
        }
    }
    
    /* Loop over the boundary triangles and assemble their contribution to the normals */
    triangle *tri;
    double n[3];
    double vec[3];
    
    int sv;
    double radius;
    
    size_t *jcD = model->jcD;
    size_t *irD = model->irD;
    
    int dof,to_dof;
    size_t k;
    double y[3];
    double alpha;
    
    /* We will need to be able determine whether a dof is on the boundary/surface or not. */
    int *onboundary;
    onboundary = (int *)calloc(Ndofs,sizeof(int));
    
    for (i=0; i<ntri; i++) {
        tri = mesh->bnd[i];
        for (j=0; j<Mspecies; j++) {
            onboundary[tri->v1p*Mspecies+j]=1;
            onboundary[tri->v2p*Mspecies+j]=1;
            onboundary[tri->v3p*Mspecies+j]=1;
        }
    }
    
    double amax;
    
    for (i=0; i<ntri; i++) {
        
        tri = mesh->bnd[i];
        
        /* Compute the unit normal of that triangle.*/
        trinormal(tri,n);
        //        printf("n1=[%g %g %g];\n",alpha*n[0],alpha*n[1],alpha*n[2]);
        /* Need to check that all normals that we assemble point in a direction that
         the species may move. */
        
        radius = norm(tri->v1);
        a = dot(vec,n);
        //        printf("a=%g\n",a);
        
        /* Add that normal to accumulated mean normals of the vertices of the triangle.*/
        sv = tri->v1p;
        
        for (j=0; j<Mspecies; j++) {
            
            alpha = 1.0;
            dof = Mspecies*sv+j;
            
            /* Find the edge with the largest dot-product with the normal. */
            a=0.0;
            amax =0.0;
            
            /* If the point has neighbours in the interior */
            for (k=jcD[dof]; k<jcD[dof+1]; k++) {
                
                to_dof=irD[k];
                
                y[0]=p[3*(to_dof-j)/Mspecies]-p[3*(dof-j)/Mspecies];
                y[1]=p[3*(to_dof-j)/Mspecies+1]-p[3*(dof-j)/Mspecies+1];
                y[2]=p[3*(to_dof-j)/Mspecies+2]-p[3*(dof-j)/Mspecies+2];
                
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
            //            printf("n2=[%g %g %g];\n",alpha*n[0],alpha*n[1],alpha*n[2]);
            boundaries[Mspecies*sv+j].n[0]+=alpha*n[0];
            boundaries[Mspecies*sv+j].n[1]+=alpha*n[1];
            boundaries[Mspecies*sv+j].n[2]+=alpha*n[2];
            boundaries[Mspecies*sv+j].isbnd=1;
        }
        
        sv = tri->v2p;
        for (j=0; j<Mspecies; j++) {
            
            alpha = 1.0;
            dof = Mspecies*sv+j;
            
            /* Find the edge with the largest dot-product with the normal. */
            a=0.0;
            amax =0.0;
            for (k=jcD[dof]; k<jcD[dof+1]; k++) {
                
                to_dof=irD[k];
                y[0]=p[3*(to_dof-j)/Mspecies]-p[3*(dof-j)/Mspecies];
                y[1]=p[3*(to_dof-j)/Mspecies+1]-p[3*(dof-j)/Mspecies+1];
                y[2]=p[3*(to_dof-j)/Mspecies+2]-p[3*(dof-j)/Mspecies+2];
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
        }
        
        sv = tri->v3p;
        for (j=0; j<Mspecies; j++) {
            
            /* Find the edge with the largest dot-product with the normal. */
            alpha = 1.0;
            dof = sv*Mspecies+j;
            a=0.0;
            amax =0.0;
            
            for (k=jcD[dof]; k<jcD[dof+1]; k++) {
                
                to_dof=irD[k];
                y[0]=p[3*(to_dof-j)/Mspecies]-p[3*(dof-j)/Mspecies];
                y[1]=p[3*(to_dof-j)/Mspecies+1]-p[3*(dof-j)/Mspecies+1];
                y[2]=p[3*(to_dof-j)/Mspecies+2]-p[3*(dof-j)/Mspecies+2];
                
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
            //			printf("n4=[%g %g %g];\n",alpha*n[0],alpha*n[1],alpha*n[2]);
            boundaries[Mspecies*sv+j].n[0]+=alpha*n[0];
            boundaries[Mspecies*sv+j].n[1]+=alpha*n[1];
            boundaries[Mspecies*sv+j].n[2]+=alpha*n[2];
            boundaries[Mspecies*sv+j].isbnd=1;
            
            
        }
        if(boundaries[Mspecies*sv+j].isbnd==1)
        printf("n=[%g %g %g];\n",boundaries[Mspecies*sv+j].n[0],boundaries[Mspecies*sv+j].n[1],boundaries[Mspecies*sv+j].n[2]);
    }
    
    /* Normalize the mean normals */
    plane *tmp;
    for (i=0; i<Ndofs; i++) {
        
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
        
        for (j=0;j<Mspecies;j++){
            
            if (boundaries[i*Mspecies+j].isbnd){ 
                vec_in_plane(boundaries[i*Mspecies+j].n,boundaries[i*Mspecies+j].v1,boundaries[i*Mspecies+j].v2);
            }
        }	
        
        for (j=0; j<Mspecies; j++) {
            /* Fulkod */
            for(int l = 0;l<Mspecies;l++)
            {	
                (boundaries[i*Mspecies+j].type).push_back(0);
            }
        }	
        
    }
    
    free(onboundary);
    
    return boundaries;	
    
}


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

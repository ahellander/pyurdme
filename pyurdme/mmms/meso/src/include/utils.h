#ifndef UTILS_H
#define UTILS_H

double dist3(double *v1,double *v2){
    return sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2));
}

void print_pos(double *pos){
    printf("(%.5g,%.5g,%.5g)",pos[0],pos[1],pos[2]);
}

void print_id(int ID){
    printf("ID=%d",ID);
}


void print_group(group *grp){
    int M = (int)(grp->particles.size());

    printf("___GROUP ID=%d___\n",grp->unique_id);
    for(int i=0;i<M;i++){
        print_pos(grp->particles[i].pos);
        printf(", ");
        print_id(grp->particles[i].unique_id);
        printf("\n");

    }
    printf("_________________\n");
}

void print_neighbors(int *n,int N){
    printf("___NEIGHBORS VECTOR___\n");
    for(int i=0;i<N;i++){
        printf("(%d,%d)\n",i,n[i]);
    }
    printf("______________________\n");
}

void print_group_sizes(vector <group>& grps){
    int N = (int)(grps.size());
    int n1=0,n2=0;
    for(int i = 0;i<N;i++){
        if((int)(grps[i].particles.size())==1){
            n1++;
        }
        else if((int)(grps[i].particles.size())==2){
            n2++;
        }
    }
    printf("Groups with 1: %d, groups with 2: %d\n",n1,n2);
    printf("Relative: %.3g, %.3g\n",(double)n1/N,(double)n2/N);
}

void print_molecules(group *grp,double dt,FILE *f_out){
    int M = (int)(grp->particles.size());
    fprintf(f_out,"%d\n",M);
    for(int i=0;i<M;i++){

        fprintf(f_out,"%.5g %.5g %.5g %d %d %d %.5g %.5g\n",grp->particles[i].pos[0],grp->particles[i].pos[1],grp->particles[i].pos[2],grp->particles[i].type,grp->particles[i].unique_id,grp->particles[i].meso_micro,grp->particles[i].clock,dt);
    }
}

void print_molecules_voxel(group *grp,double dt,FILE *f_out){
    int M = (int)(grp->particles.size());
    fprintf(f_out,"%d\n",M);
    for(int i=0;i<M;i++){

        fprintf(f_out,"%d %d %d %.5g\n",grp->particles[i].vox_id,grp->particles[i].type,grp->particles[i].unique_id,dt);
    }
}

void print_molecule(group *grp,int ID,double dt,FILE *f_out){
    int M = (int)(grp->particles.size());
    for(int i=0;i<M;i++){
        if(grp->particles[i].unique_id==ID){
            fprintf(f_out,"%.8g %.8g %.8g %d %d %.8g\n",grp->particles[i].pos[0],grp->particles[i].pos[1],grp->particles[i].pos[2],grp->particles[i].type,grp->particles[i].unique_id,dt);
        }
    }
}

void print_tent_events(vector <tent_event>& t_events){
    int N = (int)(t_events.size());
    printf("___TENTATIVE EVENTS___\n");
    for(int i=0;i<N;i++){
        printf("event type: %d, time=%.5g\n",t_events[i].type,t_events[i].t);
    }
    printf("----------------------\n");
}

void print_vec(double *v){
    printf("(%.5g,%.5g,%.5g)",v[0],v[1],v[2]);
}

void print_vecs(particle *p){
    printf("vec1=");
    print_vec(p->vec1);
    printf("\n");
    printf("vec2=");
    print_vec(p->vec2);
    printf("\n");
    printf("vec3=");
    print_vec(p->vec3);
    printf("\n");
}


void cross(double *v1,double *v2,double *out){
    out[0] = v1[1]*v2[2]-v1[2]*v2[1];
    out[1] = v1[2]*v2[0]-v1[0]*v2[2];
    out[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

/* Gets the time step for a pair of molecules, taking into account the distance to the boundary. */
double get_boundary_dist_pair(double dt,vector <particle>& particles,vector <species>& specs,double *boundary, int dimension){
//    double dt_in = dt;
    double temp = INFINITY,temp2 = INFINITY;

    for(int j=0;j<2;++j){
        for(int i=0;i<dimension;++i){
            temp2 = min(fabs(particles[j].pos[i]-boundary[2*i]),fabs(particles[j].pos[i]-boundary[2*i+1]));
            if(temp2<temp){
                temp = temp2;
            }
        }
    }
    double sigma = specs[particles[0].type].sigma+specs[particles[1].type].sigma;
    double D = specs[particles[0].type].D+specs[particles[1].type].D;
    temp = pow(temp,2)/(50*D);

    if(temp<dt){
//        printf("temp=%g, dt=%g\n",temp,dt);
        /* TODO: 20 should be defined in a parameter struct. */
        double dp = pow(dist3(particles[0].pos,particles[1].pos)-sigma,2)/(30*D);

        if(dp<dt){
//            printf("dt_in=%g, dt_out=%g\n",dt_in,max(temp,dp));
            return max(temp,dp);
        }
        else{
            return dt;
        }

    }

    return dt;
}

void check_positions(group *grp){
    int M = (int)(grp->particles.size());


    for(int i=0;i<M;i++){
        if(grp->particles[i].pos[0] != grp->particles[i].pos[0] ||
           grp->particles[i].pos[1] != grp->particles[i].pos[1] ||
           grp->particles[i].pos[2] != grp->particles[i].pos[2]){
            printf("pos nan\n");
        }
    }
}

void generate_particles(group *grp,double *boundary,int N,int type,gsl_rng *rng,int *unique_id){

    int M = (int)(grp->particles.size());
    grp->particles.resize(M+N);

    double hx = boundary[1]-boundary[0];
    double hy = boundary[3]-boundary[2];
    double hz = boundary[5]-boundary[4];

    for(int i=M;i<M+N;i++){
        grp->particles[i].pos[0] = hx*gsl_rng_uniform(rng);
        grp->particles[i].pos[1] = hy*gsl_rng_uniform(rng);
        grp->particles[i].pos[2] = hz*gsl_rng_uniform(rng);
        grp->particles[i].voxel[0] = 0;
        grp->particles[i].voxel[1] = 0;
        grp->particles[i].voxel[2] = 0;
        grp->particles[i].type = type;
        grp->particles[i].active = true;
        grp->particles[i].unique_id = *unique_id;
        *unique_id = *unique_id+1;
        grp->particles[i].clock = 0.0;
        grp->particles[i].vec1[0] = 1.0;
        grp->particles[i].vec1[1] = 0.0;
        grp->particles[i].vec1[2] = 0.0;
        grp->particles[i].vec2[0] = 0.0;
        grp->particles[i].vec2[1] = 1.0;
        grp->particles[i].vec2[2] = 0.0;
        grp->particles[i].vec3[0] = 0.0;
        grp->particles[i].vec3[1] = 0.0;
        grp->particles[i].vec3[2] = 1.0;
    }
}


void generate_particles_meso(group *grp,vector <voxel>& voxels,int N,int type,gsl_rng *rng,int *unique_id){
    
    int M = (int)(grp->particles.size());
    grp->particles.resize(M+N);
    
    int Mvox = (int)(voxels.size());
    
    for(int i=M;i<M+N;i++){
        grp->particles[i].pos[0] = 0.0;
        grp->particles[i].pos[1] = 0.0;
        grp->particles[i].pos[2] = 0.0;
        
        /*TODO: Sample according to volume of voxels. */
        grp->particles[i].vox_id = floor(Mvox*gsl_rng_uniform(rng));

        
        grp->particles[i].type = type;
        grp->particles[i].active = true;
        grp->particles[i].unique_id = *unique_id;
        *unique_id = *unique_id+1;
        grp->particles[i].clock = 0.0;
        grp->particles[i].vec1[0] = 1.0;
        grp->particles[i].vec1[1] = 0.0;
        grp->particles[i].vec1[2] = 0.0;
        grp->particles[i].vec2[0] = 0.0;
        grp->particles[i].vec2[1] = 1.0;
        grp->particles[i].vec2[2] = 0.0;
        grp->particles[i].vec3[0] = 0.0;
        grp->particles[i].vec3[1] = 0.0;
        grp->particles[i].vec3[2] = 1.0;
    }
}



double dot(double *p1,double *p2){
    return p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
}

void vec_diff(double *p1,double *p2,double *outv){
    outv[0] = p2[0]-p1[0];
    outv[1] = p2[1]-p1[1];
    outv[2] = p2[2]-p1[2];
}


/* Rotates and scales a vector v in 3D around the axis z an angle theta.
 r is the new length of the vector. */
void rotation(double *v,double *z,double theta,double r)
{
	double normv = sqrt(dot(v,v));
	double normz = sqrt(dot(z,z));
	double vn[] = {v[0]/normv,v[1]/normv,v[2]/normv};
	double zn[] = {z[0]/normz,z[1]/normz,z[2]/normz};

	double q[4];
	q[0] = cos(theta/2.0);
	q[1] = sin(theta/2.0)*zn[0];
	q[2] = sin(theta/2.0)*zn[1];
	q[3] = sin(theta/2.0)*zn[2];

	double v_rot[3];
	v_rot[0] = (q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3])*vn[0]+2*(q[1]*q[2]-q[0]*q[3])*vn[1]+2*(q[1]*q[3]+q[0]*q[2])*vn[2];
	v_rot[1] = (q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3])*vn[1]+2*(q[1]*q[2]+q[0]*q[3])*vn[0]+2*(q[2]*q[3]-q[0]*q[1])*vn[2];
	v_rot[2] = (q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3])*vn[2]+2*(q[1]*q[3]-q[0]*q[2])*vn[0]+2*(q[2]*q[3]+q[0]*q[1])*vn[1];
	v[0] = r*v_rot[0];
	v[1] = r*v_rot[1];
	v[2] = r*v_rot[2];
}




void reflect_cuboid(group *grp,double *boundary,int dimension){
    int M = (int)(grp->particles.size());
    int counter = 0;
    for(int i=0;i<M;i++){
        for(int j=0;j<dimension;j++){
            while(grp->particles[i].pos[j]<boundary[2*j] || grp->particles[i].pos[j]>boundary[2*j+1]){
                if(grp->particles[i].pos[j]<boundary[2*j]){
                    grp->particles[i].pos[j] = 2*boundary[2*j]-grp->particles[i].pos[j];
                }
                else if(grp->particles[i].pos[j]>boundary[2*j+1]){
                    grp->particles[i].pos[j] = 2*boundary[2*j+1]-grp->particles[i].pos[j];
                }
                counter++;
                if(counter>50){
                    printf("Error: Molecule outside of domain. Cannot recover.\n");
                    print_group(grp);
                    exit(EXIT_FAILURE);

                }
            }
        }
    }
}

void reflect_cuboid_p(particle *p,double *boundary){
    int counter = 0;
    for(int j=0;j<3;j++){
        while(p->pos[j]<boundary[2*j] || p->pos[j]>boundary[2*j+1]){
            if(p->pos[j]<boundary[2*j]){
                p->pos[j] = 2*boundary[2*j]-p->pos[j];
            }
            else if(p->pos[j]>boundary[2*j+1]){
                p->pos[j] = 2*boundary[2*j+1]-p->pos[j];
            }
            counter++;
            if(counter>50){
                printf("Error: Molecule outside of domain. Cannot recover.\n");
                exit(EXIT_FAILURE);

            }
        }
    }
}

//void reflect_boundary(particle *p,vector <plane>& tet_b){
//    for(int i=0;i<(int)(tet_b.size());i++){
//        double dp = tet_b[i].n[0]*(p->pos[0]-tet_b[i].p[0])+tet_b[i].n[1]*(p->pos[1]-tet_b[i].p[1])+tet_b[i].n[2]*(p->pos[2]-tet_b[i].p[2]);
//        if(dp<0){
//            p->pos[0] = 2*fabs(dp)*tet_b[i].n[0];
//            p->pos[0] = 2*fabs(dp)*tet_b[i].n[1];
//            p->pos[0] = 2*fabs(dp)*tet_b[i].n[2];
//        }
//    }
//}

//int find_closest_node(particle *p,vector <node>& nodes,int space_or_boundary){
//    double td = INFINITY;
//    double td2;
//    int cn = -1;
//    for(int i=0;i<(int)(nodes.size());i++){
//        if(tet_b[i].type==space_or_boundary){
//            td2 = dist3(p->pos,tet_b[i].p);
//            if(td2<td){
//                td = td2;
//                cn = i;
//            }
//        }
//    }
//    return cn;
//}

void copy_pos(double *v,double *src){
    v[0] = src[0];
    v[1] = src[1];
    v[2] = src[2];
}

void copy_vecs(double *vec1,double *vec2,double *vec3,particle *p){
    /* TODO: Use mem_cpy instead. */
    vec1[0] = p->vec1[0];
    vec1[1] = p->vec1[1];
    vec1[2] = p->vec1[2];

    vec2[0] = p->vec2[0];
    vec2[1] = p->vec2[1];
    vec2[2] = p->vec2[2];

    vec3[0] = p->vec3[0];
    vec3[1] = p->vec3[1];
    vec3[2] = p->vec3[2];
}

void diffuse_vec(double *pos,double *vec1,double *vec2, double *vec3,double D,double dt,gsl_rng *rng){
    double sd = sqrt(2*D*dt);
    double r1 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    double r2 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    double r3 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    pos[0] += r1*vec1[0]+r2*vec2[0]+r3*vec3[0];
    pos[1] += r1*vec1[1]+r2*vec2[1]+r3*vec3[1];
    pos[2] += r1*vec1[2]+r2*vec2[2]+r3*vec3[2];
//    if(pos[2]!=0 && fabs(pos[2])>1e-10)
//        printf("pos_diff=(%g,%g,%g)\n",pos[0],pos[1],pos[2]);
}

void diffuse(particle *p,vector <species>& specs,double dt,gsl_rng *rng){
//    printf("diffusing molecule...\n");
    double D = specs[p->type].D;
    double sd = sqrt(2*D*dt);
    double r1,r2,r3;
//    print_vecs(p);
    r1 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    r2 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    r3 = sd*gsl_ran_gaussian_ziggurat(rng,1.0);
    p->pos[0] += r1*p->vec1[0]+r2*p->vec2[0]+r3*p->vec3[0];
    p->pos[1] += r1*p->vec1[1]+r2*p->vec2[1]+r3*p->vec3[1];
    p->pos[2] += r1*p->vec1[2]+r2*p->vec2[2]+r3*p->vec3[2];
//    for(int i=0;i<3;i++){
//        p->pos[i] += sd*gsl_ran_gaussian_ziggurat(rng,1.0);
//    }
}

void divide_into_cubes(group *grp,vector <group>& grps,int Nx,int Ny,int Nz,double *boundary){

    grps.resize(Nx*Ny*Nz);

    double hx = (boundary[1]-boundary[0])/Nx;
    double hy = (boundary[3]-boundary[2])/Ny;
    double hz = (boundary[5]-boundary[4])/Nz;

    int ix,iy,iz;

    int index;
    int M = (int)(grp->particles.size());
    double *pos;

    for(int i=0;i<M;i++){
        pos = grp->particles[i].pos;
        ix = (int)(pos[0]/hx);
        iy = (int)(pos[1]/hy);
        iz = (int)(pos[2]/hz);
//        if(dimension==3){
//
//        }
//        else{
//            iz = 0;
//        }
        index = ix+Nx*iy+Nx*Ny*iz;
        if(index>Nx*Ny*Nz-1){
//            printf("(%.5g,%.5g,%.5g)\n",grp->particles[i].pos[0],grp->particles[i].pos[1],grp->particles[i].pos[2]);
            index = Nx*Ny*Nz-1;
        }
        if(index<0){
//            printf("(%.5g,%.5g,%.5g)\n",grp->particles[i].pos[0],grp->particles[i].pos[1],grp->particles[i].pos[2]);

            index = 0;
        }
        grps[index].particles.push_back(grp->particles[i]);
    }
}

void sort_molecules(group *grp,vector <group>& grps,vector <species>& specs){
    int M = (int)(grp->particles.size());
    if(M>1){

        int *n1,*n2;
        n1 = (int *)malloc(sizeof(int)*M);
        n2 = (int *)malloc(sizeof(int)*M);
        double *d1,*d2;
        d1 = (double *)malloc(sizeof(double)*M);
        d2 = (double *)malloc(sizeof(double)*M);

        for(int i=0;i<M;i++){
            d1[i] = INFINITY;
            d2[i] = INFINITY;
        }

        double d_temp,d_temp2;
        int n_temp;
        double *pos_temp1;
        for(int i=0;i<M;i++){

            pos_temp1 = grp->particles[i].pos;

            for(int j=i+1;j<M;j++){
                double D = specs[grp->particles[j].type].D+specs[grp->particles[i].type].D;
                double sigma = specs[grp->particles[j].type].sigma+specs[grp->particles[i].type].sigma;
                /* TODO: 20 should be specified in a parameter struct. */
                d_temp = pow(dist3(pos_temp1,grp->particles[j].pos)-sigma,2)/(20*D);
                /* ::::::::::::::::::::::::::::::::::::: */
                /* ::::::::::::::::::::::::::::::::::::: */
                /*      TODO: Make function of this.     */
                /* ::::::::::::::::::::::::::::::::::::: */
                if(d_temp<d1[i]){
                    d_temp2 = d1[i];

                    d1[i] = d_temp;//pow(d_temp-sigma,2)/(50*D);
                    d2[i] = d_temp2;
                    n_temp = n1[i];
                    n1[i] = j;
                    n2[i] = n_temp;
                }
                else if(d_temp<d2[i]){
                    d2[i] = d_temp;
                    n2[i] = j;
                }
                /* ::::::::::::::::::::::::::::::::::::: */
                /* ::::::::::::::::::::::::::::::::::::: */
                /* ::::::::::::::::::::::::::::::::::::: */
                if(d_temp<d1[j]){
                    d_temp2 = d1[j];
                    d1[j] = d_temp;
                    d2[j] = d_temp2;
                    n_temp = n1[j];
                    n1[j] = i;
                    n2[j] = n_temp;
                }
                else if(d_temp<d2[j]){
                    d2[j] = d_temp;
                    n2[j] = i;
                }
            }
        }
        group grp_temp;
        for(int i=0;i<M;i++){
            /* Molecules are each other's nearest neighbors. */
            if(i<n1[i] && i==n1[n1[i]] && grp->particles[i].active && grp->particles[n1[i]].active){
                grp_temp.particles.resize(2);
                grp_temp.particles[0] = grp->particles[i];
                grp_temp.particles[1] = grp->particles[n1[i]];
                if(d1[i]!=d1[n1[i]]){
                    printf("Sorted as nearest neighbors, but aren't.\n");
                    printf("pos1: (%g,%g,%g)\n",grp_temp.particles[0].pos[0],grp_temp.particles[0].pos[1],grp_temp.particles[0].pos[2]);
                    printf("pos2: (%g,%g,%g)\n",grp_temp.particles[1].pos[0],grp_temp.particles[1].pos[1],grp_temp.particles[1].pos[2]);
                    exit(EXIT_FAILURE);
                }
                grp_temp.dist = min(d2[i],d2[n1[i]]);
                grps.push_back(grp_temp);
                grp->particles[i].active = false;
                grp->particles[n1[i]].active = false;
            }
            /* Molecule will be updated as a single molecules. */
            else if(grp->particles[i].active){
                grp_temp.particles.resize(1);
                grp_temp.particles[0] = grp->particles[i];
                grp_temp.dist = d1[i];
                grps.push_back(grp_temp);
                grp->particles[i].active = false;
            }

        }
        free(n1);
        free(n2);
        free(d1);
        free(d2);
    }
    else if(M==1){
        group grp_temp;
        grp_temp.particles.resize(1);
        grp_temp.particles[0] = grp->particles[0];
        grp_temp.dist = INFINITY;
        grps.push_back(grp_temp);
        grp->particles[0].active = false;
    }
}

int cp(const char *to, const char *from)
{
    int fd_to, fd_from;
    char buf[4096];
    ssize_t nread;
    int saved_errno;

    fd_from = open(from, O_RDONLY);
    if (fd_from < 0)
        return -1;

    fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd_to < 0)
        goto out_error;

    while (nread = read(fd_from, buf, sizeof buf), nread > 0)
    {
        char *out_ptr = buf;
        ssize_t nwritten;

        do {
            nwritten = write(fd_to, out_ptr, nread);

            if (nwritten >= 0)
            {
                nread -= nwritten;
                out_ptr += nwritten;
            }
            else if (errno != EINTR)
            {
                goto out_error;
            }
        } while (nread > 0);
    }

    if (nread == 0)
    {
        if (close(fd_to) < 0)
        {
            fd_to = -1;
            goto out_error;
        }
        close(fd_from);

        /* Success! */
        return 0;
    }

out_error:
    saved_errno = errno;

    close(fd_from);
    if (fd_to >= 0)
        close(fd_to);

    errno = saved_errno;
    return -1;
}


#endif

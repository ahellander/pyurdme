#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "time.h"
#include <sys/time.h>
#include <sys/stat.h>

/* Gnu scientific library */
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_errno.h"

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#include "global_params.h"
#include "structs.h"
#include "model_parser.h"

#include "utils.h"
#include "rand_gen.h"

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>
#include "../hdf5-1.8.16/c++/src/H5Cpp.h"
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "mesh.h"
#include "coupling.h"

#include "urdmemodel.h"
//#include "read_model.h"
#include "read_matfile.h"

int dimension;
int nxsteps = 20;
int ntsteps = 30;

using namespace std;

/* TODO: In 2D: molecules at distance less than sigma. */

void main_simulator(group *grp,vector <species>& specs,vector <association>& assocs,vector <dissociation>& dissocs,vector <plane>& boundaries,double T,gsl_rng *rng,int traj_num);

bool compare_events (tent_event e1,tent_event e2) {
    return e1.t<e2.t;
}

double check_association(vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <association> associations,gsl_rng *rng,double t,int i,int j,double dt,double *micro_param){
    
    /* r_new is only used in the 2D case. In 3D it doesn't mean anything at this point. */
    double r_new = 0.0;
    
    double r_0 = micro_param[0];
    
    double D = micro_param[1];
    double sigma = micro_param[2];
    if(r_0<sigma){
        r_0 = 1.001*sigma;
    }
    int M = (int)(associations.size());
    tent_event temp_event;
    temp_event.type = 1;
    temp_event.reactants.resize(2);
    int m = 0;
    for(int k=0;k<M;k++){
        if((associations[k].reactant1==particles[i].type && associations[k].reactant2==particles[j].type) || (associations[k].reactant2==particles[i].type && associations[k].reactant1==particles[j].type)){
            
            
            m++;
            temp_event.index = k;
            temp_event.reactants[0] = i;
            temp_event.reactants[1] = j;
            micro_param[3] = associations[k].k;
            
            if(dimension==3){
                temp_event.t = random_time(dt,r_0,associations[k].k,sigma,D,rng,dimension);
            }
            else if(dimension==2){

                
                /* TODO: What to if multiple possible reactions? */
                r_new = random_time_num(associations[k].k,r_0,sigma,D,dt,rng,&temp_event.t,dimension,ntsteps,nxsteps,gsl_rng_uniform(rng));
                
            }
            
            tent_events.insert(lower_bound(tent_events.begin(),tent_events.end(),temp_event,compare_events),temp_event);
            
        }
    }
    if(dimension==2 && m==0){
        r_new = random_time_num(1e-20,r_0,sigma,D,dt,rng,&temp_event.t,dimension,ntsteps,nxsteps,gsl_rng_uniform(rng));
        int i = 0;
        while(r_new<sigma  && i<20){
            r_new = random_time_num(1e-20,r_0,sigma,D,dt,rng,&temp_event.t,dimension,ntsteps,nxsteps,gsl_rng_uniform(rng));
            i++;
        }
        if(i==20){
            printf("Unable to sample new distance between molecules.\n");
            printf("r_new=%g, sigma=%g, D=%g, r_0=%g\n",r_new,sigma,D,r_0);
            exit(EXIT_FAILURE);
        }
    }
    return r_new;
}

void update_diss(vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <dissociation> dissociations,gsl_rng *rng,double t,int i){
    
    int M = (int)(dissociations.size());
    tent_event temp_event;
    temp_event.type = 0;
    
    for(int j = 0;j<M;j++){
        if(dissociations[j].reactant==particles[i].type){
            temp_event.index = j;
            temp_event.t = -log(gsl_rng_uniform(rng))/dissociations[j].k;//+t;
            temp_event.reactants.push_back(i);
            tent_events.insert(lower_bound(tent_events.begin(),tent_events.end(),temp_event,compare_events),temp_event);
        }
    }
}


void doAssociation(vector <particle>& particles,vector <species>& specs,tent_event *event,vector <association>& assocs,gsl_rng *rng){
    double *pos0,*pos1;
    double abs_pos[3];
    /* :::::::::::::::::::::::::::::::::::::::::*/
    /* :::::::::::::::::::::::::::::::::::::::::*/
    /* TODO: Make function of this. */
    /* :::::::::::::::::::::::::::::::::::::::::*/
    pos0 = particles[0].pos;
    pos1 = particles[1].pos;
    double D0 = specs[particles[0].type].D;
    double D1 = specs[particles[1].type].D;
    double sigma0 = specs[particles[0].type].sigma;
    double sigma1 = specs[particles[1].type].sigma;
    double d0 = sqrt(D1/D0);
    double d1 = 1.0/d0;

    abs_pos[0] = d0*pos0[0]+d1*pos1[0];
    abs_pos[1] = d0*pos0[1]+d1*pos1[1];
    abs_pos[2] = d0*pos0[2]+d1*pos1[2];
    
    /* :::::::::::::::::::::::::::::::::::::::::*/
    /* :::::::::::::::::::::::::::::::::::::::::*/
    /* :::::::::::::::::::::::::::::::::::::::::*/
    
    /* We are assuming here that the molecules live on the same plane if the reaction is 2D. */

    diffuse_vec(abs_pos,particles[0].vec1,particles[0].vec2,particles[0].vec3,D0+D1,sigma0+sigma1,event->t,rng);
    double vec1[3],vec2[3],vec3[3];
    copy_vecs(vec1,vec2,vec3,&particles[0]);

    int M = (int)(assocs[event->index].products.size());
    particles.resize(M);
    if(M>0){

        /* TODO: CORRECT THIS! */
        particles[0].pos[0] = (d1*abs_pos[0])/(1+d1*d1);
        particles[0].pos[1] = (d1*abs_pos[1])/(1+d1*d1);
        particles[0].pos[2] = (d1*abs_pos[2])/(1+d1*d1);
        
        particles[0].type = assocs[event->index].products[0];
        particles[0].active = true;
        
        
        particles[0].vec1[0] = vec1[0];
        particles[0].vec1[1] = vec1[1];
        particles[0].vec1[2] = vec1[2];
        
        particles[0].vec2[0] = vec2[0];
        particles[0].vec2[1] = vec2[1];
        particles[0].vec2[2] = vec2[2];
        
        particles[0].vec3[0] = vec3[0];
        particles[0].vec3[1] = vec3[1];
        particles[0].vec3[2] = vec3[2];
        
    }
    

    
    /* TODO: This does not work in 2D yet. */
    for(int i=1;i<M;i++){
        double sigma = specs[particles[0].type].sigma+specs[assocs[event->index].products[i]].sigma;
        double new_pos[3];
        randomSphere(rng,sigma,new_pos,dimension);
        new_pos[0] += particles[0].pos[0];
        new_pos[1] += particles[0].pos[1];
        new_pos[2] += particles[0].pos[2];
        particles[i].pos[0] = new_pos[0];
        particles[i].pos[1] = new_pos[1];
        particles[i].pos[2] = new_pos[2];
        particles[i].type = assocs[event->index].products[i];
        particles[i].active = true;
      
    }
    
}

void doDissociation(vector <particle>& particles,vector <species>& specs,tent_event *event,vector <dissociation>& dissocs,vector <plane>& boundaries,gsl_rng *rng){

    int M = (int)(dissocs[event->index].products.size());
   
    double vec1[3],vec2[3],vec3[3];
    copy_vecs(vec1,vec2,vec3,&particles[0]);
    double pos[3];
    int reactant = event->reactants[0];
    pos[0] = particles[reactant].pos[0];
    pos[1] = particles[reactant].pos[1];
    pos[2] = particles[reactant].pos[2];
    particles.erase(particles.begin()+reactant);
    int N = (int)(particles.size());
    particles.resize(N+M);
    
    int index = event->index;
    int index_temp;

    
    if(M>0){
        if(N>0){
            index_temp = 1;
        }
        else{
            index_temp = 0;
        }
        particles[index_temp].pos[0] = pos[0];
        particles[index_temp].pos[1] = pos[1];
        particles[index_temp].pos[2] = pos[2];
        
//        reflect_cuboid_p(&particles[index_temp],boundary);

        particles[index_temp].type = dissocs[index].products[0];
        particles[index_temp].active = true;
        
        
        particles[index_temp].vec1[0] = vec1[0];
        particles[index_temp].vec1[1] = vec1[1];
        particles[index_temp].vec1[2] = vec1[2];
        
        particles[index_temp].vec2[0] = vec2[0];
        particles[index_temp].vec2[1] = vec2[1];
        particles[index_temp].vec2[2] = vec2[2];
        
        particles[index_temp].vec3[0] = vec3[0];
        particles[index_temp].vec3[1] = vec3[1];
        particles[index_temp].vec3[2] = vec3[2];
        
        
        /* Copy vec1,vec2,vec3. */
    }
    double new_pos[3];

    for(int i=1;i<M;i++){
        
        double sigma = specs[particles[index_temp].type].sigma+specs[dissocs[event->index].products[i]].sigma;
        
        /* TODO: This part should be a function, as we're repeating the code below. */
        randomSphere(rng,sigma,new_pos,dimension);

        new_pos[0] += particles[index_temp].pos[0];
        new_pos[1] += particles[index_temp].pos[1];
        new_pos[2] += particles[index_temp].pos[2];
        
        /* TODO: Verify that this projection is correct. */
        particles[index_temp+i].pos[0] = vec1[0]*new_pos[0]+vec2[0]*new_pos[1]+vec3[0]*new_pos[2];
        particles[index_temp+i].pos[1] = vec1[1]*new_pos[0]+vec2[1]*new_pos[1]+vec3[1]*new_pos[2];
        particles[index_temp+i].pos[2] = vec1[2]*new_pos[0]+vec2[2]*new_pos[1]+vec3[2]*new_pos[2];
        
        /* This is to make sure that the new particles don't overlap after being reflected back into the domain. */
//        reflect_cuboid_p(&particles[index_temp+i],boundary);
        while(dist3(particles[index_temp].pos,particles[index_temp+i].pos)<0.9999*sigma || (index_temp>=1 && dist3(particles[0].pos,particles[index_temp+i].pos)<0.9999*sigma)){
            
            randomSphere(rng,sigma,new_pos,dimension);
            
            new_pos[0] += particles[index_temp].pos[0];
            new_pos[1] += particles[index_temp].pos[1];
            new_pos[2] += particles[index_temp].pos[2];
            
            /* TODO: Verify that this projection is correct. */
            particles[index_temp+i].pos[0] = vec1[0]*new_pos[0]+vec2[0]*new_pos[1]+vec3[0]*new_pos[2];
            particles[index_temp+i].pos[1] = vec1[1]*new_pos[0]+vec2[1]*new_pos[1]+vec3[1]*new_pos[2];
            particles[index_temp+i].pos[2] = vec1[2]*new_pos[0]+vec2[2]*new_pos[1]+vec3[2]*new_pos[2];
            
//            reflect_cuboid_p(&particles[index_temp+i],boundary);
        }
        
        
        particles[index_temp+i].type = dissocs[event->index].products[i];
        particles[index_temp+i].active = true;
        
        
        particles[index_temp+i].vec1[0] = vec1[0];
        particles[index_temp+i].vec1[1] = vec1[1];
        particles[index_temp+i].vec1[2] = vec1[2];
        
        particles[index_temp+i].vec2[0] = vec2[0];
        particles[index_temp+i].vec2[1] = vec2[1];
        particles[index_temp+i].vec2[2] = vec2[2];
        
        particles[index_temp+i].vec3[0] = vec3[0];
        particles[index_temp+i].vec3[1] = vec3[1];
        particles[index_temp+i].vec3[2] = vec3[2];
    }


}

void simulate_group(group *grp,vector <species>& specs,double T,vector <association>& assocs,vector <dissociation>& dissocs,vector <plane>& boundaries,gsl_rng *rng,int traj_num){
    
    
    
    vector <tent_event> t_events;
    
    double t_loc = 0.0;
    double dt;
    double micro_param[4];

    while(t_loc<T){

        
        dt = T-t_loc;
        t_events.resize(0);
        
        int M = (int)(grp->particles.size());
        

        
        /* One molecule. We can diffuse or dissociate. */
        if(M==1){
            
            update_diss(grp->particles,specs,t_events,dissocs,rng,0.0,0);
            double t_reac = INFINITY;
            if((int)(t_events.size())>0){
                t_reac = t_events[0].t;
            }
            if(t_reac<dt){
                dt = t_reac;
                diffuse(&(grp->particles[0]),specs,dt,rng);
                doDissociation(grp->particles,specs,&t_events[0],dissocs,boundaries,rng);
            }
            else{
                diffuse(&(grp->particles[0]),specs,dt,rng);
            }
        }
        /* Two molecules. We can diffuse, dissociate, or associate. */
        else if(M==2){
            
            micro_param[0] = dist3(grp->particles[0].pos,grp->particles[1].pos);//r_0
            double r_new = micro_param[0];
            micro_param[1] = specs[grp->particles[0].type].D+specs[grp->particles[1].type].D;//D
            micro_param[2] = specs[grp->particles[0].type].sigma+specs[grp->particles[1].type].sigma;//sigma
            micro_param[3] = 0.0;
            
            
            double cutoff = 1.1;
            
            
            
            /* TODO: Time step should be limited by the distance to the boundary. */
//            dt = get_boundary_dist_pair(dt,grp->particles,specs,boundary,dimension);
            
            
            
            if(micro_param[0]<cutoff*micro_param[2]){
                if(dimension==2){
                    double max_points = 1e4;
                    double K_loc = 30;
                    double loc_rs = micro_param[0]-micro_param[2];
                    if(loc_rs<=0){
                        loc_rs = 0.0001*micro_param[2];
                    }
                    double max_dt = 1.0/(2*micro_param[1]*K_loc)*pow((max_points/nxsteps-1)*loc_rs,2);
                    if(max_dt<dt){
                        dt = max_dt;
                    }
                }
                
                
                
                /* Check for possible reactions. Can be associations or dissociations. */
                /* Check for associations. */
                r_new = check_association(grp->particles,specs,t_events,assocs,rng,0.0,0,1,dt,micro_param);
                
                /* Check for dissociations. */
                update_diss(grp->particles,specs,t_events,dissocs,rng,0.0,0);
                update_diss(grp->particles,specs,t_events,dissocs,rng,0.0,1);
                
                
                
                double t_reac = INFINITY;
                if((int)(t_events.size())>0){
                    t_reac = t_events[0].t;
                }
                
                /* Reaction fires. */
                if(t_reac<dt){
                    dt = t_reac;
                    /* Event is an association. */
                    if(t_events[0].type == 1){

                        doAssociation(grp->particles,specs,&t_events[0],assocs,rng);

                    }
                    /* Event is a dissociation. */
                    else if(t_events[0].type == 0){
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /*      TODO: Molecules should be updated as a pair */
                        /*            until one of them dissociated.        */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */

                        diffuse(&(grp->particles[0]),specs,dt,rng);
                        diffuse(&(grp->particles[1]),specs,dt,rng);
                        
                        doDissociation(grp->particles,specs,&t_events[0],dissocs,boundaries,rng);

                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                    }
                }
                /* No reaction. Molecules diffuse. */
                else{
                    double endr = micro_param[0]+40*sqrt(6*micro_param[1]*dt);
                    if(dimension==3){
                        r_new = random_r_irr(endr,dt,micro_param[0],micro_param[3],micro_param[2],micro_param[1],rng);

                        
                    }
                    if(r_new<micro_param[2]){
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /*      TODO: This is a hack, and this problem      */
                        /*            should be solved differently.         */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        r_new = micro_param[2]*1.001;
                    }
                    
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    /*                      Update pair of molecules.                     */
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    
                    double relativePosition[3];
                    vec_diff(grp->particles[0].pos,grp->particles[1].pos,relativePosition);
                    
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /*      TODO: Make function of this.     */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    const int N = 10;
                    double theta=0.0;
                    double psi=0.0;
                    
                    if(dimension==3){
                        
                        /* TODO: min dt should be fetched from a parameter struct. */
                        if(dt<1e-11){
                            dt = 1e-11;
                            theta = 0.0;
                        }
                        theta = gen_rnd_theta(micro_param[1],dt,micro_param[0],N,rng);
                    }
                    else if(dimension==2){
                        theta = gsl_ran_gaussian_ziggurat(rng,sqrt(2*micro_param[1]/pow(micro_param[0],2)*dt));
                        
                    }
                    
                    
                    if(dimension==3){
                        psi = 2*rand_gen_constant_pi*gsl_rng_uniform(rng);
                        double temp2[] = {relativePosition[0],relativePosition[1],relativePosition[2]};
                        double rotax[] = {relativePosition[1],-relativePosition[0],0.0};
                        rotation(relativePosition,rotax,theta,r_new);
                        rotation(relativePosition,temp2,psi,r_new);
                    }
                    else if(dimension==2){
                        /* TODO: This is a temporary hack. DOES NOT WORK FOR GENERAL PLANES. */
                        double rotax[] = {0.0,0.0,1.0};
                        rotation(relativePosition,rotax,theta,r_new);
                    }
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    
                    
                    
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /*      TODO: Make function of this.     */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    double abs_pos[3];
                    double *pos0;
                    double *pos1;
                    pos0 = grp->particles[0].pos;
                    pos1 = grp->particles[1].pos;
                    double D0 = specs[grp->particles[0].type].D;
                    double D1 = specs[grp->particles[1].type].D;
                    
                    double d0 = sqrt(D1/D0);
                    double d1 = 1.0/d0;
                    
                    abs_pos[0] = d0*pos0[0]+d1*pos1[0];
                    abs_pos[1] = d0*pos0[1]+d1*pos1[1];
                    abs_pos[2] = d0*pos0[2]+d1*pos1[2];

                    diffuse_vec(abs_pos,grp->particles[0].vec1,grp->particles[0].vec2,grp->particles[0].vec3,micro_param[1],micro_param[2],dt,rng);

                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */
                    /* ::::::::::::::::::::::::::::::::::::: */


                    
                    grp->particles[0].pos[0] = (d0*abs_pos[0]-relativePosition[0])/(1+d0*d0);
                    grp->particles[0].pos[1] = (d0*abs_pos[1]-relativePosition[1])/(1+d0*d0);
                    grp->particles[0].pos[2] = (d0*abs_pos[2]-relativePosition[2])/(1+d0*d0);
                    
                    grp->particles[1].pos[0] = grp->particles[0].pos[0]+relativePosition[0];
                    grp->particles[1].pos[1] = grp->particles[0].pos[1]+relativePosition[1];
                    grp->particles[1].pos[2] = grp->particles[0].pos[2]+relativePosition[2];
                    

                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                    /*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
                }
                

            }
            else{
                double dt_prev = dt;
                dt = pow(micro_param[0]-micro_param[2],2)/(50*micro_param[1]);
                if(dt>dt_prev){
                    dt = dt_prev;
                }
                update_diss(grp->particles,specs,t_events,dissocs,rng,0.0,0);
                update_diss(grp->particles,specs,t_events,dissocs,rng,0.0,1);
                
                
                
                double t_reac = INFINITY;
                if((int)(t_events.size())>0){
                    t_reac = t_events[0].t;
                }

                /* Reaction fires. */
                if(t_reac<dt){
                    dt = t_reac;
                    /* Event is a dissociation. */
                    if(t_events[0].type == 0){
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /*      TODO: Molecules should be updated as a pair */
                        /*            until one of them dissociated.        */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        //                        print_tent_events(t_events);
                        diffuse(&(grp->particles[0]),specs,dt,rng);
                        diffuse(&(grp->particles[1]),specs,dt,rng);
                        doDissociation(grp->particles,specs,&t_events[0],dissocs,boundaries,rng);

                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                        /* :::::::::::::::::::::::::::::::::::::::::::::::: */
                    }
                }
                else{
                    diffuse(&(grp->particles[0]),specs,dt,rng);
                    diffuse(&(grp->particles[1]),specs,dt,rng);
                }
                if(grp->particles[0].pos[0]!=grp->particles[0].pos[0]){
                    printf("grp->particles[0].pos[0]=%g, mol not close2.\n",grp->particles[0].pos[0]);
                }
                if(grp->particles[1].pos[0]!=grp->particles[1].pos[0]){
                    printf("grp->particles[0].pos[1]=%g, mol not close2.\n",grp->particles[0].pos[0]);
                }
            }
        }
        /* Maximum of two molecules in a group. */
        else if(M>2){
            main_simulator(grp,specs,assocs,dissocs,boundaries,dt,rng,traj_num);
        }
        

        
//        reflect_cuboid(grp,boundary,dimension);

        t_loc += dt;
    }
    
    /* TODO: Should restart pair if distance is less than reaction radius at this point. */
    
}


void main_simulator(group *grp,vector <species>& specs,vector <association>& assocs,vector <dissociation>& dissocs,vector <plane>& boundaries,double T,gsl_rng *rng,int traj_num){
   
    /*TODO: Re-run the reversible reaction test with Nx=Ny=Nz=1. */
    
    int Nx=1,Ny=1,Nz=1;
    
    
    double t = 0.0;
    
    vector <group> grps_final,grps;
    
    double dt;// = T;
    double dt_temp;
    
    
    
    while(t<T){
        dt = T-t;
        
        
        /* Does not work with general geometries. */
        //        divide_into_cubes(grp,grps,Nx,Ny,Nz,boundary);
        /* All molecules in grps[0]. */
        
        
        /* TODO: We should sort by voxel id instead. */
        grps.resize(0);
        grps.resize(1);
        for(int l=0;l<(int)(grp->particles.size());l++){
            grps[0].particles.push_back(grp->particles[l]);
        }
        
        
        grps_final.resize(0);
        
        
        
        
        
        dt_temp = INFINITY;
        
        for(int i=0;i<Nx*Ny*Nz;i++){
            sort_molecules(&grps[i],grps_final,specs);
        }
        
        int num_groups = (int)(grps_final.size());
        for(int i=0;i<num_groups;i++){
            if(grps_final[i].dist<dt_temp){
                dt_temp = grps_final[i].dist;
            }
        }
        if(dt_temp<dt){
            dt = dt_temp;
        }
//        printf("dt=%g\n",dt);
        
        if(dt<1e-12){
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /*      TODO: Minimum time step should be part of a  */
            /*            global parameter struct.               */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            dt = 1e-12;
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            
        }
        if(dt>1e-2){
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /*      TODO: Maximum time step should be part of a  */
            /*            global parameter struct.               */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            dt = 1e-2;
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
            /* ::::::::::::::::::::::::::::::::::::::::::::::::: */
        }
        int N_groups = (int)(grps_final.size());
        for(int i=0;i<N_groups;i++){
//            printf("t= %g\n",t);
            
            simulate_group(&grps_final[i],specs,dt,assocs,dissocs,boundaries,rng,traj_num);
            if(dimension==2){
                for(int j=0;j<(int)(grps_final[i].particles.size());j++){
                    grps_final[i].particles[j].pos[2] = 0.0;
                }
            }
        }
        
        
        
        
        
        grp->particles.resize(0);
        N_groups = (int)(grps_final.size());
        int grp_size;
        for(int i=0;i<N_groups;i++){
            grp_size = (int)(grps_final[i].particles.size());
            for(int j=0;j<grp_size;j++){
                if(grps_final[i].particles[j].active){
                    grp->particles.push_back(grps_final[i].particles[j]);
                }
            }
        }
        
        reflect_boundary(grp->particles,boundaries);
        
//        int Psize = (int)(grp->particles.size());
//        printf("[");
//        for(int q=0;q<Psize;q++){
//            printf("%g %g %g;\n",grp->particles[q].pos[0],grp->particles[q].pos[1],grp->particles[q].pos[2]);
//        }
//        printf("];");
//        printf("num_particles=%d, t=%g\n",Psize,t+dt);
        t += dt;
    }

}

void add_time_point_to_file(H5File file, group grp, int num_specs, int traj_num, int time_index){

    
    hsize_t  dims[2];        
    double position_buffer[10000][3];
    int ids_buffer[10000];
    int num_spec_type;


    for (int spec=0;spec<num_specs;spec++){

        num_spec_type=0;
        for(int j=0;j<(int)(grp.particles.size());j++){

            int type = grp.particles[j].type;

            if (type == spec){ 
               num_spec_type++; 
               position_buffer[j][0] = grp.particles[j].pos[0];
               position_buffer[j][1] = grp.particles[j].pos[1];
               position_buffer[j][2] = grp.particles[j].pos[2];
               ids_buffer[j] = grp.particles[j].unique_id;
            }

        }

        dims[0] = num_spec_type;
        dims[1] = 3;
        DataSpace position_dataspace = DataSpace(2, dims); 
        dims[1]=1;
        DataSpace ids_dataspace = DataSpace(2, dims); 

        DataSet positions = file.createDataSet("/Trajectories/"+to_string(traj_num)+"/Type_"+to_string(spec)+"/positions_"+to_string(time_index), PredType::NATIVE_DOUBLE,position_dataspace);
        DataSet ids = file.createDataSet("/Trajectories/"+to_string(traj_num)+"/Type_"+to_string(spec)+"/unique_ids_"+to_string(time_index), PredType::NATIVE_INT,ids_dataspace);

        positions.write(position_buffer, PredType::NATIVE_DOUBLE);
        ids.write(ids_buffer, PredType::NATIVE_INT);

    }
}

void read_p(fem_mesh *mesh, H5File mesh_file){

    DataSet dataset = mesh_file.openDataSet("/mesh/p");
    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    cout << "rank " << rank << ", dimensions " <<
          (unsigned long)(dims_out[0]) << " x " <<
          (unsigned long)(dims_out[1]) << endl;

    double *p;
    p=(double *)malloc(dims_out[0]*dims_out[1]*sizeof(double));
    dataset.read(p,PredType::NATIVE_DOUBLE);

    mesh->Ncells=dims_out[0];
    mesh->p=p;     

    // initialize vertices
    int nvox = dims_out[0];
    vertex **vertices = (vertex **)malloc(nvox*sizeof(vertex *));
//    vertex *vtx;
    for (int i=0; i<nvox; i++){
    
        vertices[i] = (vertex *)malloc(sizeof(vertex));
  //      vtx = vertices[i];
        //vtx->id = i;
	vertices[i]->id = i;
    }  
    mesh->vertices = vertices;
    
}

void read_t(fem_mesh *mesh, H5File mesh_file){

    DataSet dataset = mesh_file.openDataSet("/mesh/t");
    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    cout << "rank " << rank << ", dimensions " <<
          (unsigned long)(dims_out[0]) << " x " <<
          (unsigned long)(dims_out[1]) << endl;

    int *t;
    t=(int *)malloc(dims_out[0]*dims_out[1]*sizeof(int));
    dataset.read(t,PredType::NATIVE_INT);

    mesh->ntet = dims_out[0];
    mesh->t=t;     
    
}

void read_bnd(fem_mesh *mesh, H5File mesh_file){

    DataSet dataset = mesh_file.openDataSet("/mesh/boundaryfacets");
    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    cout << "rank " << rank << ", dimensions " <<
          (unsigned long)(dims_out[0]) << " x " <<
          (unsigned long)(dims_out[1]) << endl;

    int *t;
    t=(int *)malloc(dims_out[0]*dims_out[1]*sizeof(int));
    dataset.read(t,PredType::NATIVE_INT);

    /*for (int i=0;i<dims_out[0];i++){
        for (int j=0;j<dims_out[1];j++)
            cout << p[i+j*dims_out[0]] << " ";
        cout << "\n";
    }*/
    mesh->ntri = dims_out[0];
    mesh->e=t;     
    
}

void read_vertex_to_cell(fem_mesh *mesh, H5File mesh_file){
    DataSet dataset = mesh_file.openDataSet("/mesh/vertex2cells");
    DataSpace dataspace = dataset.getSpace();
    /*
    * Get the number of dimensions in the dataspace.
    */
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims_out[2];
    dataspace.getSimpleExtentDims( dims_out, NULL);
    cout << "rank " << rank << ", dimensions " <<
          (unsigned long)(dims_out[0]) << " x " <<
          (unsigned long)(dims_out[1]) << endl;



     int ttemp[(int)dims_out[0]][(int)dims_out[1]];
	
	for(int i=0;i<(int)dims_out[0];i++)
	{
		for(int j=0;j<(int)dims_out[1];j++)
		{
			ttemp[i][j] = 0.0;
		}
	}


     dataset.read(ttemp,PredType::NATIVE_INT);
  
//    vertex *vtx;
    for (int i=0; i<(int)dims_out[0];i++){
//	mesh->vertices[i]->cells.resize(0);
//	vtx = mesh->vertices[i];
        for (int j=0;j<(int)dims_out[1];j++){
            if(ttemp[i][j] >= 0)
                 mesh->vertices[i]->cells.push_back(ttemp[i][j]);
	}
    }
	

    printf("12\n");
    
    
}


int main(int argc, char* argv[]) {

// argv[1]: model.txt file
// argv[2]: model.mat file
// argv[3]: mesh.h5 file
// argv[4]: output.h5 file

    
#ifdef __MACH__
    uint64_t seed;
    seed = mach_absolute_time();
    srand48(seed);
#elif __linux__
    timespec t;
    clock_gettime(CLOCK_REALTIME,&t);
    unsigned long seed = (unsigned long)(1e9*t.tv_sec+t.tv_nsec);
    srand48(seed);
#endif

	
    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    simulation sim;
    /* Set some default values. */
    sim.ntraj = 1;
    sim.ncores = 1;
    
    vector <species> specs;
    vector <association> assocs;
    vector <dissociation> dissocs;
    vector <birth> births;
    vector <parameter> parameters;
    
    /* Read model. See documentation for full specification. */
    parse_model(argv[1],&sim,specs,assocs,dissocs,births,parameters);
    dimension = sim.dimension;

    // Read in the legacy urdme_model datastructure. This is still needed for some of the routines.
    char *urdmeinputfile;
    urdmeinputfile= argv[2]; 
    urdme_model *model;
    model = read_model(urdmeinputfile);
    model->infile = urdmeinputfile;
    
    if (model == NULL){
        printf("Fatal error. Couldn't load model file or currupt model file.");
        return(-1);
    }

    fem_mesh *mesh;
    mesh = (fem_mesh *)malloc(sizeof(fem_mesh));

    string mesh_file_name = argv[3];
    H5File mesh_file = H5File(mesh_file_name, H5F_ACC_RDONLY );
    read_p(mesh, mesh_file);
    read_t(mesh, mesh_file);
    read_bnd(mesh, mesh_file);
    printf("read_vertex_to_cell\n");
    read_vertex_to_cell(mesh, mesh_file);
    printf("read_vertex_to_cell, done\n");

    /* Initialize the primal/dual mesh format. Do we need this for pure micro?? */
    mesh_primal2dual(mesh);

    /* Compute all the planes that approximates the boundaries */

    /*for (int i=0;i<mesh->Ncells;i++){
        print_vertex(mesh->vertices[i]);
    } */  

    vector <plane> boundaries;
    boundaries = voxel_boundaries(model, mesh);

    // File where the output will be stored
    string output_filename = argv[4];

    /* Create output directory. */
  
   /* Do simulations. */
    for(int l=0;l<sim.ntraj;++l){

        /* Initialize random number generator. */
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
        long int seed_loc = lrand48();
        gsl_rng_set(rng,seed_loc);
        
        /* Initialize molecules. */
        group grp;
        for(int i=0;i<(int)(specs.size());i++){
            generate_particles(&grp,mesh,specs[i].initial_value,i,rng);
        }

        /* initialize the positions. TODO: Move inside generate particles/create particles */
        particle *part;
        for (int i=0;i<(int)(grp.particles.size());i++)
        {
            part = &(grp.particles[i]);
            part->dim=3;
            micro2meso(part, mesh);
         }

        
        double T = sim.T;
        double dt = sim.T/sim.num_intervals;
        double t = 0.0;
    
        int num_specs = (int)(specs.size());
    
    
        H5File file = H5File(output_filename, H5F_ACC_TRUNC );

        /*
        * Create the base group in the file
        */
        Group trajectory_group = file.createGroup( "/Trajectories" );
        string prefix = "/Trajectories/"+to_string(l);
        Group trajectory_0 = file.createGroup( prefix );
        for (int spec=0; spec<num_specs;spec++){
            Group trajectory_0 = file.createGroup( prefix +"/Type_"+ to_string(spec));
        }


        int time_index = 0;
        add_time_point_to_file(file,grp,num_specs,l,time_index);

        while(t<T){

            /* TODO: We should check for birth processes here. */
            main_simulator(&grp,specs,assocs,dissocs,boundaries,dt,rng,l);
            reflect_boundary(grp.particles,boundaries);
//                    int Psize = (int)(grp.particles.size());
//                    printf("[");
//                    for(int q=0;q<Psize;q++){
//                        printf("%g %g %g;\n",grp.particles[q].pos[0],grp.particles[q].pos[1],grp.particles[q].pos[2]);
//                    }
//                    printf("];");
            
            t += dt;
//            printf("t=%g\n",t);
            time_index++;
            cout << "Timestep complete, t=" << t << "\n";
            // Write solution to file.
            add_time_point_to_file(file,grp,num_specs,l,time_index);

        }

        gsl_rng_free(rng);
    }
    

    //print_model(specs,assocs,dissocs,births,&sim);

    
    gettimeofday(&end, NULL);
    printf("Run time: %.5g\n",(end.tv_sec+1e-6*end.tv_usec)-(start.tv_sec+1e-6*start.tv_usec));
    return 0;
}

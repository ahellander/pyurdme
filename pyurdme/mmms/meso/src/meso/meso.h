#ifndef MESO_MAIN_H
#define MESO_MAIN_H

#include "binaryheap.h"
#include "../include/rates.h"
#include "../include/mesh.h"
#include "../../../include/structs.h"

#define pi 3.141592653589793
int RATES = 1;


//#define DIM 3

using namespace std;




bool compare_events2 (tent_event *e1,tent_event *e2) {
    return e1->t<e2->t;
}


int reversible(dissociation *dissoc,vector <association>& associations){
    
    if((int)(dissoc->products.size())==2){
        for(int i=0;i<(int)(associations.size());i++){
            if(associations[i].products.size()==1){
                if((associations[i].reactant1==dissoc->products[0] && associations[i].reactant2==dissoc->products[1]) || (associations[i].reactant2==dissoc->products[0] && associations[i].reactant1==dissoc->products[1])){
                    if(associations[i].products[0]==dissoc->reactant){
                        return i;
                    }
                }
            }
        }
    }
    return -1;
}


void update_diss(vector <voxel>& voxels,vector <association>& assocs,vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <dissociation>& dissociations,gsl_rng *rng,double t,int i,int MICRO){
    
    int M = (int)(dissociations.size());
    tent_event temp_event;
    temp_event.type = 0;
    double k = 0.0,ktemp=0.0;
    for(int j = 0;j<M;j++){
        if(dissociations[j].reactant==particles[i].type && voxels[particles[i].voxel].isbnd==1){
            temp_event.index = j;
            if(MICRO==1){
                k = dissociations[j].k;
            }
            else{
                k = dissociations[j].kmeso;
                if(k<0){
                    /* TODO: What happens if we only have one product, e.g. A<-->B?. */
                    if((int)(dissociations[j].products.size())>1){
                    
                        if(RATES==1){
                            double htemp = pow(voxels[particles[i].voxel].vol,1.0/3.0);
                            double sigmatemp = spec[dissociations[j].products[0]].sigma+spec[dissociations[j].products[1]].sigma;
                            double Dtemp = spec[dissociations[j].products[0]].D+spec[dissociations[j].products[1]].D;
                            k = k_d_meso(assocs[dissociations[j].rev].k,dissociations[j].k,htemp,sigmatemp,Dtemp,spec[particles[i].type].dim);
                        }
                        else{
                            ktemp = 1.0/(dissociations[j].kmeso_u[0]+voxels[particles[i].voxel].vol*dissociations[j].kmeso_u[1]);
                            k = voxels[particles[i].voxel].vol*dissociations[j].k*ktemp/assocs[dissociations[j].rev].k;
                        }
                    }
                    else{
                        k = dissociations[j].k;
                    }
                }
            }
            temp_event.t = -log(gsl_rng_uniform(rng))/k+t;
            temp_event.reactants.push_back(i);
            if(MICRO==1){
                heap_insert(tent_events,temp_event,compare_events2);
            }
            else{
                heap_insert(tent_events,temp_event,compare_events2);
            }
        }
    }
}


void clear_diff(vector <tent_event>& tent_events,int index){
    int N = (int)(tent_events.size());
    for(int i = 0;i<N;i++){
        if(tent_events[i].type==-1){
            if(tent_events[i].reactants[0]==index){
                heap_delete(tent_events,i,compare_events2);
                i--;
                N--;
            }
        }
    }
}



void clear_diss(vector <tent_event>& tent_events,int index){
    int N = (int)(tent_events.size());
    for(int i = 0;i<N;i++){

        if(tent_events[i].type==0){
            if(tent_events[i].reactants[0]==index){
                heap_delete(tent_events,i,compare_events2);
                i--;
                N--;
            }
        }

    }
}

void clear_assoc(vector <tent_event>& tent_events,int index){
    int N = (int)(tent_events.size());
    for(int i = 0;i<N;i++){
        if(tent_events[i].type==1){
            if(tent_events[i].reactants[0]==index || tent_events[i].reactants[1]==index){
                heap_delete(tent_events,i,compare_events2);
                i--;
                N--;
            }
        }
    }
}

void adjust_index(vector <tent_event>& tent_events,int index){
    int N = (int)(tent_events.size());
    for(int i = 0;i<N;i++){
        if(tent_events[i].reactants[0]>index){
            tent_events[i].reactants[0]--;
        }
        if(tent_events[i].type==1){
            if(tent_events[i].reactants[1]>index){
                tent_events[i].reactants[1]--;
            }
        }
    }
}

void add_diff_event(vector <particle>& particles,vector <species>& specs,vector <tent_event>& tent_events,vector <voxel>& voxels,gsl_rng *rng,double t,int i){
    double totalD = voxels[particles[i].voxel].totalD*specs[particles[i].type].mesoD;
//    printf("D_total=%g\n",totalD);
    if(totalD>0){
        tent_event temp_event;
        temp_event.type = -1;
        temp_event.reactants.resize(1);
        
            temp_event.t = -log(gsl_rng_uniform(rng))/totalD+t;
        
        temp_event.reactants[0] = i;
        heap_insert(tent_events,temp_event,compare_events2);
    }
}


void get_diff_events(vector <particle>& particles,vector <species>& specs,vector <tent_event>& tent_events,vector <voxel>& voxels,gsl_rng *rng,double t){
    int N = (int)(particles.size());
    for(int i = 0;i<N;i++){
        if(particles[i].meso_micro==0){
            add_diff_event(particles,specs,tent_events,voxels,rng,t,i);
        }
    }
}


void get_diss_events(vector <voxel>& voxels,vector <association>& assocs,vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <dissociation>& dissociations,gsl_rng *rng,double t){

    int N = (int)(particles.size());
    tent_event temp_event;
    temp_event.type = 0;
    temp_event.reactants.resize(1);
    for(int i = 0;i<N;i++){
        if(particles[i].meso_micro==0){
            update_diss(voxels,assocs,particles,spec,tent_events,dissociations,rng,t,i,0);
        }
    }
}



void check_association(vector <voxel>& voxels,vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <association>& associations,gsl_rng *rng,double t,int i,int j,int dist){

    // if(in_same_voxel(i,j,particles)){
    int M = (int)(associations.size());
    tent_event temp_event;
    temp_event.type = 1;
    temp_event.reactants.resize(2);
    double rrate = 0.0;
    for(int k=0;k<M;k++){
        if(((associations[k].reactant1==particles[i].type && associations[k].reactant2==particles[j].type) || (associations[k].reactant2==particles[i].type && associations[k].reactant1==particles[j].type))){

            temp_event.index = k;
            temp_event.reactants[0] = i;
            temp_event.reactants[1] = j;
            
            if(dist==0){
//                printf("rate is (model)=%g\n",associations[k].kmeso);
//                printf("r is (model)=%g\n",associations[k].r);
//                printf("dist is:%d\n",dist);
                double htemp=0,sigmatemp=0,Dtemp=0;
                if(RATES==1){
                    htemp = pow(voxels[particles[i].voxel].vol,1.0/3.0);
                    sigmatemp = spec[associations[k].reactant1].sigma+spec[associations[k].reactant2].sigma;
                    Dtemp = spec[associations[k].reactant1].D+spec[associations[k].reactant2].D;
                    
                    rrate = k_r_meso(associations[k].k,htemp,sigmatemp,Dtemp,3);
//                    printf("rrate=%g, h=%g, D=%g, sigma=%g, k=%g\n",rrate,htemp,Dtemp,sigmatemp,associations[k].k);
                }
                else{
                    rrate = 1.0/(associations[k].kmeso_u[0]+associations[k].kmeso_u[1]*voxels[particles[i].voxel].vol);//associations[k].kmeso*(1-2*DIM*associations[k].r);
                }
                double q = 1.1;
                
                /* Take smallest volume that will yield a positive reaction rate. */
                while(rrate<0){
                    if(RATES!=1){
                        rrate = 1.0/(associations[k].kmeso_u[0]+associations[k].kmeso_u[1]*voxels[particles[i].voxel].vol*q);
                        q += 0.1;
                    }
                    else{
                        rrate = k_r_meso(associations[k].k,htemp*q,sigmatemp,Dtemp,spec[particles[i].type].dim);
                        q += 0.1;
                    }
                }
                
            }
            else{
                rrate = associations[k].kmeso*associations[k].r;
            }
            temp_event.t = -log(gsl_rng_uniform(rng))/rrate+t;
//            if(dist==0){
//                printf("time to meso reac: %g, current_time=%g\n\n",temp_event.t,t);
//                printf("rate is %g\n",rrate);
//            }
            
            heap_insert(tent_events,temp_event,compare_events2);
            
//            tent_events.insert(lower_bound(tent_events.begin(),tent_events.end(),temp_event,compare_events_meso),temp_event);

        }
    }
    // }
}

void get_assoc_events(vector <voxel>& voxels,vector <particle>& particles,vector <species>& specs,vector <tent_event>& tent_events,vector <association> associations,gsl_rng *rng,double t){

    int N = (int)(particles.size());
    
    int voxel;
    
    
    int dist_temp;
    
    
    for(int i=0;i<N;i++){
        
        voxel = particles[i].voxel;
        for(int j=i+1;j<N;j++){
            dist_temp = 2;
            if(particles[j].voxel==voxel){
                dist_temp = 0;
            }
            if((dist_temp==0 || dist_temp==1) && (particles[i].meso_micro==0 || particles[j].meso_micro==0)){
                check_association(voxels,particles,specs,tent_events,associations,rng,t,i,j,dist_temp);
            }
        }
    }
//    free(voxel);
}

void update_assoc(vector <voxel>& voxels,vector <particle>& particles,vector <species>& specs,vector <tent_event>& tent_events,vector <association>& associations,gsl_rng *rng,double t,int index){
    
    int N = (int)(particles.size());
    
    int voxel;
    
    int dist_temp = 0;
    
    voxel = particles[index].voxel;
    for(int i = 0;i<N;i++){
        if(i!=index){
            //Get associations
            dist_temp = 2;
            if(particles[i].voxel==voxel){
                dist_temp = 0;
            }
            if((dist_temp==0 || dist_temp==1) && (particles[i].meso_micro==0 || particles[index].meso_micro==0)){
//                printf("Checking for association...\n");
                check_association(voxels,particles,specs,tent_events,associations,rng,t,i,index,dist_temp);
            }
        }
    }
//    free(voxel);

}

//void reflect(particle *a,int *boundary,int DIM){
//    for(int i = 0;i<DIM;i++){
//        if(a->voxel[i]<boundary[i*2]){
//            a->voxel[i] += 1;
//        }
//        else if(a->voxel[i]>boundary[i*2+1]){
//            a->voxel[i] -= 1;
//        }
//    }
//}

//void diffuse(particle *a,double rnum,double cutoff,int *boundary,int DIM){
//    for(int i = 1;i<=DIM;i++){
//        if(rnum<=cutoff*(2*i-1)){
//            a->voxel[i-1]++;
//            break;
//        }
//        else if(rnum<=cutoff*(2*i)){
//            a->voxel[i-1]--;
//            break;
//        }
//    }
//    reflect(a,boundary,DIM);
//}

int diffuse(int pos,vector <voxel>& voxels,gsl_rng *rng){
    
    double totalD = voxels[pos].totalD;
  
//    *dt = -log(gsl_rng_uniform(rng))/totalD;
    double temp = voxels[pos].neighbors[0].D;
    double rtemp = gsl_rng_uniform(rng);
    int i = 0;
    
    while(rtemp>temp/totalD){
        ++i;
        temp += voxels[pos].neighbors[i].D;
    }
    return voxels[pos].neighbors[i].vox;
}

void doDissociation(vector <particle>& particles,vector <voxel>& voxels,vector <species>& specs,vector <dissociation>& dissociations,vector <tent_event>& tent_events,int reaction,int reactant,double t,vector <association>& associations,gsl_rng *rng,int *UNIQUE_ID){

    int N  = (int)(dissociations[reaction].products.size());
//    printf("1\n");
    particle new_part;
//    new_part.voxel[0] = particles[reactant].voxel[0];
//    new_part.voxel[1] = particles[reactant].voxel[1];
//    new_part.voxel[2] = particles[reactant].voxel[2];
    
    new_part.voxel = particles[reactant].voxel;
    

    new_part.pos[0] = particles[reactant].pos[0];
    new_part.pos[1] = particles[reactant].pos[1];
    new_part.pos[2] = particles[reactant].pos[2];
//printf("2\n");

    particles.erase(particles.begin()+reactant);
    clear_assoc(tent_events,reactant);
    clear_diff(tent_events,reactant);
    clear_diss(tent_events,reactant);
//printf("3\n");
    adjust_index(tent_events,reactant);
//printf("4\n");
    /*TODO: Check if a molecule is at the boundary. Product can not end up outside the domain. */
    /* NOTE: N cannot be greater than 2 here. */
    for(int i = 0;i<N;i++){
     
        new_part.clock = -t;
        new_part.meso_micro = 0;
        new_part.unique_id = *UNIQUE_ID;
        *UNIQUE_ID = *UNIQUE_ID+1;
//        printf("41\n");
        new_part.type = dissociations[reaction].products[i];
        particles.push_back(new_part);
//        printf("42\n");
        update_assoc(voxels,particles,specs,tent_events,associations,rng,t,(int)(particles.size())-1);
//        printf("43\n");
        update_diss(voxels,associations,particles,specs,tent_events,dissociations,rng,t,(int)(particles.size())-1,0);
//        printf("44\n");
        add_diff_event(particles,specs,tent_events,voxels,rng,t,(int)(particles.size())-1);
    }
//printf("5\n");
}

void doAssociation(vector <particle>& particles,vector <voxel>& voxels,vector <species>& specs,vector <dissociation>& dissociations,vector <tent_event>& tent_events,int reaction,int reactant1,int reactant2,double t,vector <association>& associations,gsl_rng *rng,int *UNIQUE_ID){

    /* TODO: Determine the voxel that the molecules react in. */
//    printf("DOING ASSOCIATION\n\n");
//    printf("a1\n");
    int N  = (int)(associations[reaction].products.size());

    particle new_part;
    new_part.voxel = particles[reactant1].voxel;

    new_part.pos[0] = INFINITY;
    new_part.pos[1] = INFINITY;
    new_part.pos[2] = INFINITY;

    if(reactant2<reactant1){
        int temp = reactant2;
        reactant2 = reactant1;
        reactant1 = temp;
    }
//printf("a2\n");
    particles.erase(particles.begin()+reactant1);
    clear_assoc(tent_events,reactant1);
    clear_diff(tent_events,reactant1);
    clear_diss(tent_events,reactant1);
//printf("a3\n");
    particles.erase(particles.begin()+reactant2-1);
    clear_assoc(tent_events,reactant2);
    clear_diff(tent_events,reactant2);
    clear_diss(tent_events,reactant2);
//printf("a4\n");
    adjust_index(tent_events,reactant2);
    adjust_index(tent_events,reactant1);

//printf("a5\n");

    for(int i = 0;i<N;i++){
        new_part.clock = -t;
        new_part.meso_micro = 0;
        new_part.unique_id = *UNIQUE_ID;
//        printf("a5_2\n");
        *UNIQUE_ID = *UNIQUE_ID+1;
//        printf("a5_21\n");
//        printf("i=%d ",i);
//        printf("N=%d\n",N);
        new_part.type = associations[reaction].products[i];
//        printf("a5_22\n");
        particles.push_back(new_part);
//        printf("a5_3\n");
        update_assoc(voxels,particles,specs,tent_events,associations,rng,t,(int)(particles.size())-1);
        update_diss(voxels,associations,particles,specs,tent_events,dissociations,rng,t,(int)(particles.size())-1,0);
        add_diff_event(particles,specs,tent_events,voxels,rng,t,(int)(particles.size())-1);
//        printf("a5_4\n");
    }
//    printf("a6\n");
}


/*

double HHP_rrate_assoc(double D,double sigma,double k_r,double L,double h){

    double kme = 4*pi*sigma*D*k_r/(4*pi*sigma*D+k_r);
    double tau_micro = pow(L,3)/kme;
    double temp = 1.5164*pow(L,3)/(6*D*h);
    double k_meso = pow(L/h,3)/(tau_micro-temp);
    return k_meso;
}

double HHP_rrate_dissoc(double D,double sigma,double k_r,double k_d,double L,double h){
    return pow(h,3)*k_d*HHP_rrate_assoc(D,sigma,k_r,L,h)/k_r;
}*/


int meso_simulator(vector <particle>& particles,vector <species>& specs,vector <association>& associations,vector <dissociation>& dissociations,vector <voxel>& voxels,double T,gsl_rng *rng,int *UNIQUE_ID){

    
    
    
    double t = 0.0;
    int MESO=0;
    
    if((int)(particles.size())==0){
//        printf("No meso particles.\n");
    }
    else{
        
        for(int i=0;i<(int)(particles.size());i++){
            if(particles[i].meso_micro==0){
                MESO++;
            }
        }
    }

    if(MESO==0){
        return -1;
    }

    vector <tent_event> tent_events;
    
    
    get_diff_events(particles,specs,tent_events,voxels,rng,t);
    get_assoc_events(voxels,particles,specs,tent_events,associations,rng,t);
    get_diss_events(voxels,associations,particles,specs,tent_events,dissociations,rng,t);

    
    
    double t_next_event;
    int type_next_event;
    int reactant1,reactant2;
    int new_pos;
    
    tent_event next_event;
    
    while(t<T){
        //Get time for next event.
        
        next_event = tent_events[0];
        t_next_event = next_event.t;
        
//        t_next_event = tent_events.back().t;
        
        if(t_next_event>T){
            break;
        }
        type_next_event = next_event.type;
        reactant1 = next_event.reactants[0];
//        printf("Time until next event: %g, type=%d\n",t_next_event,type_next_event);
        
        
        
        
        
        //Next event diffusion.
        if(type_next_event==-1){
//            diffuse(&(particles[reactant1]),gsl_rng_uniform(rng),1.0/(2.0*DIM),DIM);
            new_pos = diffuse(particles[reactant1].voxel,voxels,rng);
            particles[reactant1].voxel = new_pos;
        
//            tent_events.pop_back();
            heap_delete(tent_events,0,compare_events2);
            clear_assoc(tent_events,reactant1);
            t = t_next_event;
            update_assoc(voxels,particles,specs,tent_events,associations,rng,t,reactant1);
            
            
            /* TODO: This is not always necessary. Used for reactions in specific subdomains. */
            clear_diss(tent_events,reactant1);
            update_diss(voxels,associations,particles,specs,tent_events,dissociations,rng,t,reactant1,0);
            
            add_diff_event(particles,specs,tent_events,voxels,rng,t,reactant1);
        }
        else if(type_next_event==0){
            t = t_next_event;
            printf("%g %g %g\n",voxels[particles[reactant1].voxel].node[0],voxels[particles[reactant1].voxel].node[1],voxels[particles[reactant1].voxel].node[2]);
            doDissociation(particles,voxels,specs,dissociations,tent_events,next_event.index,reactant1,t,associations,rng,UNIQUE_ID);
            //            T_dissociation = T_global+t;
        }
        else{
//            reactant2 = tent_events.back().reactants[1];
            reactant2 = next_event.reactants[1];
            t = t_next_event;
//            printf("particles before: %d\n",(int)(particles.size()));
            
//            doAssociation(particles,voxels,specs,dissociations,tent_events,tent_events.back().index,reactant1,reactant2,t,associations,rng,UNIQUE_ID);
            
            doAssociation(particles,voxels,specs,dissociations,tent_events,next_event.index,reactant1,reactant2,t,associations,rng,UNIQUE_ID);
//            printf("particles after: %d\n",(int)(particles.size()));

        }
    }
    return 1;
}


#endif

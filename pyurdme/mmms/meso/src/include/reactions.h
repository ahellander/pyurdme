#ifndef REACTIONS_H
#define REACTIONS_H

#include "../meso/binaryheap.h"

int RATES;

void update_diss(vector <voxel>& voxels,vector <association>& assocs,vector <particle>& particles,vector <species>& spec,vector <tent_event>& tent_events,vector <dissociation>& dissociations,gsl_rng *rng,double t,int i,int MICRO){

    int M = (int)(dissociations.size());
    tent_event temp_event;
    temp_event.type = 0;
    double k = 0.0,ktemp=0.0;
    for(int j = 0;j<M;j++){
        if(dissociations[j].reactant==particles[i].type){
            temp_event.index = j;
            if(MICRO==1){
                k = dissociations[j].k;
            }
            else{
                k = dissociations[j].kmeso;
                if(k<0){
                    /* TODO: What happens if we only have one product, e.g. A<-->B?. */
                    if(RATES==1){
                        double htemp = pow(voxels[particles[i].vox_id].vol,1.0/3.0);
                        double sigmatemp = spec[dissociations[j].products[0]].sigma+spec[dissociations[j].products[1]].sigma;
                        double Dtemp = spec[dissociations[j].products[0]].D+spec[dissociations[j].products[1]].D;
                        k = k_d_meso(assocs[dissociations[j].rev].k,dissociations[j].k,htemp,sigmatemp,Dtemp,spec[particles[i].type].dim);
                    }
                    else{
                        ktemp = 1.0/(dissociations[j].kmeso_u[0]+voxels[particles[i].vox_id].vol*dissociations[j].kmeso_u[1]);
                        k = voxels[particles[i].vox_id].vol*dissociations[j].k*ktemp/assocs[dissociations[j].rev].k;
                    }
                }
            }
            temp_event.t = -log(gsl_rng_uniform(rng))/k+t;
            temp_event.reactants.push_back(i);
            if(MICRO==1){
                heap_insert(tent_events,temp_event,compare_events2);
//                tent_events.insert(lower_bound(tent_events.begin(),tent_events.end(),temp_event,compare_events),temp_event);
            }
            else{
                heap_insert(tent_events,temp_event,compare_events2);
//                tent_events.insert(lower_bound(tent_events.begin(),tent_events.end(),temp_event,compare_events_meso),temp_event);
            }
        }
    }
}



#endif

#ifndef STRUCTS_H
#define STRUCTS_H

using namespace std;

int grp_id,p_id=0;

typedef struct neighbor{
    int vox;
    double D;
}neighbor;

typedef struct voxel {
    int id;
    double node[3];
    int num_conn;
    vector <neighbor> neighbors;
    double vol;
    int sd;
    double totalD;
}voxel;

typedef struct simulation{
    int ntraj;
    int ncores;

    int dimension;

    double T;
    int num_intervals;
    double boundary[6];
    double volume;

    vector <voxel> voxels;

//    double h_x,h_y,h_z;
    int num_voxels;
    char name[MAX_CHARACTERS];
    char mesh[MAX_CHARACTERS];

}simulation;

typedef struct species{
    double D,mesoD;
    double sigma;
    int dim;
    char name[MAX_CHARACTERS];
    int initial_value;

    /* 0 is meso, 1 is micro. */
    int meso_micro;
    /* Minimum time a newly created mesoscopic molecule stays microscopic before it becomes mesoscopic. */
    double min_micro;
}species;

typedef struct particle{
    int type;
    int voxel[3];
    double pos[3];
    double vec1[3],vec2[3],vec3[3];
    bool active;
    int unique_id;
    
    int vox_id;
    
    /* 0 is meso, 1 is micro. */
    int meso_micro;

    double clock;
}particle;

typedef struct group{
    vector <particle> particles;
    double dist;
    int unique_id;
}group;

typedef struct dissociation{
    int reactant;
    vector <int> products;
    double k,kmeso;
    double kmeso_u[2];
    int rev;
    double r;
}dissociation;

typedef struct association{
    int reactant1,reactant2;
    vector <int> products;
    double k,kmeso;
    double kmeso_u[2];
    double r;
}association;

typedef struct birth{
    int product;
    double k;
}birth;

typedef struct tent_event{
    int type; //-1 diffusion, 0 dissociation, 1 association
    int index;
    double t;
    vector <int> reactants;

}tent_event;

typedef struct parameter{
    double value;
    char name[MAX_CHARACTERS];
}parameter;

typedef struct node{
    double p[3];

}node;

typedef struct tetrahedron{
    int type;
    int nodes[4];
}tetrahedron;

typedef struct plane{
    double p[3];
    double v1[3];
    double v2[3];
    double n[3];
    int cn;
}plane;

bool compare_events (tent_event e1,tent_event e2) {
    return e1.t<e2.t;
}

bool compare_events2 (tent_event *e1,tent_event *e2) {
    return e1->t<e2->t;
}

bool compare_events_meso (tent_event e1,tent_event e2) {
    return e1.t>e2.t;
}

#endif

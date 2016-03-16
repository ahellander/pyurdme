#ifndef STRUCTS_H
#define STRUCTS_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sstream>

/* Gnu scientific library */
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_errno.h"

#include "global_params.h"
#include "mesh.h"

using namespace std;

//int grp_id,p_id=0;

typedef struct simulation{
    int ntraj;
    int ncores;
    
    int dimension;
    
    double T;
    int num_intervals;
    double boundary[6];
    double volume;
    
    int num_voxels;
    char name[MAX_CHARACTERS];
    
}simulation;

typedef struct species{
    double D;
    double sigma;
    int dim;
    char name[MAX_CHARACTERS];
    int initial_value;
}species;

typedef struct particle{
    int type;
    int voxel[3];
    double pos[3];
    double vec1[3],vec2[3],vec3[3];
    bool active;
    int unique_id;
}particle;

typedef struct group{
    vector <particle> particles;
    double dist;
    int unique_id;
}group;

typedef struct dissociation{
    int reactant;
    vector <int> products;
    double k;
    double r;
}dissociation;

typedef struct association{
    int reactant1,reactant2;
    vector <int> products;
    double k;
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

typedef struct plane{
    int id;
    double p[3];
    double v1[3];
    double v2[3];
    double n[3];
    int cn;
    int isbnd; 
    vector <int> type;

}plane;


#endif

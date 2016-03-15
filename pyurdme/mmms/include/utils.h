#ifndef UTILS_H
#define UTILS_H

#include "structs.h"
#include <vector>
#include <math.h>
using namespace std;


double dist3(double *v1,double *v2);
void print_pos(double *pos);
void print_id(int ID);
void print_group(group *grp);
void print_neighbors(int *n,int N);
void print_group_sizes(vector <group>& grps);
void print_molecules(group *grp,double dt,FILE *f_out);
void print_molecule(group *grp,int ID,double dt,FILE *f_out);
void print_tent_events(vector <tent_event>& t_events);
void print_vec(double *v);
void print_vecs(particle *p);
void cross(double *v1,double *v2,double *out);
double get_boundary_dist_pair(double dt,vector <particle>& particles,vector <species>& specs,double *boundary, int dimension);
void check_positions(group *grp);
void generate_particles(group *grp,double *boundary,int N,int type,gsl_rng *rng);
double dot(double *p1,double *p2);
void vec_diff(double *p1,double *p2,double *outv);
void rotation(double *v,double *z,double theta,double r);
void reflect_cuboid(group *grp,double *boundary,int dimension);
void reflect_cuboid_p(particle *p,double *boundary);
void copy_pos(double *v,double *src);
void copy_vecs(double *vec1,double *vec2,double *vec3,particle *p);
void diffuse_vec(double *pos,double *vec1,double *vec2, double *vec3,double D,double sigma,double dt,gsl_rng *rng);
void diffuse(particle *p,vector <species>& specs,double dt,gsl_rng *rng);
void divide_into_cubes(group *grp,vector <group>& grps,int Nx,int Ny,int Nz,double *boundary);
void sort_molecules(group *grp,vector <group>& grps,vector <species>& specs);
int cp(const char *to, const char *from);


#endif

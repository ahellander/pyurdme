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
#include "gsl/gsl_fit.h"

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#include "include/global_params.h"
#include "include/structs.h"
#include "include/model_parser.h"
#include "include/utils.h"
#include "include/reactions.h"
#include "include/mesh.h"
#include "include/est_rates.h"

#include "meso/meso.h"

#include "omp.h"

#define pi 3.141592653589793




int main(int argc, char* argv[]){

//    if(argc==1){
//        printf("No input model set.\n");
//        exit(EXIT_FAILURE);
//    }
//    RATES = 0;
//    if(argc==3){
//        RATES = atoi(argv[2]);
//    }
//    printf("Num args=%d\n",argc);
//    printf("RATES=%d\n",RATES);

    
    RATES = 1;
    
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

    simulation sim;
    /* Set some default values. */
    sim.ntraj = 1;
    sim.ncores = 1;
    /* Some initialization. */
    sim.dimension = -1;
    sim.num_voxels = -1;
    sim.volume = -1;

    vector <species> specs;
    vector <association> assocs;
    vector <dissociation> dissocs;
    vector <birth> births;
    vector<parameter> parameters;

    parse_model(argv[1],&sim,specs,assocs,dissocs,births,parameters);
    for(int i=0;i<(int)(specs.size());++i){
        specs[i].mesoD = specs[i].D;
    }

//    FILE *fmesh;
//    fmesh = fopen(sim.mesh,"r");
//    if(fmesh==NULL){
//        printf("Failed to open mesh file %s: %s\n",sim.mesh,strerror(errno));
//        exit(0);
//    }
//    printf("Reading mesh: %s...\n",sim.mesh);
//    read_mesh(fmesh,sim.voxels,&sim);
    
    sim.volume = 6.7181e-20;
    cube_cartesian(sim.voxels,&sim,1.0,8);
    
    
    printf("Number of voxels: %d\n",(int)(sim.voxels.size()));
    printf("Done reading mesh.\n");

    printf("Computing rates...\n");


    gettimeofday(&start, NULL);

    gsl_rng *rng_rates = gsl_rng_alloc(gsl_rng_taus);
    long int seed_loc_rates = lrand48();
    gsl_rng_set(rng_rates,seed_loc_rates);




    int num_samples = 500;
    int Nest = 500;
    double std1 = 8.0;
    double std2 = 8.0;

    printf("num_samples=%d\n",num_samples);
    printf("Nest=%d\n",Nest);
    printf("std1=%g\n",std1);
    printf("std2=%g\n",std2);

    int *vox_id_out = (int *)malloc(num_samples*sizeof(int));
    double *vol_out = (double *)malloc(num_samples*sizeof(double));
    double *est_rebind_out = (double *)malloc(num_samples*sizeof(double));
    double *est_inf_out = (double *)malloc(num_samples*sizeof(double));
    double *rates_est_out = (double *)malloc(num_samples*sizeof(double));
    
    
    
    
    double Dtemp,katemp,sigma_temp;
    double kmeso;
    double tau_micro;
    double c0,c1,cov00,cov01,cov11,sumsq;

    char filenamerates[MAX_CHARACTERS];
    memset(filenamerates,'\0',sizeof(filenamerates));
    strcpy(filenamerates,sim.mesh);
    strcat(filenamerates,".rates");
    
    if(RATES==0){

        FILE *frates;
        printf("Reading rates from file...\n");
        if((frates = fopen(filenamerates,"r"))){
            //Read reaction rates from file.
            fscanf(frates,"%d\n",&num_samples);
    //        double temp1,temp2,vol_temp;
            double foo;
            for(int i=0;i<num_samples;i++){
                fscanf(frates,"%lf %lf %lf %lf %d\n",&(vol_out[i]),&foo,&(est_rebind_out[i]),&(est_inf_out[i]),&(vox_id_out[i]));
            }
            fclose(frates);
            printf("Finished reading rates.\n");
        }
        else{
            //Compute reaction rates.
            printf("Computing rates: %s\n",filenamerates);
            frates = fopen(filenamerates,"w");
            estimate_rates(sim.voxels,&sim,Nest,rng_rates,vol_out,est_rebind_out,est_inf_out,vox_id_out,num_samples,std1,std2);

            fprintf(frates,"%d\n",num_samples);
            for(int i=0;i<num_samples;i++){
                fprintf(frates,"%.5g %.5g %.5g %.5g %d\n",vol_out[i],sim.voxels[vox_id_out[i]].totalD,est_rebind_out[i],est_inf_out[i],vox_id_out[i]);
            }
            fclose(frates);
            printf("Finished writing rates to file.\n");
            
        }
        

        for(int i=0;i<(int)(assocs.size());++i){
            katemp = assocs[i].k;
            Dtemp = specs[assocs[i].reactant1].D+specs[assocs[i].reactant2].D;

            sigma_temp = specs[assocs[i].reactant1].sigma+specs[assocs[i].reactant2].sigma;
            kmeso = 4*pi*sigma_temp*Dtemp*katemp/(4*pi*sigma_temp*Dtemp+katemp);
            tau_micro = sim.volume/kmeso;
            
            for(int j=0;j<num_samples;j++){

                rates_est_out[j] = sim.voxels[vox_id_out[j]].totalD*est_rebind_out[j]/(tau_micro-est_inf_out[j]/Dtemp);
                rates_est_out[j] = 1.0/rates_est_out[j];
            }
            gsl_fit_linear(vol_out,1,rates_est_out,1,num_samples,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
            printf("Best fit: 1/VOXEL_RATE=%.5g+%.5g*VOXEL_VOL\n",c0,c1);
            assocs[i].kmeso_u[0] = c0;
            assocs[i].kmeso_u[1] = c1;
        }
    }
    else if(RATES==1){
        // This should be for HHP rates.
        // I think we can set assocs[i].kmeso to HHP rate, make RATES a global variable and then change the relevant functions.
    }
    else if(RATES==2){
        double Dtemp,katemp,sigma_temp,kmeso;
        for(int i=0;i<(int)(assocs.size());++i){
            katemp = assocs[i].k;
            Dtemp = specs[assocs[i].reactant1].D+specs[assocs[i].reactant2].D;
            sigma_temp = specs[assocs[i].reactant1].sigma+specs[assocs[i].reactant2].sigma;
            kmeso = 4*pi*sigma_temp*Dtemp*katemp/(4*pi*sigma_temp*Dtemp+katemp);
            assocs[i].kmeso_u[0] = 0.0;
            assocs[i].kmeso_u[1] = 1.0/kmeso;
        }

    }


    for(int i=0;i<(int)(dissocs.size());++i){
        int rev = reversible(&(dissocs[i]),assocs);
        if(rev>=0){
            dissocs[i].kmeso = -1;
            dissocs[i].kmeso_u[0] = assocs[rev].kmeso_u[0];
            dissocs[i].kmeso_u[1] = assocs[rev].kmeso_u[1];
            dissocs[i].rev = rev;
        }
        else{
            dissocs[i].kmeso = dissocs[i].k;
        }
    }

    printf("Done computing rates.\n");
    gettimeofday(&end, NULL);
    printf("\nComputed rates in: %.5gs\n",(end.tv_sec+1e-6*end.tv_usec)-(start.tv_sec+1e-6*start.tv_usec));

    print_model(specs,assocs,dissocs,births,&sim);
    omp_set_num_threads(sim.ncores);

    /* Create output directory. */

    struct tm *tm;
    char str_date[100];
    time_t curr_time;
    curr_time = time(NULL);
    tm = localtime(&curr_time);
    strftime(str_date,sizeof(str_date),"%Y-%m-%d-%H-%M-%S",tm);

    char filename1[MAX_CHARACTERS];
    memset(filename1,'\0',sizeof(filename1));
    strcat(filename1,"out/");
    mkdir(filename1,0777);
    strcat(filename1,sim.name);
    mkdir(filename1,0777);
    strcat(filename1,"/");
    strcat(filename1,str_date);
    char filename_temp[MAX_CHARACTERS];
    memset(filename_temp,'\0',sizeof(filename_temp));
    strcpy(filename_temp,filename1);
    int FF = 2;
    char num_temp[MAX_CHARACTERS];
    while(mkdir(filename_temp,0777)==-1){
        strcpy(filename_temp,filename1);
        sprintf(num_temp,"_%d",FF);
        strcat(filename_temp,num_temp);
    }



    /* Copy model file. We assume that the model file has not been changed since start of execution. */
    char model_file[MAX_CHARACTERS];
    char model_file_out[MAX_CHARACTERS];
    sprintf(model_file,"/%s.txt","model");
    strcpy(model_file_out,filename_temp);
    strcat(model_file_out,model_file);
    cp(model_file_out,argv[1]);


    /* Copy mesh file. This could be a bad idea if the mesh is fine-grained. */
    char mesh_file[MAX_CHARACTERS];
    char mesh_file_out[MAX_CHARACTERS];
    sprintf(mesh_file,"/%s.txt","mesh");
    strcpy(mesh_file_out,filename_temp);
    strcat(mesh_file_out,mesh_file);
    cp(mesh_file_out,sim.mesh);

    omp_lock_t writelock;
    omp_init_lock(&writelock);




    /* Do simulations. */

    //double TOTAL_TIME = sim.ntraj*sim.T;
    double TOTAL_T = 0.0;

    gettimeofday(&start, NULL);
    
    
    fflush(stdout);

#pragma omp parallel for
    for(int l=0;l<sim.ntraj;++l){

        int UNIQUE_ID = 0;


         /* Initialize random number generator. */
        omp_set_lock(&writelock);
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
        long int seed_loc = lrand48();
        gsl_rng_set(rng,seed_loc);
        omp_unset_lock(&writelock);


        /* Open appropriate output streams. */
        char traj_num[30];
        sprintf(traj_num,"/%d_total.txt",l);
        char outputfile[MAX_CHARACTERS];
        char outputpos[MAX_CHARACTERS];
        char outputvoxel[MAX_CHARACTERS];
        memset(outputfile,'\0',sizeof(outputfile));
        strcpy(outputfile,filename_temp);
        strcpy(outputpos,filename_temp);
        strcpy(outputvoxel,filename_temp);
        strcat(outputfile,traj_num);
        sprintf(traj_num,"/%d_pos.txt",l);
        strcat(outputpos,traj_num);
        sprintf(traj_num,"/%d_vox.txt",l);
        strcat(outputvoxel,traj_num);

        FILE *file_total;
        file_total = fopen(outputfile,"w");

        FILE *file_pos;
        file_pos = fopen(outputpos,"w");

        FILE *file_vox;
        file_vox = fopen(outputvoxel,"w");


        group grp;


        /* Initialize molecules. */
        for(int i=0;i<(int)(specs.size());i++){
            generate_particles_meso(&grp,sim.voxels,specs[i].initial_value,i,rng,&UNIQUE_ID);
        }



        /* Initialize particles to either mesoscopic or microscopic. */
        for(int i=0;i<(int)(grp.particles.size());i++){
            grp.particles[i].meso_micro = 0;
        }
        /* ************************************************** */




        double T = sim.T;
        double dt;
        dt = sim.T/sim.num_intervals;

        double t = 0.0;

        int num_specs = (int)(specs.size());
        int *num_mol;
        num_mol = (int *)malloc(sizeof(int)*num_specs);

        for(int j=0;j<num_specs;j++){
            num_mol[j] = 0;
        }

        for(int j=0;j<(int)(grp.particles.size());j++){
            num_mol[grp.particles[j].type]++;
        }
        for(int j=0;j<num_specs;j++){
            fprintf(file_total,"%d ",num_mol[j]);
        }
        for(int j=0;j<num_specs;j++){
            num_mol[j] = 0;
        }

        fprintf(file_total,"%.5g ",0.0);
        fprintf(file_total,"\n");


        print_molecules(&grp,0.0,file_pos);
        fflush(file_pos);
        print_molecules_voxel(&grp,0.0,file_vox);
        fflush(file_vox);

        group meso_temp;

        while(t<T){

            // if(1==1){
            //     gettimeofday(&end, NULL);
            //     double time_exec = (end.tv_sec+1e-6*end.tv_usec)-(start.tv_sec+1e-6*start.tv_usec);
            //     printf("\rEstimated time left: %g",time_exec*(TOTAL_TIME/TOTAL_T-1.0));
            //     fflush(stdout);
            // }



            meso_simulator(grp.particles,specs,assocs,dissocs,sim.voxels,dt,rng,&UNIQUE_ID);


            t += dt;
            TOTAL_T += dt;

//            print_molecules(&grp,t,file_pos);
            print_molecules_voxel(&grp,0.0,file_vox);

            for(int j=0;j<(int)(grp.particles.size());j++){
                num_mol[grp.particles[j].type]++;
            }
            for(int j=0;j<num_specs;j++){
                fprintf(file_total,"%d ",num_mol[j]);
            }
            fprintf(file_total,"%.5g ",t);
            fprintf(file_total,"\n");
            for(int j=0;j<num_specs;j++){
                num_mol[j] = 0;
            }
            fflush(file_total);
            fflush(file_pos);
            fflush(file_vox);
        }
        gsl_rng_free(rng);
        free(num_mol);
        fclose(file_total);
        fclose(file_pos);
    }


    gettimeofday(&end, NULL);
    printf("\nRun time: %.5g\n",(end.tv_sec+1e-6*end.tv_usec)-(start.tv_sec+1e-6*start.tv_usec));

    return 0;
}

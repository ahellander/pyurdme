#ifndef EST_RATES_H
#define EST_RATES_H

int diffuse_rrate(int pos,vector <voxel>& voxels,double *dt,gsl_rng *rng){
    double totalD = voxels[pos].totalD;

    *dt = 1.0/totalD;
    double temp = voxels[pos].neighbors[0].D;
    double rtemp = gsl_rng_uniform(rng);
    while(rtemp==1.0){
        rtemp = gsl_rng_uniform(rng);
    }
    int i = 0;
    
    while(rtemp>temp/totalD){
        ++i;
        temp += voxels[pos].neighbors[i].D;
    }
    return voxels[pos].neighbors[i].vox;
}

void estimate_rates(vector <voxel>& voxels,simulation *sys,int N,gsl_rng *rng,double *vol_out,double *est_rebind_out,double *est_inf_out,int *vox_id_out,int num_samples,double std1,double std2){

    int counter = 0;
    //int total = N*num_samples;

#pragma omp parallel for
    for(int ll=0;ll<num_samples;++ll){//(int)(voxels.size());++ll){

        int *res = (int *)malloc(N*sizeof(int));
        double *times = (double *)malloc(N*sizeof(double));

        int vox_init;
        int pos;

        double D = 1.0;

        // printf("\r%.3g",100*(double)counter/total);
        // fflush(stdout);
        vox_init = floor((int)(voxels.size()*gsl_rng_uniform(rng)));//ll;

        /* ********************************************** */
        /* ********************************************** */
        /* *********** Estimate tau_infinity. *********** */
        /* ********************************************** */
        /* ********************************************** */
        double max_t1 = std1*pow(sys->volume,2.0/3.0)/(6*D);
        double Pr = 0.0;
        int num_reac = 0;
        double times_total = 0.0;
        for(int i=0;i<N;i++){
            counter = counter+1;
            int num_steps = 0;
            double t = 0.0;
            pos = floor((int)(voxels.size())*gsl_rng_uniform(rng));
//            printf("%d\n",pos);
            double dt = 0.0;
            while(pos!=vox_init && t<max_t1){
                pos = diffuse_rrate(pos,voxels,&dt,rng);
//                printf("%d\n",pos);

                t += dt;
                ++num_steps;
            }//while(pos!=vox_init && t<max_t1);
            if(pos==vox_init){//t<max_t1){
                times_total += t;
                ++num_reac;
            }
            res[i] = num_steps;
            //            times[i] = t;
//            printf("%d\n",num_steps);
        }
        Pr = (double)num_reac/(double)N;
//        printf("Pr=%g\n",Pr);
//        int total = 0;
//        double time_total = 0.0;
//        for(int i=0;i<N;++i){
//            total += res[i];
//        }
        double est_inf = 0.0;

        est_inf = (double)times_total/num_reac+(1.0-Pr)/Pr*max_t1;
        // printf("Pr=%g\n",Pr);
        fflush(stdout);
//        double cart_inf = sys->volume/6*1.5164*pow(voxels[vox_init].vol,-1.0/3.0);


        /* ********************************************** */
        /* ********************************************** */
        /* *********** Estimate tau_0^0. **************** */
        /* ********************************************** */
        /* ********************************************** */
        double max_t2 = std2*pow(sys->volume,2.0/3.0)/(6*D);
        for(int i=0;i<N;i++){
            int num_steps = 0;
            double t = 0.0;
            pos = vox_init;
            double dt = 0.0;
            do{
                pos = diffuse_rrate(pos,voxels,&dt,rng);
                t += dt;
                ++num_steps;
            }while(t<max_t2 && pos!=vox_init);
            if(t>=max_t2){
                t += est_inf;
            }
            res[i] = num_steps;
            times[i] = t;
        }

        // total = 0;
        double time_total = 0.0;
        for(int i=0;i<N;++i){
            // total += res[i];
            time_total += times[i];
        }
        double est_rebind = (double)time_total/N;
//        double cart_rebind = sys->volume/6*pow(voxels[vox_init].vol,-1.0/3.0);
        /* ********************************************** */
        /* ********************************************** */
        /* ********************************************** */
        /* ********************************************** */
        /* ********************************************** */
//        double kr = 1.0;
//        double sigma = 2e-3;
//        printf("****************************************\n");
//        printf("****************************************\n");
//        printf("Number of voxels: %d\n",(int)(voxels.size()));
//        printf("System volume=%g, voxel_volume=%g\n",sys->volume,voxels[vox_init].vol);
//        printf("Pr=%.5g\n",Pr);
//        printf("tau_inf =%.5g, tau_rebind =%.5g, q =%g\n",est_inf,est_rebind,est_inf/est_rebind);
//        printf("cart_inf=%.5g, cart_rebind=%.5g, q =%g\n",cart_inf,cart_rebind,cart_inf/cart_rebind);
//        double kmeso = 4*pi*sigma*D*kr/(4*pi*sigma*D+kr);
//        double tau_micro = sys->volume/kmeso;
//        double rho = voxels[vox_init].totalD*est_rebind/(tau_micro-est_inf);
//        printf("tau_micro=%g, kmeso = %.5g, rho=%.5g\n",tau_micro,kmeso/voxels[vox_init].vol,rho);
//        printf("voxel_number=%d, node=(%.5g,%.5g,%.5g)\n",vox_init,voxels[vox_init].node[0],voxels[vox_init].node[1],voxels[vox_init].node[2]);
//        printf("****************************************\n");
//        printf("****************************************\n");

        vox_id_out[ll] = vox_init;
        vol_out[ll] = voxels[vox_init].vol;
        est_rebind_out[ll] = est_rebind;
        est_inf_out[ll] = est_inf;
//        rates_est_out[ll] = 1.0/rho;
        free(res);
        free(times);

    }

//    double c0,c1,cov00,cov01,cov11,sumsq;

//    gsl_fit_linear(vol_out,1,rates_est_out,1,num_samples,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
//    printf("Best fit: 1/VOXEL_RATE=%.5g+%.5g*VOXEL_VOL\n",c0,c1);
}


#endif

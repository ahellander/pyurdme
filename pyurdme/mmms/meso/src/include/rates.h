#ifndef RATES_H
#define RATES_H

#define pi 3.141592653589793
#define C_alpha_3 1.5164
#define C_alpha_2 0.1951




double bisection(double (*f)(double,void**),void **args,double a,double b){
    if(f(a,args)*f(b,args)>0){
        //No zero in interval [a,b].
     //  printf("fa=%g, fb=%g\n",f(a,args),f(b,args));
        return 0.0;
    }

    int M = 1500;
    double err_tol = 1e-11;
    double p = (a+b)/2.0;
    double err = f(p,args);
    int i = 0;


    while(fabs(err)>err_tol && i<M){

        if(err*f(a,args)<0){
            b = p;
        }
        else{
            a = p;
        }
        p = (a+b)/2.0;

        err = f(p,args);
        i++;
    }
    return p;
}

double tau_micro_rebind(double kr,double L,int DIM){
    return pow(L,DIM)/kr;
}

double tau_micro(double kr,double D,double sigma,double L,int DIM){
    if(DIM==2){
        double lambda = pow(4*pi/3,1.0/3.0)*sigma/L;
        double alpha = kr/(2*pi*D);
        double F = log(1.0/lambda)/pow(1-pow(lambda,2),2)-(3-pow(lambda,2))/(4*(1-pow(lambda,2)));
        return (1+alpha*F)/kr*pow(L,2);
    }
    else if(DIM==3){
        double kCK = 4*pi*sigma*D*kr/(4*pi*sigma*D+kr);
        return pow(L,3)/kCK;
    }
    return -1;
}

double Nsteps(int N,int DIM){
    if(DIM==2){
        return 1.0/pi*N*log(N)-0.8049*N;
    }
    else if(DIM==3){
        return 0.5164*N;
    }
    return -1;
}

double meso_tau1(double kmeso,double r,double D,double L,int N,int DIM){
    double h = L/pow(N,1.0/(double)DIM);

    double tauJ = pow(h,2)/(2*DIM*D);

    double kmeso0 = (1.0-2.0*DIM*r)*kmeso;
    /* TODO: Shouldn't this be 2*DIM*r*kmeso? */
    double kmeso1 = r*kmeso;

    double p0 = 1.0/kmeso0;
    double p1 = 1.0/kmeso1;

    return N*(p0+tauJ)*p1/(2*DIM*(p0+tauJ)+p1);
}

double meso_tau0(double kmeso,double r,double D,double L,int N,double tau1,int DIM){

    double h = L/pow(N,1.0/(double)DIM);
    double tauJ = pow(h,2)/(2*DIM*D);
    double kmeso0 = (1.0-2*DIM*r)*kmeso;
    double p0 = 1.0/kmeso0;
    double te0 = 1.0/(kmeso0+1.0/tauJ);
    return te0+p0/(p0+tauJ)*tau1;
}

double tau_meso(double kmeso,double r,double D,double L,int N,int DIM){

    double Ns1 = Nsteps(N,DIM);
    double h = L/pow(N,1.0/(double)DIM);
    double tauJ = pow(h,2)/(2*DIM*D);


    double tau1 = meso_tau1(kmeso,r,D,L,N,DIM);

    return Ns1*tauJ+tau1;
}

double tau_meso_rebind(double kmeso,double r,double D,double L,int N,int DIM){
    double tau1 = meso_tau1(kmeso,r,D,L,N,DIM);
    double tau0 = meso_tau0(kmeso,r,D,L,N,tau1,DIM);
    return tau0;
    double kmeso0 = (1.0-2*DIM*r)*kmeso;
    double kmeso1 = r*kmeso;
    double P0 = kmeso0/(kmeso0+2*DIM*kmeso1);
    double P1 = 1-P0;
    return P0*tau0+P1*tau1;
}

double tau_diff(double kmeso,void **args){//double r,double kr,double D,double sigma,double L,int N){
    double *r = (double*)args[0];
    double *kr = (double*)args[1];
    double *D = (double*)args[2];
    double *sigma = (double*)args[3];
    double *L = (double*)args[4];
    int *N = (int*)args[5];
    int *DIM = (int*)args[6];
    //printf("diff=%g\n",tau_meso(kmeso,*r,*D,*sigma,*L,*N)-tau_micro(*kr,*D,*sigma,*L));
    return tau_meso(kmeso,*r,*D,*L,*N,*DIM)-tau_micro(*kr,*D,*sigma,*L,*DIM);
}

double tau_rebind_diff(double r,void **args){//double r,double kr,double D,double sigma,double L,int N){
    //double *kmeso = (double*)args[0];
    args[0] = (void*)&r;
    double kmeso = bisection(tau_diff,args,0.00001,1e28);
    //    printf("kmeso=%g\n",kmeso);
    double *kr = (double*)args[1];
    double *D = (double*)args[2];
    // double *sigma = (double*)args[3];
    double *L = (double*)args[4];
    int *N = (int*)args[5];
    int *DIM = (int*)args[6];
    //    printf("kr=%g, D=%g, sigma=%g, L=%g, N=%d, r=%g\n",*kr,*D,*sigma,*L,*N,r);
    //printf("diff=%g\n",tau_meso(kmeso,*r,*D,*sigma,*L,*N)-tau_micro(*kr,*D,*sigma,*L));
    return tau_meso_rebind(kmeso,r,*D,*L,*N,*DIM)-tau_micro_rebind(*kr,*L,*DIM);
}





double G(double h,double sigma,int dim){
    if(dim==3){
        return C_alpha_3/(6.0*h)-1.0/(4*pi*sigma);
    }
    else if(dim==2){
        return 0.25*(3/(2*pi)+C_alpha_2)-1/(2*pi)*log(1/sqrt(pi)*h/sigma);
    }
    return -1;
}

double k_r_meso(double k_r,double h,double sigma,double D,int dim){
    return k_r/pow(h,dim)*pow(1-k_r/D*G(h,sigma,dim),-1);
}

double k_d_meso(double k_r,double k_d,double h,double sigma,double D,int dim){
    return pow(h,dim)*k_d*k_r_meso(k_r,h,sigma,D,dim)/k_r;
}

double k_r_meso_CK(double k_r,double h,double sigma,double D,int dim){
    if(dim==3){
        return 4*pi*sigma*D*k_r/(4*pi*sigma*D+k_r)/pow(h,dim);
    }
    return -1;
}

double k_d_meso_CK(double k_r,double k_d,double h,double sigma,double D,int dim){
    if(dim==3){
        return pow(h,dim)*k_d*k_r_meso_CK(k_r,h,sigma,D,dim)/k_r;
    }
    return -1;
}

void non_local_k_r(double L,int N,double sigma,double D,double k_micro,double *krmeso,double *r,int DIM){

    double r_temp = 0.00000001;
    void *args[7];
    args[0] = (void*)&r_temp;
    args[1] = (void*)&k_micro;
    args[2] = (void*)&D;
    args[3] = (void*)&sigma;
    args[4] = (void*)&L;
    args[5] = (void*)&N;
    args[6] = (void*)&DIM;
    double upper_l = 1.0/4.0-r_temp;

    double h = L/pow(N,1.0/DIM);

    double hstar = 5.0979*sigma;

    if(DIM==3){
        hstar = 3.1759*sigma;
        upper_l = 1.0/6.0-r_temp;
    }

    if(h>=hstar){
        *krmeso = k_r_meso(k_micro,h,sigma,D,DIM);
        *r = 0.0;
    }
    else{
        double r_est = bisection(tau_rebind_diff,args,r_temp,upper_l);
        //printf("r_est=%g\n",r_est);
        args[0] = (void*)&r_est;
        double kmeso_nl = bisection(tau_diff,args,0.01,1e22);
        if(kmeso_nl==-1){
            printf("No zero in interval\n");
        }
        else{
            *r = r_est;
            *krmeso = kmeso_nl;
        }
    }
}

double non_local_k_d(double L,int N,double sigma,double D,double k_r,double k_d,double *r_temp,int DIM){

    double krmeso_temp;
    non_local_k_r(L,N,sigma,D,k_r,&krmeso_temp,r_temp,DIM);

    double h = L/pow(N,1.0/DIM);
    return pow(h,DIM)*k_d*krmeso_temp/k_r;
}

#endif

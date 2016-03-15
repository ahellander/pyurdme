#ifndef RAND_GEN_H
#define RAND_GEN_H

#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_linalg.h"
#include <gsl/gsl_integration.h>
#include "gsl/gsl_poly.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_sf_gamma.h"


#define rand_gen_constant_pi 3.141592653589793



void randomSphere(gsl_rng *rng,double r,double *rand,int dimension)
{
    if(dimension==3){
        double r1 = gsl_ran_gaussian_ziggurat(rng,1.0);
        double r2 = gsl_ran_gaussian_ziggurat(rng,1.0);
        double r3 = gsl_ran_gaussian_ziggurat(rng,1.0);
        double norm = 1/sqrt(r1*r1+r2*r2+r3*r3);
        rand[0] = r*norm*r1;
        rand[1] = r*norm*r2;
        rand[2] = r*norm*r3;

    }
    else{
        double r1 = 2*gsl_rng_uniform(rng)-1;
        double r2 = 2*gsl_rng_uniform(rng)-1;
        double p2 = r1*r1+r2*r2;
        while(p2>=1){
            r1 = 2*gsl_rng_uniform(rng)-1;
            r2 = 2*gsl_rng_uniform(rng)-1;
            p2 = r1*r1+r2*r2;
        }
        
        rand[0] = r*(r1*r1-r2*r2)/p2;
        rand[1] = r*2*r1*r2/p2;
        rand[2] = 0.0;
//        printf("randnorm=%g\n",sqrt(pow(rand[0],2)+pow(rand[1],2)+pow(rand[2],2)));
    }
}

double gen_rnd_theta(double D,double t,double r, const int N,gsl_rng *rng)
{
	
	double rand_r = gsl_rng_uniform(rng);
	double dtr = D*t/(r*r);
	int M = round(sqrt(8/dtr));
	if(M==0)
	{
		M = 1;
	}
	/*
	 if(M>10000)
	 {
	 printf("M = %d\n",M);
	 }
	 */
	double *theta = (double *)malloc(N*sizeof(double));
	double right = (15*sqrt(2*D*t))/(rand_gen_constant_pi*r);//sqrt(2*D*t/(r_0*r_0))/(rand_gen_constant_pi*r_0*r_0);
	if(right>rand_gen_constant_pi)
		right = rand_gen_constant_pi;
	double h = right/N;
	for(int i = 0;i<N;i++)
	{
		theta[i] = (i+1)*h;
	}
	//double Dt = D*t;
	
	
	//double sum = 0.0;
	//double sumprev = sum;
	
	double *E1 = (double *)malloc(M*sizeof(double));
	double *T = (double *)malloc(M*sizeof(double));
	
	double cons;
	double cdf = 0;
	double cdf_prev;
	for(int i = 0;i<M;i++)
	{
		E1[i] = exp(-i*(i+1)*dtr)*(2*i+1);
	}
	
	double res,res_prev;
	
	int j = 0;
	res = 0.0;
	cons = 0.25*h*sin(theta[j]);
	/*for(int i = 0;i<M;i++)
	 {
	 res += gsl_sf_legendre_Pl(i,cos(theta[j]))*E1[i];
	 }
	 res *= cons;*/
	
	for(j = 1;j<N;j++)
	{
		res_prev = res;
		res = 0.0;
		cons = 0.25*h*sin(theta[j]);
		
		gsl_sf_legendre_Pl_array(M-1,cos(theta[j]),T);
		for(int i = 0;i<M;i++)
		{
			res += T[i]*E1[i];
		}
		res *= cons;
		cdf_prev = cdf;
		cdf += res_prev+res;
		if(cdf>rand_r)
		{
			free(theta);
			free(E1);
			free(T);
			return h*(j-1)+(cdf-rand_r)/(cdf-cdf_prev)*h;
		}
	}
	free(theta);
	free(E1);
	free(T);
	return right;
	return h*(j-1)+(cdf-rand_r)/(cdf-cdf_prev)*h;
	return right-h+(rand_r-cdf)/(cdf-cdf_prev)*h;
}


double random_time_num(double k_r,double r_0,double sigma,
					   double D,double t,gsl_rng *rng,double *t_reac,int dimension,int time_steps,int space_steps,double ran)
{
	if(t!=t)
	{
		printf("t is nan\n");
	}
    
    
	double K = 15;
    
	
	double k = t/(double)time_steps;
	double endr = r_0+sqrt(2*D*t*K);
	space_steps = round(space_steps*(r_0-sigma)/(endr-sigma));
	if(space_steps<=1){
		space_steps = 1;
	}
	double h = (r_0-sigma)/space_steps;
	
	
    if(r_0<sigma*(1+1/1000.0))
    {
        r_0 = sigma+sigma/1000.0;
    
        h = (r_0-sigma)/space_steps;
    
    }
    if(r_0+h>=endr)
	{
		endr = r_0+10*h;
	}
	
	int J = (int)(round((endr-sigma)/h));
	//printf("J=%d, space_steps=%d, r_0=%g, endr=%g\n",J,space_steps,r_0,endr);
	if(J<=1)
	{
    
		endr = r_0+10*h;
		J = (int)(round((endr-sigma)/h));
    
	}
	/*if(J>1e3)
     {
     printf("J = %d\n",J);
     
     h = (r_0-sigma)/20;
     J = (int)(round((endr-sigma)/h));
     printf("New J = %d\n",J);
     
     }*/
    
	
	double *pdf = (double *)malloc(J*sizeof(double));
	double *S = (double *)malloc(time_steps*sizeof(double));
	
	
	int d = dimension;
	
	double const_d = 0.0;
	if(d==3)
		const_d = 4*rand_gen_constant_pi*sigma*sigma*D;
	else if (d==2)
		const_d = 2*rand_gen_constant_pi*sigma*D;
	else if (d==1)
		const_d = D;
	
//	double C1 = 1/(1-(h*k_r)/(const_d));
	
    
	double *r = (double *)malloc(J*sizeof(double));
	for(int i = 0;i<J;i++)
	{
		r[i] = sigma+(i+1)*h;
	}
	
	double t_init = pow(r_0-sigma,2)/(20*D);
    
	t_init = 0.0;
    
	int T = (int)round((t-t_init)/k);
	double *p = (double *)calloc(J,sizeof(double));
	double *p_prev = (double *)malloc(J*sizeof(double));
    
	int INIT = (int)(round((r_0-sigma-h)/h));
    if(INIT>=J || INIT<0){
        printf("Initial point wrong in 2D solver. INIT=%d, h=%g, sigma=%g\n",INIT,h,sigma);
        exit(EXIT_FAILURE);
    }
    
	if(d == 3)
		p[INIT] = 1.0/(4*rand_gen_constant_pi*r[INIT]*r[INIT]*h);
	else if(d == 2)
		p[INIT] = 1.0/(2*rand_gen_constant_pi*r[INIT]*h);
	else if(d == 1)
		p[INIT] = 1.0/h;
    
	for(int i = 0;i<J;i++) {
		p_prev[i] = p[i];
	}
    
	double *diag = (double *)malloc(J*sizeof(double));//[J];
	double *low= (double *)malloc((J-1)*sizeof(double));//[J-1];
	double *up= (double *)malloc((J-1)*sizeof(double));//[J-1];
	
    double Cconst = D*k/pow(h,2);
	diag[0] = 1+2*Cconst+(-Cconst+(d-1)*D*k/(2*h*sigma))/(1+h*k_r/const_d);//(1-D*k*(1/pow(h,2)*(-2+C1)-1/r[0]*C1*(d-1)/(2*h)));
    
	for(int i = 0;i<J-1;i++)
	{
		diag[i+1] = 1+2*Cconst;
		low[i] = -Cconst+(d-1)*D*k/(2*h*r[i+1]);//-D*k*(1.0/pow(h,2)-1/r[i+1]*(d-1)/(2*h));
		up[i] = -Cconst-(d-1)*D*k/(2*h*r[i]);//-D*k*(1.0/pow(h,2)+1/r[i]*(d-1)/(2*h));
	}
    //-Cconst-(d-1)*D/(2*h*r[i+1]);
    //Cconst+(d-1)*D/(2*h*r[i]);
	
	gsl_vector *dndiag = gsl_vector_alloc(J);
	gsl_vector *dnlow = gsl_vector_alloc(J-1);
	gsl_vector *dnup = gsl_vector_alloc(J-1);
	gsl_vector *gslp = gsl_vector_alloc(J);
    
	for(int i = 0;i<J-1;i++) {
		gsl_vector_set(dndiag,i,diag[i]);
		gsl_vector_set(dnup,i,up[i]);
		gsl_vector_set(dnlow,i,low[i]);
		gsl_vector_set(gslp,i,p[i]);
	}
	gsl_vector_set(dndiag,J-1,diag[J-1]);
	gsl_vector_set(gslp,J-1,p[J-1]);
    
	
	
    
	
	
	double *newp = (double *)malloc(J*sizeof(double));
	memcpy(newp,p,J*sizeof(double));
    
	
	
	double con1 = 0.0;
	double con2 = 0.0;
	double con3 = 0.0;
	for(int i = 0;i<T;i++)
	{
		gsl_linalg_solve_tridiag(dndiag,dnup,dnlow, gslp,gslp);
        
		for(int j = 0;j<J;j++) {
			p[j] = gsl_vector_get(gslp,j);
		}
		
		if(d==3){
			con1 = 4*rand_gen_constant_pi*sigma*sigma;
			con2 = 4*rand_gen_constant_pi*r[0]*r[0];
		}
		else if(d==2){
			con1 = 2*rand_gen_constant_pi*sigma;
			con2 = 2*rand_gen_constant_pi*r[0];
		}
		else if(d==1){
			con1 = 1;
			con2 = 1;
		}
		S[i] = (1/(1+(h*k_r)/const_d)*con1*p[0]+con2*p[0]);
		
		
		for(int j = 0;j<J-1;j++)
		{
			if(d == 3) {
				con2 = 4*rand_gen_constant_pi*r[j]*r[j];
				con3 = 4*rand_gen_constant_pi*r[j+1]*r[j+1];
			}
			else if(d == 2) {
				con2 = 2*rand_gen_constant_pi*r[j];
				con3 = 2*rand_gen_constant_pi*r[j+1];
			}
			else if(d == 1) {
				con2 = 1;
				con3 = 1;
			}
			S[i] += (con2*p[j]+con3*p[j+1]);
		}
		
		S[i] = 0.5*h*S[i];
        
		if(S[i]<1-ran && k_r>0)
		{
			double S_prev;
			if(i==0)
				S_prev = 1;
			else
				S_prev = S[i-1];
			
			(*t_reac) = i*k+(S[i]-(1-ran))/(S[i]-S_prev)*k;
			
			free(S);
			free(r);
			free(pdf);
			free(p);
			free(newp);
			free(p_prev);
			free(diag);
			free(up);
			free(low);
			
			gsl_vector_free(dndiag);
			gsl_vector_free(dnup);
			gsl_vector_free(dnlow);
			gsl_vector_free(gslp);
            
			return 0;
		}
		
	}
    
    
	
    
	for(int i = 0;i<J;i++)
	{
		if(d == 3)
			con2 = 4*rand_gen_constant_pi*r[i]*r[i];
		else if(d == 2)
			con2 = 2*rand_gen_constant_pi*r[i];
		else if(d == 1)
			con2 = 1;
		pdf[i] = con2*p[i];
	}
    
	ran = 1-ran;
	/* The case of no reaction. We determine r_new from the distribution function pdf. */
	(*t_reac) = t+1;
	double p_dist = 0.0;
	double p_prev2 = 0.0;
	int i = 1;
    
	p_dist = h/2*(1/(1+(h*k_r)/(const_d))*con1*p[0]+pdf[0]);
	while(p_dist<ran && i<J)
	{
		p_prev2 = p_dist;
		p_dist += h/2*(pdf[i]+pdf[i-1]);
		i++;
	}
	
	
	double r_new = r[i-1]-h+(ran-p_prev2)*h/(p_dist-p_prev2);
//	if(r_new-sigma<0.0) {
//        
//		printf("r_new = %g\n",r_new);
//	}
	free(pdf);
	free(S);
	free(r);
	free(p);
	free(newp);
	free(p_prev);
	free(diag);
	free(up);
	free(low);
	
	gsl_vector_free(dndiag);
	gsl_vector_free(dnup);
	gsl_vector_free(dnlow);
	gsl_vector_free(gslp);
    
	return r_new;
}


//double random_time_num(double k_r,double r_0,double sigma,
//                                 double D,double t,gsl_rng *rng,double *t_reac,int dimension,int time_steps,int space_steps,double ran)
//{
//	if(t!=t)
//	{
//		printf("t is nan\n");
//	}
//	//printf("r_0 = %g, sigma = %g, t = %g\n",r_0,sigma,t);
//    //	int space_steps = 40;
//	double K = 30;
//	//double ran = gsl_rng_uniform(rng);
//	
//	double k = t/(double)time_steps;
//	double h = (r_0-sigma)/space_steps;
//	
//	/*if(h<-1e-20)
//	 {
//	 printf("h = %e,  r_0 = %e,  sigma = %e\n",h,r_0,sigma);
//	 if(k_r>1e-30)
//	 {
//	 (*t_reac) = 1e-9;
//	 return 0;
//	 }
//	 (*t_reac) = t+1;
//	 return sigma+r_0;
//	 }*/
//	double endr = r_0+sqrt(2*D*t*K);
//	if(r_0+h>=endr)
//	{
//		endr = r_0+10*h;
//	}
//	if(r_0<sigma*(1+1/10000.0))
//	{
//		r_0 = sigma+sigma/10000.0;
//		//h = r_0-sigma;
//		h = (r_0-sigma)/space_steps;
//        //	printf("h = %e,  r_0 = %e,  sigma = %e,  k_r = %e,  t = %e\n",h,r_0,sigma,k_r,t);
//	}
//	
//	int J = (int)(round((endr-sigma)/h));
//	if(J<=1)
//	{
//		//printf("J<=1\n");
//		//h = (endr-r_0)/10.0;
//		endr = r_0+10*h;
//		J = (int)(round((endr-sigma)/h));
//		//	printf("New J = %d\n",J);
//	}
//	if(J>1e4)
//	{
//        	printf("J = %d\n",J);
//        
//		h = (r_0-sigma)/3;
//		J = (int)(round((endr-sigma)/h));
//		printf("New J = %d\n",J);
//		
//	}
//	
//	//printf("ran = %f\n",ran);
//	
//	double *pdf = (double *)malloc(J*sizeof(double));
//	double *S = (double *)malloc(time_steps*sizeof(double));
//	
//	
//	int d = dimension;
//	
//	double const_d = 0.0;
//	if(d==3)
//    const_d = 4*rand_gen_constant_pi*sigma*sigma*D;
//	else if (d==2)
//    const_d = 2*rand_gen_constant_pi*sigma*D;
//	else if (d==1)
//    const_d = D;
//	
//	double C1 = 1/(1+(h*k_r)/(const_d));
//	
//	/*if(J<=1)
//     {
//     printf("endr = %e\n",endr);
//     printf("r_0 = %e, h = %e, T = %d\n",r_0,h,T);
//     printf("J<=1\n");
//     }*/
//	double *r = (double *)malloc(J*sizeof(double));
//	for(int i = 0;i<J;i++)
//	{
//		r[i] = sigma+(i+1)*h;
//	}
//	
//	double t_init = pow(r_0-sigma,2)/(20*D);
//    
//	t_init = 0.0;
//    //	printf("t_init = %g\n",t_init);
//	int T = (int)round((t-t_init)/k);
//	double *p = (double *)calloc(J,sizeof(double));
//	double *p_prev = (double *)malloc(J*sizeof(double));
//    
//	int INIT = (int)(round((r_0-sigma-h)/h));
//    
//    /*
//     if(d==3) {
//     for(int i = 0;i<J;i++) {
//     p[i] = 1.0/(4*rand_gen_constant_pi*pow(r[i],2))*gsl_ran_gaussian_pdf(r[i]-h-r_0,sqrt(2*D*t));
//     }
//     }
//     else if(d==2) {
//     for(int i = 0;i<J;i++) {
//     p[i] = 1.0/(2*rand_gen_constant_pi*r[i])*gsl_ran_gaussian_pdf(r[i]-h-r_0,sqrt(2*D*t));
//     }
//     }
//     else if(d==1) {
//     for(int i = 0;i<J;i++) {
//     p[i] = gsl_ran_gaussian_pdf(r[i]-r_0,sqrt(2*D*t));
//     }
//     }*/
//    
//	if(d == 3)
//    p[INIT] = 1.0/(4*rand_gen_constant_pi*r[INIT]*r[INIT]*h);
//	else if(d == 2)
//    p[INIT] = 1.0/(2*rand_gen_constant_pi*r[INIT]*h);
//	else if(d == 1)
//    p[INIT] = 1.0/h;
//    
//	for(int i = 0;i<J;i++) {
//		p_prev[i] = p[i];
//	}
//    
//	double *diag = (double *)malloc(J*sizeof(double));//[J];
//	double *low= (double *)malloc((J-1)*sizeof(double));//[J-1];
//	double *up= (double *)malloc((J-1)*sizeof(double));//[J-1];
//	
//	diag[0] = (1-D*k*(1/pow(h,2)*(-2+C1)-1/r[0]*C1*(d-1)/(2*h)));
//	for(int i = 0;i<J-1;i++)
//	{
//		diag[i+1] = 1+D*k*2.0/pow(h,2);
//		low[i] = -D*k*(1.0/pow(h,2)-1/r[i+1]*(d-1)/(2*h));
//		up[i] = -D*k*(1.0/pow(h,2)+1/r[i]*(d-1)/(2*h));
//	}
//	
//	/*gsl_vector_view dndiag = gsl_vector_view_array(diag,J);
//     gsl_vector_view dnlow = gsl_vector_view_array(low,J-1);
//     gsl_vector_view dnup = gsl_vector_view_array(up,J-1);*/
//	
//	gsl_vector *dndiag = gsl_vector_alloc(J);
//	gsl_vector *dnlow = gsl_vector_alloc(J-1);
//	gsl_vector *dnup = gsl_vector_alloc(J-1);
//	gsl_vector *gslp = gsl_vector_alloc(J);
//	/*dndiag->data = diag;
//     dnup->data = up;
//     dnlow->data = low;*/
//	for(int i = 0;i<J-1;i++) {
//		gsl_vector_set(dndiag,i,diag[i]);
//		gsl_vector_set(dnup,i,up[i]);
//		gsl_vector_set(dnlow,i,low[i]);
//		gsl_vector_set(gslp,i,p[i]);
//	}
//	gsl_vector_set(dndiag,J-1,diag[J-1]);
//	gsl_vector_set(gslp,J-1,p[J-1]);
//	//gsl_vector_view gslprev = gsl_vector_view_array(p,J);
//	
//	
//	/*diag[J-1] = 0;*/
//	
//	
//	double *newp = (double *)malloc(J*sizeof(double));
//	memcpy(newp,p,J*sizeof(double));
//	//gsl_vector_view gslp = gsl_vector_view_array(newp,J);
//	
//	
//	double con1 = 0.0;
//	double con2 = 0.0;
//	double con3 = 0.0;
//	for(int i = 0;i<T;i++)
//	{
//		gsl_linalg_solve_tridiag(dndiag,dnup,dnlow, gslp,gslp);//(&dndiag.vector,&dnup.vector,&dnlow.vector, &gslp.vector, &gslp.vector);
//        
//		for(int j = 0;j<J;j++) {
//			p[j] = gsl_vector_get(gslp,j);
//		}
//		
//		
//        /*	double p_0 = 1/(1+(h*k_r)/(const_d))*p[0];
//         double pdf_0 = 2*rand_gen_constant_pi*sigma*p_0;// *sigma*p_0;
//         double simp = h/3.0*(pdf_0+4*2*rand_gen_constant_pi*r[0]*p[0]+2*rand_gen_constant_pi*r[1]*p[1]);
//         for(int j = 2;j<J-1;j+=2)
//         {
//         double temp1,temp2,temp3;
//         temp1 = 2*rand_gen_constant_pi*r[j-1];
//         temp2 = 2*rand_gen_constant_pi*r[j];
//         temp3 = 2*rand_gen_constant_pi*r[j+1];
//         simp += h/3.0*(temp1*p[j-1]+4*temp2*p[j]+temp3*p[j+1]);
//         }*/
//		/*
//		 if(fabs(simp-S[i])>1e-2) {
//		 printf("r_0 = %g, t = %g\n",r_0,t);
//		 }*/
//		
//		//printf("simp = %g\n",simp);
//		
//		
//		
//		
//		if(d==3){
//			con1 = 4*rand_gen_constant_pi*sigma*sigma;
//			con2 = 4*rand_gen_constant_pi*r[0]*r[0];
//		}
//		else if(d==2){
//			con1 = 2*rand_gen_constant_pi*sigma;
//			con2 = 2*rand_gen_constant_pi*r[0];
//		}
//		else if(d==1){
//			con1 = 1;
//			con2 = 1;
//		}
//		S[i] = (1/(1+(h*k_r)/(const_d))*con1*p[0]+con2*p[0]);
//        
//        
//		for(int j = 0;j<J-1;j++)
//		{
//			if(d == 3) {
//				con2 = 4*rand_gen_constant_pi*r[j]*r[j];
//				con3 = 4*rand_gen_constant_pi*r[j+1]*r[j+1];
//			}
//			else if(d == 2) {
//				con2 = 2*rand_gen_constant_pi*r[j];
//				con3 = 2*rand_gen_constant_pi*r[j+1];
//			}
//			else if(d == 1) {
//				con2 = 1;
//				con3 = 1;
//			}
//			S[i] += (con2*p[j]+con3*p[j+1]);
//		}
//		
//		S[i] = 0.5*h*S[i];
//        //	printf("S_%d = %g\n",i,S[i]);
//        
//		
//		/*if(r_0<=sigma+sigma/10000.0 && i==T-1){//fabs(S[i]-simp)/simp>0.002) {
//         printf("Error: %g\n",fabs(S[i]-simp)/simp);
//         printf("S_%d = %g, simp = %g\n",i,S[i],simp);
//         printf("t = %g\n",t);
//         }*/
//        //	S[i] = simp;
//        
//        
//		if(S[i]<1-ran && k_r>0)
//		{
//			double S_prev;
//			if(i==0)
//            S_prev = 1;
//			else
//            S_prev = S[i-1];
//			
//			(*t_reac) = i*k+(S[i]-(1-ran))/(S[i]-S_prev)*k;
//			
//			free(S);
//			free(r);
//			free(pdf);
//			free(p);
//			free(newp);
//			free(p_prev);
//			free(diag);
//			free(up);
//			free(low);
//			
//			gsl_vector_free(dndiag);
//			gsl_vector_free(dnup);
//			gsl_vector_free(dnlow);
//			gsl_vector_free(gslp);
//            
//			return 0;
//		}
//		
//	}
//	//double con2;
//    
//	
//    
//	for(int i = 0;i<J;i++)
//	{
//		if(d == 3)
//        con2 = 4*rand_gen_constant_pi*r[i]*r[i];
//		else if(d == 2)
//        con2 = 2*rand_gen_constant_pi*r[i];
//		else if(d == 1)
//        con2 = 1;
//		pdf[i] = con2*p[i];
//	}
//	/*double p_0 = 1/(1+(h*k_r)/(const_d))*p[0];
//     double pdf_0 = 4*rand_gen_constant_pi*sigma*sigma*p_0;
//     double simp = h/3.0*(pdf_0+4*pdf[0]+pdf[1]);
//     for(int i = 3;i<J-1;i+=2) {
//     simp += h/3.0*(pdf[i-1]+4*pdf[i]+pdf[i+1]);
//     }
//     printf("simp = %g\n",simp);
//     */
//	//printf("ran = %g\n",ran);
//	ran = 1-ran;
//	/* The case of no reaction. We determine r_new from the distribution function pdf. */
//	(*t_reac) = t+1;
//	double p_dist = 0.0;
//	double p_prev2 = 0.0;
//	int i = 1;
//    /*
//     if(d == 3)
//     p_dist = h/2*(1/(1+(h*k_r)/(const_d))*(4*rand_gen_constant_pi*sigma*sigma*p[0])+pdf[0]);
//     else if(d==2)
//     p_dist = h/2*(1/(1+(h*k_r)/(const_d))*(2*rand_gen_constant_pi*sigma*p[0])+pdf[0]);
//     else if(d==1)
//     p_dist = h/2*(1/(1+(h*k_r)/(const_d))*(p[0])+pdf[0]);
//     */
//	p_dist = h/2*(1/(1+(h*k_r)/(const_d))*con1*p[0]+pdf[0]);
//	while(p_dist<ran && i<J)
//	{
//		p_prev2 = p_dist;
//		p_dist += h/2*(pdf[i]+pdf[i-1]);
//		i++;
//	}
//	
//	
//	double r_new = r[i-1]-h+(ran-p_prev2)*h/(p_dist-p_prev2);
//	if(r_new-sigma<0.0) {
//        
//		printf("r_new = %g\n",r_new);
//	}
//	free(pdf);
//	free(S);
//	free(r);
//	free(p);
//	free(newp);
//	free(p_prev);
//	free(diag);
//	free(up);
//	free(low);
//	
//	gsl_vector_free(dndiag);
//	gsl_vector_free(dnup);
//	gsl_vector_free(dnlow);
//	gsl_vector_free(gslp);
//    
//	return r_new;
//}

double W2_real(double a,double b)
{
	
	if(a+b>4)
	{
		//printf("a=%10g,b=%g\n",a,b);
		double x = a+b;
		double temp = exp(-a*a)/(x*sqrt(rand_gen_constant_pi));
		double sum = 0;
		int i = 1;
		double term = 0;
        int NT = round(pow(a+b,2)-0.5);
        if(NT>4){
            NT = 5;
        }
		while(i<NT)
		{
            // printf("q = %.10g\n",fabs(term/sum));
            // printf("sum=%.10g\n",sum);
            //printf("i=%i\n",i);
			term = pow(-1,i)*gsl_sf_doublefact(2*i-1)/(pow(2*x*x,i));
			sum += term;
			i++;
		}
		return temp*(1+sum);
	}
	else
	{
        ////     //   printf(":::\n");
		double res = exp(2*a*b+b*b)*erfc(a+b);
		return res;
	}
}


double cdf_exact_refl(double t,double r_0,double r,double sigma,double D)
{
	double term1 = -(sqrt(D*t))/(r_0*sqrt(rand_gen_constant_pi))*(exp(-pow(r-r_0,2)/(4*D*t))-
                                                exp(-pow(r+r_0-2*sigma,2)/(4*D*t)))+
    0.5*(erf((r-r_0)/sqrt(4*D*t))+erf((r+r_0-2*sigma)/sqrt(4*D*t)));
	
	
	double alpha_irr = sqrt(D)/sigma;
	double K = (r+r_0-2*sigma)/sqrt(4*D*t);
	double term2 = 1/r_0*(r-sqrt(D)/alpha_irr)*W2_real(K,alpha_irr*sqrt(t));
	return term1-term2;
}


double p_irr(double t,double r_0,double r,double k_r,double sigma,double D){
    double a = 1.0/(4*rand_gen_constant_pi*r*r_0*sqrt(D));
    double b = 1.0/sqrt(4*rand_gen_constant_pi*t);
    double k_D = 4*rand_gen_constant_pi*sigma*D;
    double alpha_irr = (1+k_r/k_D)*sqrt(D)/sigma;
    double arg1 = -pow((r-r_0),2)/(4*D*t);
    double arg2 = -pow((r+r_0-2*sigma),2)/(4*D*t);
    double p_irr = b*(exp(arg1)+exp(arg2))-alpha_irr*W2_real((r+r_0-2*sigma)/sqrt(4*D*t),alpha_irr*sqrt(t));
    return a*p_irr;
}

double p_irr_helper(double h,void *params){
	double *temp = (double *)params;
    
	return 4*rand_gen_constant_pi*h*h*p_irr(temp[0],temp[1],h,temp[2],temp[3],temp[4]);
}

double cdf_exact_irr(double t,double r_0,double r,double k_r,
                     double sigma,double D)
{
	double term1 = -(sqrt(D*t))/(r_0*sqrt(rand_gen_constant_pi))*(exp(-pow(r-r_0,2)/(4*D*t))-exp(-pow(r+r_0-2*sigma,2)/(4*D*t)))+0.5*(erf((r-r_0)/sqrt(4*D*t))+erf((r+r_0-2*sigma)/sqrt(4*D*t)));
	//printf("term1_irr = %.15f\n",term1);
	double k_D = 4*rand_gen_constant_pi*sigma*D;
	double alpha_irr = (1+k_r/k_D)*sqrt(D)/sigma;
	//printf("alpha_irr = %.15f\n",alpha_irr);
	double K = (r+r_0-2*sigma)/sqrt(4*D*t);
	
	//double w2_test = W2_real(K,alpha_irr*sqrt(t));
	//complex temp = Complex(alpha_irr*sqrt(t),0.0);
	//complex w2_test2 = W2(K,&temp);
	//printf("W2_real = %.15f	 W2_complex = %.15f\n",w2_test,w2_test2.real);
	double term2 = 1/r_0*(r-sqrt(D)/alpha_irr)*W2_real(K,alpha_irr*sqrt(t));
	double term3 = k_r/(4*rand_gen_constant_pi*r_0*sigma*sqrt(D))*(-1.0/alpha_irr*W2_real((r_0-sigma)/sqrt(4*D*t),alpha_irr*sqrt(t)));
    //	printf("term2_irr = %.15f\n",term2);
	return term1-term2-term3;
}


double cdf_irr(double t,double r_0,double r,double k_r,double sigma,double D){
    gsl_set_error_handler_off();
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
  	double params[5] = {t,r_0,k_r,sigma,D};
	gsl_function F;
  	F.function = &p_irr_helper;
  	F.params = &params;
  	gsl_integration_qag(&F, sigma, r, 0, 1e-4, 1000,6,w, &result, &error);
    gsl_integration_workspace_free (w);
  	//printf("result=%.10g, error=%.10g\n",result,error);
  	return result;
}

double survivalA_irr(double t,double r_0,double k_r,double sigma,double D)
{
	double k_D = 4*rand_gen_constant_pi*sigma*D;
	double alpha_irr = (1+k_r/k_D)*sqrt(D)/sigma;
	//printf("x=%.15g\n",(r_0-sigma)/sqrt(4*D*t));
    double erfc_term = 0.0;
    if((r_0-sigma)/sqrt(4*D*t)>4){
        double x = (r_0-sigma)/sqrt(4*D*t);
        double sum = 0;
		int i = 1;
		double term = 0;
        int NT = round(pow(x,2)-0.5);
        if(NT>3){
            NT = 4;
        }
		while(i<NT)
		{
            // printf("q = %.10g\n",fabs(term/sum));
            // printf("sum=%.10g\n",sum);
            //printf("i=%i\n",i);
			term = pow(-1,i)*gsl_sf_doublefact(2*i-1)/(pow(2*x*x,i));
			sum += term;
			i++;
		}
        erfc_term = exp(-x*x)/(x*sqrt(rand_gen_constant_pi))*(1+sum);
    }
    else{
        erfc_term = erfc((r_0-sigma)/sqrt(4*D*t));
    }
    
	double c1 = erfc_term-W2_real((r_0-sigma)/sqrt(4*D*t),sqrt(t)*alpha_irr);
	return 1-(sigma/r_0)*(k_r/(k_r+k_D))*c1;
}

double survivalA_1D(double x_0,double k_r,double k_d,double t,double D)
{
	double K = x_0/sqrt(4*D*t);
	double Delta = sqrt(k_r*k_r-4*D*k_d);
	double lambda_plus = (k_r+Delta)/2.0;
	double lambda_minus = (k_r-Delta)/2.0;
	double sqdt = sqrt(t/D);
	double B_minus = lambda_minus*sqdt;
	double B_plus  = lambda_plus*sqdt;
	double k1=W2_real(K,B_minus);
	double k2=W2_real(K,B_plus);
	//printf("k1: %.16f k2: %.16f\n",k1,k2);
	//(W2_real(K,B_minus)-W2_real(K,B_plus))
	double temp = k_r/Delta*(k1-k2);
	if(k_r<1e-30)
	{
		temp = 0.0;
	}
	if(1-temp<0)
		printf("S<0\n");
	return 1-temp;
}

double random_time(double endt,double r_0,double k_r,double sigma,double D, gsl_rng *rn, int dimension)
{
	
	int J = 1000;
	int i = 0;
	double TOL = 1e-13;
	double Sexact = 0.0;
	if(dimension == 3)
	{
		Sexact = survivalA_irr(endt,r_0,k_r,sigma,D);//cdf_irr(endt,r_0,r_0+10*sqrt(2*D*endt),k_r,sigma,D);//
		/*if(k_r < 1e-30)
		 {
		 printf("Sexact = %.15f\n",Sexact);
		 }*/
	}
	else if(dimension == 1)
	{
		Sexact = survivalA_1D(r_0,k_r,0.0,endt,D);
	}
	else
	{
		printf("Invalid dimension.");
	}
	double r = gsl_rng_uniform(rn);
	if(r>(1-Sexact))
	{
		return endt+1;
	}
	double T = 0.0;
	
	double h = endt/4.0;
	double tt = endt/2.0;
	double ttprev = tt;
	
	double T_prev;
	
	
	while(fabs(r-(1.0-T))>TOL && i<J)
	{
		ttprev = tt;
		T_prev = T;
		if(dimension == 3)
			T = survivalA_irr(tt,r_0,k_r,sigma,D);//cdf_irr(tt,r_0,r_0+10*sqrt(2*D*tt),k_r,sigma,D);//
		else if(dimension == 1)
			T = survivalA_1D(r_0,k_r,0.0,tt,D);
		if(r-(1.0-T)<0)
		{
			tt -= h;
		}
		else
		{
			tt += h;
		}
		h = h/2.0;
		i++;
	}
	if(tt<(1e-10))
	{
		return 1e-10;
	}
	return tt;
	
	
	if(i==J)
	{
		return tt+(r-Sexact-T)/(T-T_prev)*h;
		return (tt+ttprev)/2.0;
	}
	
	return ttprev+h*((r-Sexact)-T_prev)/(T-T_prev);
}

double random_r_irr(double endr,double t,double r_0,double k_r,double sigma,
                              double D,gsl_rng *rn)
{
	double r;
	double R=0.0;
	double leftB,rightB;
	double TOL;
	int i;
	double surv_exact;
	double integ,F;
	
	
	leftB = sigma;
	rightB = endr;
	i = 0;
	TOL = 1e-10;//sigma/1e6;
	if(k_r<(1e-30))
	{
		surv_exact = 1;
	}
	else
	{
        //	surv_exact = cdf_irr(t,r_0,rightB,k_r,sigma,D);//
        
        surv_exact = survivalA_irr(t,r_0,k_r,sigma,D);
        /*  if(surv_exact<0.9){
         printf("surv_exact=%.10g, test=%.10g\n",surv_exact,test);
         }*/
        
	}
	r = gsl_rng_uniform(rn)*surv_exact;
	int counter = 0;
	while(fabs(rightB-leftB)/rightB>TOL && i<150)
	{
		R = (leftB+rightB)/2.0;
		if(k_r>(1e-30))
		{
            //SHOULD TRY THIS WITH CDF_IRR!!!
            
            
			integ = cdf_irr(t,r_0,R,k_r,sigma,D);
//            integ = cdf_exact_irr(t,r_0,R,k_r,sigma,D);
//            printf("integ=%.10g, exact=%.10g\n",integ,integ2);
            //printf("surv_exact=%.5g, integ=%.10g,  r=%.10g, iter=%d\n",surv_exact,integ,r,counter);
            //   printf("r_0=%.5g, r_new=%.5g\n",r_0,R);
		}
		else
		{
            //This case does not apply here.
			integ = cdf_exact_refl(t,r_0,R,sigma,D);//-(1-surv_exact);
            
		}
        counter++;
		F = integ;//-(1-surv_exact);
		
		//Intervall-halvering....
		if(F>r)
		{
			rightB = R;
		}
		else
		{
			leftB = R;
		}
		i++;
	}
	return R;
	
}

#endif
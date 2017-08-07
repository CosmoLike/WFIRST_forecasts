#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
/*
#include "../../theory/basics.c"
#include "../../theory/structs.c"
#include "../../theory/parameters.c"
#include "../../emu13/emu.c"
#include "../../theory/recompute.c"
#include "../../theory/cosmo3D.c"
#include "../../theory/redshift.c"
#include "../../theory/halo.c"
#include "../../theory/HOD.c"*/
#include "../../theory/GRS.c"
void init_GRS(int n_trade);

double log_like_GRS(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, 
	double B1, double B2, double B3, double B4, double B5, double B6, double B7, 
	 double SIGMAP1, double SIGMAP2, double SIGMAP3, double SIGMAP4, double SIGMAP5, double SIGMAP6, double SIGMAP7,
	 double SIGMAZ, double PSHOT, double KSTAR);
/*.......................................*/
/*.......................................*/
/*.......................................*/
/*.......................................*/
void fill_GRS_reference_cosmology(){
  if (GRS.N_z > 0 && GRS.N_z <= 10){
	int i;
	for (i = 0; i < GRS.N_z; i++){
		if (GRS.z[i]> 0. && GRS.z[i]< 1./limits.a_min-1.){
			GRS.H_ref[i] = cosmology.h0*hoverh0(1./(1+GRS.z[i]));
			GRS.DA_ref[i] = f_K(1./(1.+GRS.z[i]))/cosmology.h0;
		}
		else{
			printf("fill_GRS_reference_cosmology: GRZ.z[%d] = %e outside supported redshift range\nEXIT!\n",i,GRS.z[i]);
			exit(1);
		}	
	}
  }
  else{
	printf("fill_GRS_reference_cosmology: GRZ.N_z = %d outside supported range\nEXIT!\n",GRS.N_z);
	exit(1);
  }
}
void fill_GRS_default_parameters(){
	GRS.N_k = 25;
	GRS.N_mu = 10;
	// default value from Wang et al. 2013 (p.4)
	GRS.k_star = 0.24; //in h/Mpc
	GRS.k_min = 0.01; //in h/Mpc
	GRS.k_max = 0.3; //in h/Mpc
}

void init_GRS_WFIRST_trade1(){

	/*** now redshift bin related quantities***/
	GRS.N_z = 7;
	// representative redshift - density weighted over redshift bin
	GRS.z[0] = 0.72;
	GRS.z[1] = 1.00;
	GRS.z[2] = 1.34;
	GRS.z[3] = 1.67;
	GRS.z[4] = 2.01;
	GRS.z[5] = 2.29;
	GRS.z[6] = 2.57;
	// comoving survey volume in (Mpc/h)^3 for each redshift bin 
	GRS.V_z[0] = 1.181e+9;
	GRS.V_z[1] = 2.078e+9;
	GRS.V_z[2] = 2.588e+9;
	GRS.V_z[3] = 2.879e+9;
	GRS.V_z[4] = 2.583e+9;
	GRS.V_z[5] = 2.620e+9;
	GRS.V_z[6] = 2.614e+9;
	// galaxy parameters
	// galaxy density in (h/Mpc)^3
	GRS_gal.n_g[0] = 2.016e-3; 
	GRS_gal.n_g[1] = 2.603e-3; 
	GRS_gal.n_g[2] = 1.672e-3; 
	GRS_gal.n_g[3] = 7.560e-4; 
	GRS_gal.n_g[4] = 1.041e-4; 
	GRS_gal.n_g[5] = 3.247e-5; 
	GRS_gal.n_g[6] = 4.813e-6; 
	GRS_gal.b_g[0] = 1.19;
	GRS_gal.b_g[1] = 1.30;
	GRS_gal.b_g[2] = 1.44;
	GRS_gal.b_g[3] = 1.57;
	GRS_gal.b_g[4] = 1.70;
	GRS_gal.b_g[5] = 1.82;
	GRS_gal.b_g[6] = 1.93;
	for (int i = 0; i < GRS.N_z; i++){
		GRS_gal.sigma_p[i] = 290.0; // in km/s
		GRS_gal.sigma_z[i] = 1.e-3; // fractional accuracy
		GRS_gal.P_shot[i] = 0.0; // in (Mpc/h)^3
	}
	fill_GRS_default_parameters();
	fill_GRS_reference_cosmology();

}

void init_GRS_WFIRST_trade2(){

	/*** now redshift bin related quantities***/
	GRS.N_z = 7;
	// representative redshift - density weighted over redshift bin
	GRS.z[0] = 0.72;
	GRS.z[1] = 1.00;
	GRS.z[2] = 1.34;
	GRS.z[3] = 1.65;
	GRS.z[4] = 1.96;
	GRS.z[5] = 2.25;
	GRS.z[6] = 2.51;
	// comoving survey volume in (Mpc/h)^3 for each redshift bin 
	GRS.V_z[0] = 1.189e+9;
	GRS.V_z[1] = 2.078e+9;
	GRS.V_z[2] = 2.588e+9;
	GRS.V_z[3] = 2.457e+9;
	GRS.V_z[4] = 2.571e+9;
	GRS.V_z[5] = 2.617e+9;
	GRS.V_z[6] = 2.183e+9;
	// galaxy parameters
	// galaxy density in (h/Mpc)^3
	GRS_gal.n_g[0] = 2.551e-3; 
	GRS_gal.n_g[1] = 2.973e-3; 
	GRS_gal.n_g[2] = 1.947e-3; 
	GRS_gal.n_g[3] = 9.839e-4; 
	GRS_gal.n_g[4] = 1.485e-4; 
	GRS_gal.n_g[5] = 5.684e-5; 
	GRS_gal.n_g[6] = 1.214e-5; 
	GRS_gal.b_g[0] = 1.19;
	GRS_gal.b_g[1] = 1.30;
	GRS_gal.b_g[2] = 1.44;
	GRS_gal.b_g[3] = 1.56;
	GRS_gal.b_g[4] = 1.69;
	GRS_gal.b_g[5] = 1.80;
	GRS_gal.b_g[6] = 1.91;
	for (int i = 0; i < GRS.N_z; i++){
		GRS_gal.sigma_p[i] = 290.0; // in km/s
		GRS_gal.sigma_z[i] = 1.e-3; // fractional accuracy
		GRS_gal.P_shot[i] = 0.0; // in (Mpc/h)^3
	}
	fill_GRS_default_parameters();
	fill_GRS_reference_cosmology();

}
void init_GRS_WFIRST_trade3(){

	/*** now redshift bin related quantities***/
	GRS.N_z = 7;
	// representative redshift - density weighted over redshift bin
	GRS.z[0] = 0.73;
	GRS.z[1] = 1.00;
	GRS.z[2] = 1.34;
	GRS.z[3] = 1.65;
	GRS.z[4] = 1.96;
	GRS.z[5] = 2.23;
	GRS.z[6] = 2.49;
	// comoving survey volume in (Mpc/h)^3 for each redshift bin 
	GRS.V_z[0] = 2.059e+9;
	GRS.V_z[1] = 4.157e+9;
	GRS.V_z[2] = 5.176e+9;
	GRS.V_z[3] = 4.913e+9;
	GRS.V_z[4] = 5.142e+9;
	GRS.V_z[5] = 5.235e+9;
	GRS.V_z[6] = 2.622e+9;
	// galaxy parameters
	// galaxy density in (h/Mpc)^3
	GRS_gal.n_g[0] = 6.760e-4; 
	GRS_gal.n_g[1] = 8.166e-4; 
	GRS_gal.n_g[2] = 4.411e-4; 
	GRS_gal.n_g[3] = 1.717e-4; 
	GRS_gal.n_g[4] = 1.143e-5; 
	GRS_gal.n_g[5] = 2.425e-6; 
	GRS_gal.n_g[6] = 3.397e-7; 
	GRS_gal.b_g[0] = 1.19;
	GRS_gal.b_g[1] = 1.30;
	GRS_gal.b_g[2] = 1.44;
	GRS_gal.b_g[3] = 1.56;
	GRS_gal.b_g[4] = 1.68;
	GRS_gal.b_g[5] = 1.79;
	GRS_gal.b_g[6] = 1.90;
	for (int i = 0; i < GRS.N_z; i++){
		GRS_gal.sigma_p[i] = 290.0; // in km/s
		GRS_gal.sigma_z[i] = 1.e-3; // fractional accuracy
		GRS_gal.P_shot[i] = 0.0; // in (Mpc/h)^3
	}
	fill_GRS_default_parameters();
	fill_GRS_reference_cosmology();
}

void init_GRS_WFIRST_trade4(){

	/*** now redshift bin related quantities***/
	GRS.N_z = 7;
	// representative redshift - density weighted over redshift bin
	GRS.z[0] = 0.72;
	GRS.z[1] = 1.00;
	GRS.z[2] = 1.34;
	GRS.z[3] = 1.66;
	GRS.z[4] = 1.97;
	GRS.z[5] = 2.25;
	GRS.z[6] = 2.52;
	// comoving survey volume in (Mpc/h)^3 for each redshift bin 
	GRS.V_z[0] = 5.148e+8;
	GRS.V_z[1] = 1.039e+9;
	GRS.V_z[2] = 1.294e+9;
	GRS.V_z[3] = 1.228e+9;
	GRS.V_z[4] = 1.286e+9;
	GRS.V_z[5] = 1.309e+9;
	GRS.V_z[6] = 1.091e+9;
	// galaxy parameters
	// galaxy density in (h/Mpc)^3
	GRS_gal.n_g[0] = 4.907e-3; 
	GRS_gal.n_g[1] = 5.605e-3; 
	GRS_gal.n_g[2] = 4.013e-3; 
	GRS_gal.n_g[3] = 2.303e-4; 
	GRS_gal.n_g[4] = 4.839e-4; 
	GRS_gal.n_g[5] = 2.360e-4; 
	GRS_gal.n_g[6] = 7.448e-5; 
	GRS_gal.b_g[0] = 1.19;
	GRS_gal.b_g[1] = 1.30;
	GRS_gal.b_g[2] = 1.44;
	GRS_gal.b_g[3] = 1.56;
	GRS_gal.b_g[4] = 1.69;
	GRS_gal.b_g[5] = 1.80;
	GRS_gal.b_g[6] = 1.91;
	for (int i = 0; i < GRS.N_z; i++){
		GRS_gal.sigma_p[i] = 290.0; // in km/s
		GRS_gal.sigma_z[i] = 1.e-3; // fractional accuracy
		GRS_gal.P_shot[i] = 0.0; // in (Mpc/h)^3
	}
	fill_GRS_default_parameters();
	fill_GRS_reference_cosmology();

}
int set_cosmology_params_GRS(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=S8;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  cosmology.h0=H0;

  if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
  if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  return 1;
}


int set_nuisance_GRS_gbias(double B1, double B2, double B3, double B4, double B5, double B6, double B7)
{
  if (B1 == -42.){return 1;}
  GRS_gal.b_g[0] = B1;
  GRS_gal.b_g[1] = B2;
  GRS_gal.b_g[2] = B3;
  GRS_gal.b_g[3] = B4;
  GRS_gal.b_g[4] = B5;
  GRS_gal.b_g[5] = B6;
  GRS_gal.b_g[6] = B7;
//  GRS_gal.b_g[7] = B8;
//  GRS_gal.b_g[8] = B9;
//  GRS_gal.b_g[9] = B10;
  for (int i = 0; i < GRS.N_z; i++){
    if (GRS_gal.b_g[i] < 1.0 || GRS_gal.b_g[i] > 2.0) return 0;
  }
  return 1;
} 

int set_nuisance_GRS_sigmas(double SIGMAZ, double SIGMAP1,double SIGMAP2,double SIGMAP3,double SIGMAP4,double SIGMAP5,double SIGMAP6,double SIGMAP7)
{
  if (SIGMAP1 == -42.){return 1;}
  int i;
   if (SIGMAZ < 1.e-06 || SIGMAZ > 0.1) return 0;
   GRS_gal.sigma_p[0] = SIGMAP1;
   GRS_gal.sigma_p[1] = SIGMAP2;
   GRS_gal.sigma_p[2] = SIGMAP3;
   GRS_gal.sigma_p[3] = SIGMAP4;
   GRS_gal.sigma_p[4] = SIGMAP5;
   GRS_gal.sigma_p[5] = SIGMAP6;
   GRS_gal.sigma_p[6] = SIGMAP7;
   for (i = 0; i < GRS.N_z; i++){
  		GRS_gal.sigma_z[i] = SIGMAZ;
   		if (GRS_gal.sigma_p[i] < 150.0 || GRS_gal.sigma_p[i] > 500.0) return 0;
   }
  return 1;
} 

int set_nuisance_GRS_pshot(double PSHOT)
{
  if (PSHOT == -42.){return 1;}	
  for (int i = 0; i < GRS.N_z; i++){
  	GRS_gal.P_shot[i] = PSHOT;
    if (GRS_gal.P_shot[i] < -1.e3 || GRS_gal.P_shot[i] > 1.e3) return 0;
  }
  return 1;
} 
int set_nuisance_GRS_kstar(double kstar)
{
  if (kstar == -42.){return 1;}	
  GRS.k_star=kstar;
  if (kstar < 0.01 || kstar > 1.) return 0;
  return 1;
} 

void set_data_GRS(double *k, double *mu,double *data){
	int nz,nk,nm,i;
	for (nz = 0; nz < GRS.N_z; nz++){
		for (nk = 0; nk < GRS.N_k; nk++){
			for(nm = 0; nm < GRS.N_mu; nm ++){
				i = nz*GRS.N_k*GRS.N_mu + nk*GRS.N_mu +nm;
				data[i] = P_obs(k[nk],mu[nm],nz);
				//printf("%d %e %e %e\n",i, k[nk],mu[nm],data[i]);
			}
		}
	}
}
//Gaussian 3D power spectrum "covariance" = variance, assuming linear binning in k and mu
void set_variance_GRS(double *k, double *mu,double dk, double dmu, double *var){
	int nz,nk,nm,i;
	double P, n, SN =0;
	for (nz = 0; nz < GRS.N_z; nz++){
		n = GRS_gal.n_g[nz];
		for (nk = 0; nk < GRS.N_k; nk++){
			for(nm = 0; nm < GRS.N_mu; nm ++){
				i = nz*GRS.N_k*GRS.N_mu + nk*GRS.N_mu +nm;
				P = P_obs(k[nk],mu[nm],nz);
				var[i] = pow(P+1./n,2.0)
					*1./(pow(n*P/(1.+n*P),2.0)*GRS.V_z[nz])
					*1./(k[nk]*k[nk]/pow(2.*M_PI,2.0)*dk*dmu);
				//SN +=P*P/var[i];
			}
		}
	}
//	printf("total S/N: %e\n",SN);
}

void init_GRS(int nt){
	set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
	switch(nt){
		case 1:
			printf("initialize GRS trade study case 1\n");
			init_GRS_WFIRST_trade1();
			break;
		case 2:
			printf("initialize GRS trade study case 2\n");
			init_GRS_WFIRST_trade2();
			break;
		case 3:
			printf("initialize GRS trade study case 3 (wide)\n");
			init_GRS_WFIRST_trade3();
			break;
		case 4:
			printf("initialize GRS trade study case 4 (deep)\n");
			init_GRS_WFIRST_trade4();
			break;
		default:
			printf("like_grs.c:init: GRS trade study number %d not specified\nEXIT\n",nt);
			exit(1);
	}
	double dk, dmu;
	int i;
	GRS.k = create_double_vector(0,GRS.N_k-1);
	GRS.mu = create_double_vector(0,GRS.N_mu-1);
	GRS.datav = create_double_vector(0,GRS.N_z*GRS.N_k*GRS.N_mu-1);
	GRS.var = create_double_vector(0,GRS.N_z*GRS.N_k*GRS.N_mu-1);
	dk = (GRS.k_max-GRS.k_min)/GRS.N_k;
	dmu= 1./GRS.N_mu;
	for (i =0; i < GRS.N_k; i++){
		// use linear spacing in k!
		GRS.k[i] = GRS.k_min +(i+0.5)*dk;
	}
	for (i =0; i < GRS.N_mu; i++){
		// use linear spacing in mu!
		GRS.mu[i] = (i+0.5)*dmu;
	}
	set_data_GRS(GRS.k,GRS.mu,GRS.datav);
	set_variance_GRS(GRS.k,GRS.mu, dk, dmu, GRS.var);
	like.GRS = 1;
	printf("GRS initalized\n");
}


double log_like_GRS(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, 
	double B1, double B2, double B3, double B4, double B5, double B6, double B7, 
	 double SIGMAP1, double SIGMAP2, double SIGMAP3, double SIGMAP4, double SIGMAP5, double SIGMAP6, double SIGMAP7,
	 double SIGMAZ, double PSHOT, double KSTAR)
{
//	printf("%le %le %le %le\n",B1, B2, B3, B4);
//	printf("%le %le %le %le\n",OMM, S8, NS, W0);
  int i,j,k,m=0,l;
  static double *pred;
  if (pred==0) pred= create_double_vector(0, GRS.N_z*GRS.N_k*GRS.N_mu-1);
  double chisqr,a,log_L_prior=0.0;
  
  if (set_cosmology_params_GRS(OMM,S8,NS,W0,WA,OMB,H0)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e15;
  }
  if (set_nuisance_GRS_gbias(B1,B2,B3,B4,B5,B6,B7)==0){
    printf("GRS Gbias out of range\n");
    return -1.0e15;
  } 
  if (set_nuisance_GRS_sigmas(SIGMAZ,SIGMAP1,SIGMAP2,SIGMAP3,SIGMAP4,SIGMAP5,SIGMAP6,SIGMAP7)==0){
    printf("GRS sigmas out of range\n");
    return -1.0e15;
  } 
  if (set_nuisance_GRS_pshot(PSHOT)==0){
    printf("GRS PSHOT out of range\n");
    return -1.0e15;
  } 
  if (set_nuisance_GRS_kstar(KSTAR)==0){
    printf("GRS KSTAR out of range\n");
    return -1.0e15;
  } 

  log_L_prior=0.0;
  // add Gaussian prior on SIGMAZ
  log_L_prior -= 0.5*pow((SIGMAZ-1.e-3)/1.e-4,2.);
  // add Gaussian prior on k_star
  log_L_prior -= 0.5*pow((KSTAR-0.24)/0.024,2.);

  set_data_GRS(GRS.k,GRS.mu,pred);

  chisqr=0.0;
  for (i=0; i<GRS.N_z*GRS.N_k*GRS.N_mu; i++){
//  	printf("%i %e %e\n",i,pred[i],GRS.datav[i]);
    a=(pred[i]-GRS.datav[i])*(pred[i]-GRS.datav[i])/GRS.var[i];
    chisqr=chisqr+a;
  }
  if (chisqr<0.0){
    printf("errror: chisqr < 0\n");
  }
  if (chisqr<-1.0) exit(EXIT_FAILURE);
  
//	printf("%le %e\n",-0.5*chisqr,log_L_prior);
  return -0.5*chisqr+log_L_prior;
}


// int main(void){
// 	init(1);
// //	double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, 
// //	double B1, double B2, double B3, double B4, double B5, double B6, double B7, 
// //	 double SIGMAP1, double SIGMAP2, double SIGMAP3, double SIGMAP4, double SIGMAP5, double SIGMAP6, double SIGMAP7,
// //	 double SIGMAZ, double PSHOT, double KSTAR)

// 	log_multi_like(0.315,0.831,0.965,-1.0,0.0,0.049,0.673,1.19,1.30,1.44,1.57,1.70,1.82,1.93,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.094);
// 	double dk, dmu, *k, *mu, *data, *variance;
// 	int i;
// 	k = create_double_vector(0,GRS.N_k-1);
// 	mu = create_double_vector(0,GRS.N_mu-1);
// 	data = create_double_vector(0,GRS.N_z*GRS.N_k*GRS.N_mu-1);
// 	variance = create_double_vector(0,GRS.N_z*GRS.N_k*GRS.N_mu-1);
// 	dk = (GRS.k_max-GRS.k_min)/GRS.N_k;
// 	dmu= 1./GRS.N_mu;
// 	for (i =0; i < GRS.N_k; i++){
// 		// use linear spacing in k!
// 		k[i] = GRS.k_min +(i+0.5)*dk;
// 	}
// 	for (i =0; i < GRS.N_mu; i++){
// 		// use linear spacing in mu!
// 		mu[i] = (i+0.5)*dmu;
// 	}
// 	set_data_GRS(k,mu,data);
// 	set_variance_GRS(k,mu, dk, dmu, variance);
// 	return 0;
// }

/*int main(void){
	set_cosmological_parameters_to_Planck_WP();
	init_GRS_single_zbin(0.16,0.4);
	double z;
	GRS_gal.b_g[0] = 1.0;
	for (z = 0.1; z < 2.5; z+=0.05){
		GRS.z[0] = z;
		printf("%.3f  %e %e %e\n", z, 1./P_obs(0.14,0.6,0), 1./P_obs(0.5,0.6,0),1./P_obs(1.0,0.6,0));
	}
	return 0;
}
*/
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <fftw3.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
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
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/pt.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "../cosmolike_core/theory/GRS.c"
#include "like_grs.c"


double log_like_wrapper(input_cosmo_params ic, input_nuisance_params_grs in);

double log_like_wrapper(input_cosmo_params ic, input_nuisance_params_grs in)
{
	// printf("%le %le %le %le %le %le %le\n",ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0);
	// printf("%le %le %le %le %le %le %le\n",in.grsbias[0], in.grsbias[1], in.grsbias[2], in.grsbias[3],in.grsbias[4], in.grsbias[5], in.grsbias[6]);
	//printf("%le %le %le %le %le %le %le %le %le %le\n",in.grssigmap[0], in.grssigmap[1], in.grssigmap[2], in.grssigmap[3], in.grssigmap[4], in.grssigmap[5], in.grssigmap[6],in.grssigmaz, in.grspshot, in.grskstar);
  double like = log_like_GRS(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,in.grsbias[0], in.grsbias[1], in.grsbias[2], in.grsbias[3],in.grsbias[4], in.grsbias[5], in.grsbias[6], in.grssigmap[0], in.grssigmap[1], in.grssigmap[2], in.grssigmap[3], in.grssigmap[4], in.grssigmap[5], in.grssigmap[6],in.grssigmaz, in.grspshot, in.grskstar);
  
  return like;
}

int main(void){
	double res;
	init_GRS();
	res=log_like_GRS(0.3156,0.831,0.9645,-1.0,0.0,0.0491685,0.6727,0.0,0.0,1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,2.975011712138650,3.376705680190931,3.725882076395691,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.24);
	printf("loglike=%le\n",res);
	res=log_like_GRS(0.2325384451023471,1.0310887691615935,0.8595022066924793, -1.137588053600136,  0.14993639190238117,  0.05165678142188447, 0.8247432633482439,0.0,0.0,1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,2.975011712138650,3.376705680190931,3.725882076395691,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.24);
	printf("loglike=%le\n",res);
return 0;
}	

//	double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, 
//	double B1, double B2, double B3, double B4, double B5, double B6, double B7, 
//	 double SIGMAP1, double SIGMAP2, double SIGMAP3, double SIGMAP4, double SIGMAP5, double SIGMAP6, double SIGMAP7,
//	 double SIGMAZ, double PSHOT, double KSTAR)

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
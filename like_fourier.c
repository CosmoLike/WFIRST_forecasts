#include <math.h>
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


#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/HOD.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"
#include "../cosmolike_core/theory/IA.c"
#include "../cosmolike_core/theory/cluster.c"
#include "../cosmolike_core/theory/BAO.c"
#include "../cosmolike_core/theory/external_prior.c"
#include "init_WFIRST_forecasts.c"
#include "like_grs.c"

double C_shear_tomo_sys(double ell,int z1,int z2);
double C_cgl_tomo_sys(double ell_Cluster,int zl,int nN, int zs);
double C_gl_tomo_sys(double ell,int zl,int zs);
void set_data_shear(int Ncl, double *ell, double *data, int start);
void set_data_ggl(int Ncl, double *ell, double *data, int start);
void set_data_clustering(int Ncl, double *ell, double *data, int start);
void set_data_cluster_N(double *data, int start);
void set_data_cgl(double *ell_Cluster, double *data, int start);
void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q,double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope);
double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q,double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope);
void write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in);
double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in);
int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_ell(void);

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}
int get_N_ell(void){
  return like.Ncl;
}


double C_shear_tomo_sys(double ell, int z1, int z2)
{
  double C;
  // C= C_shear_tomo_nointerp(ell,z1,z2);
  // if(like.IA==1) C+=C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  
  if(like.IA!=1) C= C_shear_tomo_nointerp(ell,z1,z2);
  //if(like.IA==1) C= C_shear_shear_IA(ell,z1,z2);
  if(like.IA==1) C = C_shear_tomo_nointerp(ell,z1,z2)+C_II_nointerp(ell,z1,z2)+C_GI_nointerp(ell,z1,z2);
  if(like.IA==2) C += C_II_lin_nointerp(ell,z1,z2)+C_GI_lin_nointerp(ell,z1,z2);  
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
  //printf("%le %d %d %le\n",ell,z1,z2,C_shear_tomo_nointerp(ell,z1,z2)+C_II_JB_nointerp(ell,z1,z2)+C_GI_JB_nointerp(ell,z1,z2));
return C;
}

double C_gl_tomo_sys(double ell,int zl,int zs)
{
  double C;
  // C=C_gl_tomo_nointerp(ell,zl,zs); 
  // if(like.IA==1) C += C_gI_nointerp(ell,zl,zs);
  
  if(like.IA!=1) C=C_gl_tomo_nointerp(ell,zl,zs);
  if(like.IA==1) C = C_ggl_IA(ell,zl,zs);
  if(like.IA==2) C += C_gI_lin_nointerp(ell,zl,zs);
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}

double C_cgl_tomo_sys(double ell_Cluster, int zl,int nN, int zs)
{
  double C;
  C=C_cgl_tomo_nointerp(ell_Cluster,zl,nN,zs);
  //if(like.IA!=0) C += 
  if(like.shearcalib==1) C *=(1.0+nuisance.shear_calibration_m[zs]);
return C;
}      

void set_data_shear(int Ncl, double *ell, double *data, int start)
{
  int i,z1,z2,nz;
  double a;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < Ncl; i++){
      if (ell[i] < like.lmax_shear){ data[Ncl*nz+i] = C_shear_tomo_sys(ell[i],z1,z2);}
      else {data[Ncl*nz+i] = 0.;}
    }
  }
}

void set_data_ggl(int Ncl, double *ell, double *data, int start)
{
  int i, zl,zs,nz;  
  for (nz = 0; nz < tomo.ggl_Npowerspectra; nz++){
    zl = ZL(nz); zs = ZS(nz);
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],zl)){
        data[start+(Ncl*nz)+i] = C_gl_tomo_sys(ell[i],zl,zs);
      }
      else{
        data[start+(Ncl*nz)+i] = 0.;
      }
    } 
  }
}

void set_data_clustering(int Ncl, double *ell, double *data, int start){
  int i, nz;
  for (nz = 0; nz < tomo.clustering_Npowerspectra; nz++){
    //printf("%d %e %e\n",nz, gbias.b[nz][1],pf_photoz(gbias.b[nz][1],nz));
    for (i = 0; i < Ncl; i++){
      if (test_kmax(ell[i],nz)){data[start+(Ncl*nz)+i] = C_cl_tomo_nointerp(ell[i],nz,nz);}
      else{data[start+(Ncl*nz)+i] = 0.;}
      //printf("%d %d %le %le\n",nz,nz,ell[i],data[Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + nz)+i]);
    }
  }
}

void set_data_cluster_N(double *data, int start){
  int nN, nz;
  for (nz = 0; nz < tomo.cluster_Nbin; nz++){
    for (nN = 0; nN < Cluster.N200_Nbin; nN++){
      data[start+Cluster.N200_Nbin*nz+nN] = N_N200(nz, nN);
    }
  }
}


void set_data_cgl(double *ell_Cluster, double *data, int start)
{
  int zl,zs,nN,nz,i,j;
  for(nN = 0; nN < Cluster.N200_Nbin; nN++){
    for (nz = 0; nz < tomo.cgl_Npowerspectra; nz++){
      zl = ZC(nz); zs = ZSC(nz);
      for (i = 0; i < Cluster.lbin; i++){
        j = start;
        j += (nz*Cluster.N200_Nbin+nN)*Cluster.lbin +i;
        data[j] = C_cgl_tomo_sys(ell_Cluster[i],zl,nN,zs);
      }
    }
  }
}


int set_cosmology_params(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  cosmology.sigma_8=S8;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  cosmology.h0=H0;
  cosmology.MGSigma=MGSigma;
  cosmology.MGmu=MGmu;

  if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0;
  if (cosmology.omb < 0.04 || cosmology.omb > 0.055) return 0;
  if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0;
  if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0;
  if (cosmology.w0 < -2.1 || cosmology.w0 > -0.0) return 0;
  if (cosmology.wa < -2.6 || cosmology.wa > 2.6) return 0;
  if (cosmology.h0 < 0.4 || cosmology.h0 > 0.9) return 0;
  //CH BEGINS 
  //CH: to use for running planck15_BA0_w0_wa prior alone) 
  //printf("like_fourier.c from WFIRST_forecasts: cosmology bounds set for running with planck15_BA0_w0_wa prior\n");
  //if (cosmology.Omega_m < 0.05 || cosmology.Omega_m > 0.6) return 0; 
  //if (cosmology.omb < 0.01 || cosmology.omb > 0.1) return 0; 
  //if (cosmology.sigma_8 < 0.5 || cosmology.sigma_8 > 1.1) return 0; 
  //if (cosmology.n_spec < 0.84 || cosmology.n_spec > 1.06) return 0; 
  //if (cosmology.w0 < -2.1 || cosmology.w0 > 1.5) return 0; 
  //if (cosmology.wa < -5.0 || cosmology.wa > 2.6) return 0; 
  //if (cosmology.h0 < 0.3 || cosmology.h0 > 0.9) return 0; 
  //CH ENDS
  return 1;
}

void set_nuisance_shear_calib(double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10)
{
  nuisance.shear_calibration_m[0] = M1;
  nuisance.shear_calibration_m[1] = M2;
  nuisance.shear_calibration_m[2] = M3;
  nuisance.shear_calibration_m[3] = M4;
  nuisance.shear_calibration_m[4] = M5;
  nuisance.shear_calibration_m[5] = M6;
  nuisance.shear_calibration_m[6] = M7;
  nuisance.shear_calibration_m[7] = M8;
  nuisance.shear_calibration_m[8] = M9;
  nuisance.shear_calibration_m[9] = M10;
}

int set_nuisance_shear_photoz(double SP1,double SP2,double SP3,double SP4,double SP5,double SP6,double SP7,double SP8,double SP9,double SP10,double SPS1)
{
  int i;
  nuisance.bias_zphot_shear[0]=SP1;
  nuisance.bias_zphot_shear[1]=SP2;
  nuisance.bias_zphot_shear[2]=SP3;
  nuisance.bias_zphot_shear[3]=SP4;
  nuisance.bias_zphot_shear[4]=SP5;
  nuisance.bias_zphot_shear[5]=SP6;
  nuisance.bias_zphot_shear[6]=SP7;
  nuisance.bias_zphot_shear[7]=SP8;
  nuisance.bias_zphot_shear[8]=SP9;
  nuisance.bias_zphot_shear[9]=SP10;
  
  for (i=0;i<tomo.shear_Nbin; i++){ 
    nuisance.sigma_zphot_shear[i]=SPS1;
    if (nuisance.sigma_zphot_shear[i]<0.001) return 0;
  }
  return 1;
}

int set_nuisance_clustering_photoz(double CP1,double CP2,double CP3,double CP4,double CP5,double CP6,double CP7,double CP8,double CP9,double CP10,double CPS1)
{
  int i;
  nuisance.bias_zphot_clustering[0]=CP1;
  nuisance.bias_zphot_clustering[1]=CP2;
  nuisance.bias_zphot_clustering[2]=CP3;
  nuisance.bias_zphot_clustering[3]=CP4;
  nuisance.bias_zphot_clustering[4]=CP5;
  nuisance.bias_zphot_clustering[5]=CP6;
  nuisance.bias_zphot_clustering[6]=CP7;
  nuisance.bias_zphot_clustering[7]=CP8;
  nuisance.bias_zphot_clustering[8]=CP9;
  nuisance.bias_zphot_clustering[9]=CP10;
  
  for (i=0;i<tomo.clustering_Nbin; i++){ 
    nuisance.sigma_zphot_clustering[i]=CPS1;
    if (nuisance.sigma_zphot_clustering[i]<0.001) return 0;
  }
  return 1;
}

int set_nuisance_ia(double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q)
{
  nuisance.A_ia=A_ia;  
  nuisance.beta_ia=beta_ia;
  nuisance.eta_ia=eta_ia;
  nuisance.eta_ia_highz=eta_ia_highz;
  nuisance.LF_alpha=LF_alpha;
  nuisance.LF_P=LF_P;
  nuisance.LF_Q=LF_Q;
  nuisance.LF_red_alpha=LF_red_alpha;
  nuisance.LF_red_P=LF_red_P;
  nuisance.LF_red_Q=LF_red_Q;
  if (nuisance.A_ia < 0.0 || nuisance.A_ia > 10.0) return 0;
  if (nuisance.beta_ia < -1.0 || nuisance.beta_ia > 3.0) return 0;
  if (nuisance.eta_ia < -3.0 || nuisance.eta_ia> 3.0) return 0;
  if (nuisance.eta_ia_highz < -1.0 || nuisance.eta_ia_highz> 1.0) return 0;
  // if(like.IA!=0){
  //  if (check_LF()) return 0;
  // }
return 1;
}

int set_nuisance_cluster_Mobs(double cluster_Mobs_lgN0,  double cluster_Mobs_alpha, double cluster_Mobs_beta, double cluster_Mobs_sigma0, double cluster_Mobs_sigma_qm, double cluster_Mobs_sigma_qz)
{
  //  nuisance.cluster_Mobs_lgM0 = mass_obs_norm;  //fiducial : 1.72+log(1.e+14*0.7); could use e.g. sigma = 0.2 Gaussian prior
  //  nuisance.cluster_Mobs_alpha = mass_obs_slope; //fiducial: 1.08; e.g. sigma = 0.1 Gaussian prior
  //  nuisance.cluster_Mobs_beta = mass_z_slope; //fiducial: 0.0; e.g. sigma = 0.1 Gaussian prior
  //  nuisance.cluster_Mobs_sigma = mass_obs_scatter; //fiducial 0.25; e.g. sigma = 0.05 Gaussian prior

  // fiducial values and priors from Murata et al. (2018) except for redshift-related parameters
  nuisance.cluster_Mobs_lgN0 = cluster_Mobs_lgN0; //fiducial: 3.207, flat prior [0.5, 5.0]
  nuisance.cluster_Mobs_alpha = cluster_Mobs_alpha; //fiducial: 0.993, flat prior [0.0, 2.0]
  nuisance.cluster_Mobs_beta = cluster_Mobs_beta; //fiducial: 0.0, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma0 = cluster_Mobs_sigma0; //fiducial: 0.456, flat prior [0.0, 1.5]
  nuisance.cluster_Mobs_sigma_qm = cluster_Mobs_sigma_qm; //fiducial: -0.169, flat prior [-1.5, 1.5]
  nuisance.cluster_Mobs_sigma_qz = cluster_Mobs_sigma_qz; //fiducial: 0.0, flat prior [-1.5, 1.5]

  if (nuisance.cluster_Mobs_lgN0 < 0.5 || nuisance.cluster_Mobs_lgN0 > 5.0) return 0;
  if (nuisance.cluster_Mobs_alpha < 0.0 || nuisance.cluster_Mobs_alpha > 2.0) return 0;
  if (nuisance.cluster_Mobs_beta < -1.5 || nuisance.cluster_Mobs_beta > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma0 < 0.0|| nuisance.cluster_Mobs_sigma0 > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma_qm < -1.5 && nuisance.cluster_Mobs_sigma_qm > 1.5) return 0;
  if (nuisance.cluster_Mobs_sigma_qz < -1.5 && nuisance.cluster_Mobs_sigma_qz > 1.5)return 0;

return 1;
}


int set_nuisance_gbias(double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8,double B9, double B10)
{
  int i;
  gbias.b[0] = B1;
  gbias.b[1] = B2;
  gbias.b[2] = B3;
  gbias.b[3] = B4;
  gbias.b[4] = B5;
  gbias.b[5] = B6;
  gbias.b[6] = B7;
  gbias.b[7] = B8;
  gbias.b[8] = B9;
  gbias.b[9] = B10;
  if(like.bias==1){
    for (i = 0; i < 10; i++){
      if (gbias.b[i] < 0.8 || gbias.b[i] > 3.0) return 0;
    }
  }
  return 1;
} 

double log_multi_like(double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5, double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q,double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope)
{
  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0, log_L_GRS=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
    ell_Cluster= create_double_vector(0, Cluster.lbin-1);
    darg=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
    for (l=0;l<Cluster.lbin;l++){
      ell_Cluster[l]=exp(log(Cluster.l_min)+(l+0.5)*darg);
    }
  }
  if (set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu)==0){
    printf("Cosmology out of bounds\n");
    return -1.0e8;
  }
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  if (set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1)==0){
    printf("Shear photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1)==0){
    printf("Clustering photo-z sigma too small\n");
    return -1.0e8;
  }
  if (set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q)==0){
    printf("IA parameters out of bounds\n");
    return -1.0e8; 
  }
  if (set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10)==0){
    printf("Bias out of bounds\n");
    return -1.0e8;
  }
  if (set_nuisance_cluster_Mobs(mass_obs_norm, mass_obs_slope, mass_z_slope, mass_obs_scatter_norm, mass_obs_scatter_mass_slope, mass_obs_scatter_z_slope)==0){
    printf("Mobs out of bounds\n");
    return -1.0e8;
  }
       
  // printf("like %le %le %le %le %le %le %le %le\n",cosmology.Omega_m, cosmology.Omega_v,cosmology.sigma_8,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,cosmology.h0); 
  // printf("like %le %le %le %le\n",gbias.b[0][0], gbias.b[1][0], gbias.b[2][0], gbias.b[3][0]);    
  // for (i=0; i<10; i++){
  //   printf("nuisance %le %le %le\n",nuisance.shear_calibration_m[i],nuisance.bias_zphot_shear[i],nuisance.sigma_zphot_shear[i]);
  // }

  log_L_prior=0.0;
  if(like.Aubourg_Planck_BAO_SN==1) log_L_prior+=log_L_Planck_BAO_SN();
  if(like.SN==1) log_L_prior+=log_L_SN();
  if(like.BAO==1) log_L_prior+=log_L_BAO();
  if(like.Planck==1) log_L_prior+=log_L_Planck();
  if(like.Planck15_BAO_w0wa==1) log_L_prior+=log_L_Planck15_BAO_w0wa();//CH
  if(like.Planck15_BAO_H070p6_JLA_w0wa==1) log_L_prior+=log_L_Planck15_BAO_H070p6_JLA_w0wa();//CH
  if(like.IA!=0) log_L_prior+=log_L_ia();
  if(like.IA!=0) log_L_prior+=log_like_f_red();
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.clphotoz!=0) log_L_prior+=log_L_clphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();
  if(like.clusterMobs==1) log_L_prior+=log_L_clusterMobs();
 
  // printf("%d %d %d %d\n",like.BAO,like.wlphotoz,like.clphotoz,like.shearcalib);
  // printf("logl %le %le %le %le\n",log_L_shear_calib(),log_L_wlphotoz(),log_L_clphotoz(),log_L_clusterMobs());
  int start=0;  
  
  if(like.shear_shear==1) {
    set_data_shear(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    set_data_ggl(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    set_data_clustering(like.Ncl,ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }
  if(like.clusterN==1){ 
    set_data_cluster_N(pred,start);
    start= start+tomo.cluster_Nbin*Cluster.N200_Nbin;
  }
  if(like.clusterWL==1){
    set_data_cgl(ell_Cluster,pred, start);
  }
  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
    for (j=0; j<like.Ndata; j++){
      a=(pred[i]-data_read(1,i))*invcov_read(1,i,j)*(pred[j]-data_read(1,j));
      chisqr=chisqr+a;
    }
   //  if (fabs(data_read(1,i)/pred[i]-1.0) >1.e-4){
   // printf("%d %le %le %le\n",i,data_read(1,i),pred[i],data_read(1,i)/pred[i]);
   //  }
  }
  if (chisqr<0.0){
    printf("errror: chisqr < 0\n");
  }
  if (chisqr<-1.0) exit(EXIT_FAILURE);
  if (like.GRS == 1){
    log_L_GRS = log_like_GRS(cosmology.Omega_m, cosmology.sigma_8, cosmology.n_spec, cosmology.w0,cosmology.wa, cosmology.omb, cosmology.h0,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.,-42.);
  }
//printf("%le\n",chisqr);
  return -0.5*chisqr+log_L_prior+ log_L_GRS;
}

void compute_data_vector(char *details, double OMM, double S8, double NS, double W0,double WA, double OMB, double H0, double MGSigma, double MGmu, double B1, double B2, double B3, double B4,double B5, double B6, double B7, double B8, double B9, double B10, double SP1, double SP2, double SP3, double SP4, double SP5,double SP6, double SP7, double SP8, double SP9, double SP10, double SPS1, double CP1, double CP2, double CP3, double CP4, double CP5, double CP6, double CP7, double CP8, double CP9, double CP10, double CPS1, double M1, double M2, double M3, double M4, double M5, double M6, double M7, double M8, double M9, double M10, double A_ia, double beta_ia, double eta_ia, double eta_ia_highz, double LF_alpha, double LF_P, double LF_Q, double LF_red_alpha, double LF_red_P, double LF_red_Q, double mass_obs_norm, double mass_obs_slope, double mass_z_slope, double mass_obs_scatter_norm, double mass_obs_scatter_mass_slope, double mass_obs_scatter_z_slope)
{

  int i,j,k,m=0,l;
  static double *pred;
  static double *ell;
  static double *ell_Cluster;
  static double darg;
  double chisqr,a,log_L_prior=0.0;
  
  if(ell==0){
    pred= create_double_vector(0, like.Ndata-1);
    ell= create_double_vector(0, like.Ncl-1);
    darg=(log(like.lmax)-log(like.lmin))/like.Ncl;
    for (l=0;l<like.Ncl;l++){
      ell[l]=exp(log(like.lmin)+(l+0.5)*darg);
    }
    ell_Cluster= create_double_vector(0, Cluster.lbin-1);
    darg=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
    for (l=0;l<Cluster.lbin;l++){
      ell_Cluster[l]=exp(log(Cluster.l_min)+(l+0.5)*darg);    
    }
  }
// for (l=0;l<like.Ncl;l++){
//   printf("%d %le\n",i,ell[l]);
// }
  set_cosmology_params(OMM,S8,NS,W0,WA,OMB,H0,MGSigma,MGmu);
  set_nuisance_shear_calib(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10);
  set_nuisance_shear_photoz(SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,SPS1);
  set_nuisance_clustering_photoz(CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,CPS1);
  set_nuisance_ia(A_ia,beta_ia,eta_ia,eta_ia_highz,LF_alpha,LF_P,LF_Q,LF_red_alpha,LF_red_P,LF_red_Q);
  set_nuisance_gbias(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10);
  set_nuisance_cluster_Mobs(mass_obs_norm,mass_obs_slope,mass_z_slope,mass_obs_scatter_norm,mass_obs_scatter_mass_slope,mass_obs_scatter_z_slope);
  
  int start=0;  
  if(like.shear_shear==1) {
    set_data_shear(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.shear_Npowerspectra;
  }
  if(like.shear_pos==1){
    //printf("ggl\n");
    set_data_ggl(like.Ncl, ell, pred, start);
    start=start+like.Ncl*tomo.ggl_Npowerspectra;
  } 
  if(like.pos_pos==1){
    //printf("clustering\n");
    set_data_clustering(like.Ncl,ell,pred, start);
    start=start+like.Ncl*tomo.clustering_Npowerspectra;
  }

  if(like.clusterN==1){ 
    set_data_cluster_N(pred,start);
    start= start+tomo.cluster_Nbin*Cluster.N200_Nbin;
  }
  if(like.clusterWL==1){
    set_data_cgl(ell_Cluster,pred, start);
  }
  FILE *F;
  char filename[300];
  if (strstr(details,"FM") != NULL){
    sprintf(filename,"%s",details);
  }
  else {sprintf(filename,"datav/%s_%s_%s",survey.name,like.probes,details);}
  F=fopen(filename,"w");
  for (i=0;i<like.Ndata; i++){  
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);
}
void write_vector_wrapper(char *details, input_cosmo_params ic, input_nuisance_params in)
{
  compute_data_vector(details, ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,
    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,
    in.lf[0], in.lf[1], in.lf[2], in.lf[3], in.lf[4], in.lf[5],
    in.m_lambda[0], in.m_lambda[1], in.m_lambda[2], in.m_lambda[3],
    in.m_lambda[4], in.m_lambda[5]);
}

double log_like_wrapper(input_cosmo_params ic, input_nuisance_params in)
{
  double like = log_multi_like(ic.omega_m, ic.sigma_8, ic.n_s, ic.w0, ic.wa, ic.omega_b, ic.h0, ic.MGSigma, ic.MGmu,
    in.bias[0], in.bias[1], in.bias[2], in.bias[3],in.bias[4], in.bias[5], in.bias[6], in.bias[7],in.bias[8], in.bias[9], 
    in.source_z_bias[0], in.source_z_bias[1], in.source_z_bias[2], in.source_z_bias[3], in.source_z_bias[4], 
    in.source_z_bias[5], in.source_z_bias[6], in.source_z_bias[7], in.source_z_bias[8], in.source_z_bias[9], 
    in.source_z_s, 
    in.lens_z_bias[0], in.lens_z_bias[1], in.lens_z_bias[2], in.lens_z_bias[3], in.lens_z_bias[4], 
    in.lens_z_bias[5], in.lens_z_bias[6], in.lens_z_bias[7], in.lens_z_bias[8], in.lens_z_bias[9], 
    in.lens_z_s, 
    in.shear_m[0], in.shear_m[1], in.shear_m[2], in.shear_m[3], in.shear_m[4], 
    in.shear_m[5], in.shear_m[6], in.shear_m[7], in.shear_m[8], in.shear_m[9], 
    in.A_ia, in.beta_ia, in.eta_ia, in.eta_ia_highz,
    in.lf[0], in.lf[1], in.lf[2], in.lf[3], in.lf[4], in.lf[5], 
    in.m_lambda[0], in.m_lambda[1], in.m_lambda[2], in.m_lambda[3],
    in.m_lambda[4], in.m_lambda[5]);
  
  return like;
}



void save_zdistr_sources(int zs){
  double z,dz =(redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) for source redshift bin %d\n",zs);
  
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_sources_bin%d.txt",zs);
   F1 = fopen(filename,"w");
   for (z =redshift.shear_zdistrpar_zmin; z< redshift.shear_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, zdistr_photoz(z,zs));
   }
}


void save_zdistr_lenses(int zl){
   double z,dz =(redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/300.0;
  printf("Printing redshift distribution n(z) and bias b(z) for lens redshift bin %d\n",zl);
   
   FILE *F1;
   char filename[300];
   sprintf(filename,"zdistris/zdist_lenses_bin%d.txt", zl);
   F1 = fopen(filename,"w");
   for (z =redshift.shear_zdistrpar_zmin; z< redshift.shear_zdistrpar_zmax; z+= dz){
      fprintf(F1,"%e %e\n", z, pf_photoz(z,zl));
   }
}


 int main(void)
{
  clock_t begin, end;
  double time_spent;
  int i;
/* here, do your time-consuming job */

  init_cosmo();
  init_binning_fourier(15,20.0,15000.0,3000.0,10.0,10);
  init_survey("WFIRST");
  init_galaxies("zdistris/zdistri_WFIRST_LSST_lensing_fine_bin","zdistris/zdistri_WFIRST_LSST_clustering_fine_bin", "gaussian", "gaussian", "source");
  init_clusters();
  init_IA("none", "GAMA");
  init_priors("none","none","PhotoBAO","none");
  //init_priors("none","none","none","Planck15_BAO_w0wa"); //CH: for running with Planck15_BAO_w0wa prior
  init_probes("all_2pt_clusterN_clusterWL");
  // init_probes("all_2pt");
  init_Pdelta("emu",0.8,0.35);
  // init_Pdelta("linear",0.8,0.35);

  for (i =0; i< 10; i++){
    save_zdistr_sources(i);
    save_zdistr_lenses(i);
  }

  // compute_data_vector("fid",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.207,0.993,0.0,0.456,0.0,0.0);
  // compute_data_vector("mu1_Sigma0",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,0.,1.,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  // compute_data_vector("mu1_Sigma1",0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.,1.,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);

  //init_data_inv("cov/cov_LSST_2.600000e+01_1.800000e+04_Rmin10_Ncl25_Ntomo10_2pt_clusterN_clusterWL_inv","datav/LSST_all_2pt_clusterN_clusterWL_fid");
  

  // begin = clock();
  //log_multi_like(0.3156 ,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.35,1.5,1.65,1.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  // printf("knonlin %le\n",nonlinear_scale_computation(1.0));
  // printf("knonlin %le\n",nonlinear_scale_computation(0.5));
  // end = clock();
  // time_spent = (double)(end - begin) / CLOCKS_PER_SEC;      
  // printf("timespent %le\n",time_spent);
  
  //CH BEGINS
  //for testing Planck15_BAO_w0wa prior alone
  //compute_data_vector("fid",3.50989e-01,8.04675e-01,9.64061e-01,-5.05518e-01,-1.46884e+00,5.46245e-02,6.39839e-01,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  //expect 0.0 for the following
  //log_multi_like(3.50989e-01,8.04675e-01.1,9.64061e-01,-5.05518e-01,-1.46884e+00,5.46245e-02,6.39839e-01,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9);
  //expect -13.195605 for the following
  //log_multi_like(0.35449914, 0.81272201,0.97370111, -0.51057289,-1.48353327,0.05517077,0.64623714,0.,0.,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,5.92,1.1,-0.47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.72+log(1.e+14*0.7),1.08,0.0,0.25,0.9,0.9,0.9,0.9); 
  //CH ENDS
  
  return 0;
}



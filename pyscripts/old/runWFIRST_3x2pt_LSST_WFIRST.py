#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_forecasts')

from cosmolike_libs import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin")
data_file = os.path.join(dirname, "datav/WFIRST_all_2pt_fid")
cov_file = os.path.join(dirname, "cov/WFIRST_LSST_3x2pt_inv")
chain_file = os.path.join(dirname, "like/like_WFIRST_LSST_3x2pt")

initcosmo("halofit")
initbins(25,30.0,15000.0,4000.0,21.0,10,10)
initpriors("photo_opti","shear_opti","none","none")
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","none")
# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("all_2pt")
initdatainv(cov_file ,data_file)

sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,1000,560,1,chain_file, blind=False, pool=MPIPool())


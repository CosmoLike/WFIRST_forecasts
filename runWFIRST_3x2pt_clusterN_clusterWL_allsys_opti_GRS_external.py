#!/usr/bin/env python

import sys
# sys.path.append('/home/u17/timeifler/WFIRST_forecasts')
# sys.path.append('/Users/timeifler/WFIRST_forecasts')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin")
data_file = os.path.join(dirname, "datav/WFIRST_3x2pt_clusterN_clusterWL_fid_opti")
cov_file = os.path.join(dirname, "cov/WFIRST_3x2pt_clusterN_clusterWL_inv")
chain_file = "/extra/timeifler/WFIRST_forecasts/chains/like_WFIRST_3x2pt_clusterN_clusterWL_sys_opti_GRS_ext"

initcosmo("halofit")
initbins(25,30.0,15000.0,4000.0,21.0,10,10)
initpriors("photo_opti","shear_opti","Photo_BAO","Planck15_BAO_H070p6_JLA_w0wa")
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","SN10")
initclusters()
initia("none","none")

initGRS=lib.init_GRS
initGRS.argtypes=[ctypes.c_int,ctypes.c_int]
initGRS(2, 1)

# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("3x2pt_clusterN_clusterWL")
initdatainv(cov_file ,data_file)

#sample_params=sample_LCDM_only()
#sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_clusterN_clusterWL_nuisance(get_N_tomo_shear()) 

sample_main(sample_params,2000,560,1,chain_file, blind=False, pool=MPIPool())

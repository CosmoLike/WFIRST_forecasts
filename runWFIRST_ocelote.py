	#!/usr/bin/env /home/u17/timeifler/bin/python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_forecasts')

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering")
data_file = os.path.join(dirname, "datav/WFIRST_all_2pt_clusterN_clusterWL_fid")
cov_file = os.path.join(dirname, "cov/WFIRST_3x2pt_inv")
chain_file = os.path.join(dirname, "like/like_WFIRST_ocelote")

initcosmo()
initbins(25,30.0,15000.0,4000.0,21.0,10)
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","source")
initclusters()
initia("none","none")
initpriors("none","none","none","none") 
# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("all_2pt_clusterN_clusterWL")
initdatainv(cov_file ,data_file)

sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 

sample_main(sample_params,5,14,1,chain_file, blind=False)


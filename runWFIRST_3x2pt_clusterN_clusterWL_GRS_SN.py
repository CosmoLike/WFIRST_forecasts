#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_forecasts')

from cosmolike_libs_opti import * 
from schwimmbad import MPIPool

file_source_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_lensing_fine_bin")
file_lens_z = os.path.join(dirname, "zdistris/zdistri_WFIRST_LSST_clustering_fine_bin")
data_file = os.path.join(dirname, "datav/WFIRST_3x2pt_opti")
cov_file = os.path.join(dirname, "cov/WFIRST_3x2pt_inv")
#chain_file = "/extra/timeifler/WFIRST_forecasts/chains/like_WFIRST_3x2pt_sys_opti"
chain_file = "like/like_3x2pt_clusterN_clusterWL_GRS_SN"
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
initprobes("3x2pt")
initdatainv(cov_file ,data_file)

sample_params=sample_cosmology_only()

sample_main(sample_params,5000,400,1,chain_file, blind=False, pool=MPIPool())


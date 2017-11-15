#!/usr/bin/env /home/heinrich/bin/python

import sys
sys.path.append('/home/heinrich/Projects/WFIRST_forecasts/git/WFIRST_forecasts')

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
file_lens_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
data_file = os.path.join(dirname, "datav/WFIRST_shear_shear_fid_All_0")
cov_file = os.path.join(dirname, "cov/cov_WL_ALL_0cut_4.428000e+01_2.000000e+03_WFIRST_Ncl20_Ntomo10_shear_inv")
chain_file = os.path.join(dirname, "./like/like_WFIRST_ALL_0cut_4.428000e+01_2.000000e+03_Rmin10_Ncl20_Ntomo10_no_sys")

initcosmo()
initbins(20,20.0,5000.0,5000.0,10.0,7,10)
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","redmagic")
initclusters()
initia("none","DEEP2")
initpriors("none","none","none","Planck15_BAO_w0wa") 
# test also with
#initpriors("none","none","none","Planck")
#initpriors("none","none","none","random")
initprobes("shear_shear")
initdatainv(cov_file ,data_file)

sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,4000,64,32,chain_file, blind=False)

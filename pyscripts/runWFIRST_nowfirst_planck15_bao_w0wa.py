#!/usr/bin/env /home/heinrich/anaconda2/bin/python

import sys
sys.path.append('/home/heinrich/Projects/WFIRST_forecasts/git/WFIRST_forecasts')

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
file_lens_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
data_file = os.path.join(dirname, "datav/WFIRST_all_2pt_clusterN_clusterWL_fid")
cov_file = os.path.join("/home/teifler/CosmoLike/WFIRST_forecasts/", "cov/cov_3x2pt_cluster_4.500000e+01_2.000000e+03_WFIRST_2pt_inv")
chain_file = os.path.join(dirname, "./like/like_no_wfirst_planck15_BAO_w0wa_20000steps_512walkers")

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

sample_main(sample_params,10000,512,32,chain_file, blind=False)

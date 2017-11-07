#!/usr/bin/env /home/teifler/anaconda2/bin/python 

import sys
sys.path.append('/home/teifler/CosmoLike/WFIRST_forecasts/')

from cosmolike_libs import * 

file_source_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
file_lens_z = os.path.join(dirname, "zdistris/redshifts_All_0.txt")
data_file = os.path.join(dirname, "datav/WFIRST_all_2pt_fid")
cov_file = os.path.join(dirname, "cov/cov_3x2pt_5.400000e+01_2.000000e+03_WFIRST_Ncl15_Ntomo10_lmax5000_lmin_20_2pt_inv")
chain_file = os.path.join(dirname, "like/like_WFIRST_5.400000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys")

initcosmo()
initbins(15,20.0,5000.0,5000.0,10.0,7,10)
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","source")
initclusters()
initia("none","DEEP2")
initpriors("none","none","PhotoBAO","none")
initprobes("all_2pt")
initdatainv(cov_file ,data_file)

sample_params= sample_cosmology_only(MG=False)
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,10000,64,32,chain_file, blind=False)

#!/usr/bin/env /home/teifler/Python-2.7.8/python

import sys
sys.path.append('/home/teifler/Dropbox/cosmolike/top-level/WFIRST/')
from cosmolike_libs import * 

############ parameters ##################
NTRADE = 1
file_source_z = os.path.join(dirname, "../../zdistris/zdistribution_WFIRST")
file_lens_z = os.path.join(dirname, "../../zdistris/zdistribution_const_comoving")
data_file = os.path.join(dirname, "datav_mg_musigma/WFIRST_all_2pt_fid")
cov_file = os.path.join(dirname, "cov/cov_Multi_Probe_WFIRST_4.500000e+01_2.200000e+03_Rmin10_Ncl20_Ntomo10_2pt_inv")
chain_file = os.path.join(dirname, "./like/like_WFIRST_all_2pt+GRS_trade%d"%(NTRADE))

############standard code ###############
initcosmo()
init_GRS(NTRADE)
initbins(20,20.0,5000.0,5000.0,10.0,7)
initsurvey("WFIRST")
initgalaxies(file_source_z,file_lens_z,"gaussian","gaussian","redmagic")
initpriors("none","none","PhotoBAO","none")
initprobes("all_2pt")
initdatainv(cov_file ,data_file)
sample_params= sample_cosmology_only()
#sample_params = sample_cosmology_shear_nuisance(get_N_tomo_shear())
#sample_params = sample_cosmology_2pt_nuisance(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_nuisance_IA_marg(get_N_tomo_shear(),get_N_tomo_clustering())
#sample_params = sample_cosmology_2pt_cluster_nuisance(get_N_tomo_shear(),get_N_tomo_clustering()) 
sample_main(sample_params,2000,64,32,chain_file, blind=False)

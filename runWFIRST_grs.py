#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_forecasts')

from cosmolike_libs_grs import * 
from schwimmbad import MPIPool


chain_file = "/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/new/GRS_nuisance_kmax025" 


init_GRS()
#sample_params = sample_cosmology_grs()
sample_params = sample_cosmology_grs_nuisance()

sample_main(sample_params,5000,560,1,chain_file, blind=False, pool=MPIPool())


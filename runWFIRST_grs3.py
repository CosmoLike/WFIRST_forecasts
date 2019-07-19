#!/usr/bin/env python

import sys
sys.path.append('/home/u17/timeifler/WFIRST_forecasts')

from cosmolike_libs_grs import * 
from schwimmbad import MPIPool

n_trade = 3

chain_file = "/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade%d_kmin0001_std_corrected_Planck18_BAO_nuisance" %(n_trade)


init_GRS(n_trade, 0)
#sample_params = sample_cosmology_grs()
sample_params = sample_cosmology_grs_nuisance()

sample_main(sample_params,2000,560,1,chain_file, blind=False, pool=MPIPool())


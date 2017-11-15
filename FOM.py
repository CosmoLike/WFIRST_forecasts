import numpy as np
from chainconsumer import ChainConsumer

filename=["/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_4.500000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_3.300000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG","/Users/teifler/Dropbox/cosmolike_store/WFIRST_forecasts/like/like_WFIRST_5.400000e+01_2.000000e+03_Rmin10_Ncl15_Ntomo10_no_sys_MG"]

d1 = np.genfromtxt(filename[0])
d2 = np.genfromtxt(filename[1])
d3 = np.genfromtxt(filename[2])

cov =np.cov(d1)
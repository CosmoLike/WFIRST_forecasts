#!/usr/bin/env python

import numpy as np
import emcee
import ctypes
n_trade = 4
chainfile = "/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade%d_kmin0001_nss8only" %(n_trade)
print chainfile

lib=ctypes.cdll.LoadLibrary("/Users/timeifler/CosmoLike/WFIRST_forecasts/like_grs.so")

init=lib.init_GRS
init.argtypes=[ctypes.c_int,ctypes.c_int]

loglike=lib.log_like_GRS
loglike.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
loglike.restype=ctypes.c_double

init(n_trade, 0)

def likelihood(p):
	like=loglike(0.3156,p[0],p[1],-1.,0.,0.0491685,0.6727,1.19,1.3,1.44,1.56,1.69,1.80,1.91,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.24)
	# like=loglike(p[0],p[1],p[2],p[3],p[4],p[5],p[6],1.19,1.3,1.44,1.56,1.69,1.80,1.91,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.24)
	if like < -1.0e+10:
		return -np.inf
	return like
  
ndim=2
nwalker=200
sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood,threads=32)
# starting_point=np.array([0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727])
starting_point=np.array([0.831,0.9645])
std=np.array([0.05,0.02])
p0 = emcee.utils.sample_ball(starting_point, std, size=nwalker)
f=open(chainfile,'w')
#define varied parameters here
# varied_parameters = ['omega_m']
varied_parameters = ['sigma_8']
# varied_parameters.append('sigma_8')
varied_parameters.append('n_s')
# varied_parameters.append('w0')
# varied_parameters.append('wa')
# varied_parameters.append('omega_b')
# varied_parameters.append('h0')

#write header here
f.write('# ' + '    '.join(varied_parameters)+" log_like\n")

for (p, loglike, state) in sampler.sample(p0,iterations=2000):
    for row,logl in zip(p,loglike):
        p_text = '  '.join(str(r) for r in row)
        f.write('%s %e\n' % (p_text,logl))
    f.flush()
f.close()
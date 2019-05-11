#!/usr/bin/env python

import numpy as np
import emcee
import ctypes
n_trade = 3
chainfile = "/Users/timeifler/Dropbox/cosmolike_store/WFIRST_forecasts/chains/trade%d_sys" %(n_trade)
print chainfile

lib=ctypes.cdll.LoadLibrary("/Users/timeifler/CosmoLike/WFIRST_forecasts/like_grs.so")

init=lib.init_GRS
init.argtypes=[ctypes.c_int,ctypes.c_int]

loglike=lib.log_like_GRS
loglike.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
loglike.restype=ctypes.c_double

init(n_trade, 1)

def likelihood(p):
	like=loglike(p[0],p[1],p[2],p[3],p[4],p[5],p[6],1.19,1.3,1.44,1.56,1.68,1.79,1.9,290.,290.,290.,290.,290.,290.,290.,0.001,0.0,0.24)
	if like < -1.0e+10:
		return -np.inf
	return like
  
ndim=23
nwalker=500
sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood,threads=32)
starting_point=np.array([0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.2,1.3,1.45,1.6,1.70,1.8,1.9,290.,290.,290.,290.,290.,290.,290.,0.001,0.24])
std=np.array([0.05,0.05,0.02,0.1,0.2,0.002,0.03,0.2,0.2,0.2,0.2,0.2,0.2,0.2,30.,30.,30.,30.,30.,30.,30.,0.001,0.01])
p0 = emcee.utils.sample_ball(starting_point, std, size=nwalker)
f=open(chainfile,'w')
for (p, loglike, state) in sampler.sample(p0,iterations=10000):
  #print p
  for row in p:
    #print row
    f.write('%s\n' % ('  '.join(list([str(r) for r in row]))))

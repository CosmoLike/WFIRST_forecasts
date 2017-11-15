#!/usr/bin/env python

import numpy as np
import emcee
import ctypes
n_trade = 3
chainfile = "trade%d_sys" %(n_trade)
print chainfile
lib=ctypes.cdll.LoadLibrary("./like_grs.so")

init=lib.init
init.argtypes=[ctypes.c_int]

loglike=lib.log_multi_like
loglike.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]
loglike.restype=ctypes.c_double

init(n_trade)

def likelihood(p):
#   print p
   like=loglike(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],p[17],p[18],p[19],p[20],p[21],0.0,p[22])
#   print "like = ", like
   return like
  
ndim=23
nwalker=50
sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood,threads=4)
starting_point=np.array([0.3156,0.831,0.9645,-1.,0.,0.0491685,0.6727,1.2,1.3,1.45,1.6,1.70,1.8,1.9,290.,290.,290.,290.,290.,290.,290.,0.001,0.094])
std=np.array([0.05,0.05,0.03,0.5,1,0.005,0.05,0.1,0.1,0.1,0.1,0.1,0.1,0.1,10.,10.,10.,10.,10.,10.,10.,0.001,0.01])/10.
p0 = emcee.utils.sample_ball(starting_point, std, size=nwalker)
f=open(chainfile,'w')
for (p, loglike, state) in sampler.sample(p0,iterations=5000):
  #print p
  for row in p:
    #print row
    f.write('%s\n' % ('  '.join(list([str(r) for r in row]))))
  
    
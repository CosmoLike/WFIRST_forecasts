#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm#, SymLogNorm
# import matplotlib.image as mpimg
import sys
import math
from numpy import linalg as LA

###########################################################

pathOut = "./"
Nell = 25
Ng = 10
Ns = 10

###########################################################
###########################################################
# full dn/dz

path = "zdistris/zdistribution_WFIRST_LSST"
data = np.genfromtxt(path)
zmin = data[:,0]
z = data[:,1]
zmax = data[:,2]
nz = data[:,3]

# nzsum=np.cumsum(nz)
# nz=nz/nzsum
fig=plt.figure(1)
ax=fig.add_subplot(111)

ax.plot(z, nz*300/4*45.0, 'blue', lw=2, label=r'1 tomo bin, true redshift distribution')

# ###########################################################
# ###########################################################
# # full dn/dz no ifc

# path2 = "zdistris/zdistribution_WFIRST_nonifc_0"
# data2 = np.genfromtxt(path2)
# zmin2 = data2[:,0]
# z2 = data2[:,1]
# zmax2 = data2[:,2]
# nz2 = data2[:,3]

# ax.plot(z2, nz2*35.43, 'red', lw=2, label=r'1 tomo bin, true redshift, no IFC')
# ###########################################################
# # source bins

path = "zdistris/zdist_sources_bin"

for iS in range(Ns):
   data = np.genfromtxt(path+str(iS)+".txt")
   z = data[:,0]
   #nz_norm = data[:,1]
   nz = data[:,1]
   #print "normalized to", np.trapz(nz,z)
   
   ax.fill_between(z, 0., nz*4.5, facecolor=plt.cm.winter(1.*iS/Ns), edgecolor='', alpha=0.8)
ax.plot([], [], color=plt.cm.winter(0.), lw=2, alpha=0.7, label=r'10 tomo bins, Gaussian photo-z')

# ###########################################################

# path = pathOut+"datav/LSST_cmbs4_LSSxCMB_Nell25_Ns10_Ng4_zdist_lenses_bin"

# for ig in range(Ng):
#    data = np.genfromtxt(path+str(ig)+".txt")
#    z = data[:,0]
#    nz_norm = data[:,1]
#    nz = data[:,2]
#    bz = data[:,3]
#    print "normalized to", np.trapz(nz,z)

#    ax.fill_between(z, 0., 10.*nz, facecolor=plt.cm.autumn(1.*ig/Ng), edgecolor='', alpha=1)
# ax.plot([], [], color=plt.cm.autumn(0.), lw=2, alpha=0.7, label=r'lenses$\times$10')


# # lens bins
# path2 = "zdistris/zdist_sources_noifc_bin"
# for ig in range(Ng):
#    data = np.genfromtxt(path2+str(ig)+".txt")
#    z = data[:,0]
#    #nz_norm = data[:,1]
#    nz = data[:,1]
#    #print "normalized to", np.trapz(nz,z)

#    ax.fill_between(z, 0., nz*3.543, facecolor=plt.cm.autumn(1.*ig/Ng), edgecolor='', alpha=0.7)
# ax.plot([], [], color=plt.cm.autumn(0.), lw=2, alpha=0.7, label=r'10 tomo bins, Gaussian photo-z, without IFC')

# ###########################################################
# # lensing efficiency

# path = pathOut+"zdist/lensing_efficiency.txt"

# data = np.genfromtxt(path)
# z_cmblens = data[:,0]
# w_cmblens = data[:,1]

# ax.plot(z_cmblens, w_cmblens * 35., 'm-', lw=2, label=r'CMB lensing efficiency')

###########################################################

ax.plot([], [], 'k', lw=2)
#
ax.set_xlim((0., 4))
ax.set_ylim((0., 40.))
ax.legend(loc=1)
ax.set_xlabel(r'$z$', fontsize=24)
ax.set_ylabel(r'$\frac{dN}{d\Omega \, dz}$ [arcmin$^{-2}$]', fontsize=22)

fig.savefig(pathOut+"zbins_LSST.pdf", bbox_inches='tight')

